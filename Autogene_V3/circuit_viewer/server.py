"""
server.py — lightweight Flask API wrapping EcoCyc graph logic from Autogene.py
Run:  python circuit_viewer/server.py [--pgdb-dir 29.6] [--port 5000]
"""

import argparse
import os
import sys
import threading
import time

from flask import Flask, request, jsonify
from flask_cors import CORS

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, ".."))


try:
    from Autogene import (
        EcoCycClient,
        GENE_NODE_PREFIX,
        UNIPROT_NODE_PREFIX,
        resolve_gene_input,
        build_graph,
        get_node_label,
        _node_kind,
        seed_graph_node_id,
        uniprot_url_for_graph_node,
        uniprot_entry_url,
        ecocyc_node_uniprot_entry_url,
        normalize_legacy_uniprot_gene_search_url,
        analyze_sequences_for_graph_hints,
    )
except ImportError as e:
    print(f"ERROR: Could not import Autogene.py - {e}")
    sys.exit(1)

app = Flask(__name__, static_folder=".", static_url_path="")
CORS(app)


@app.after_request
def _no_store_viewer_assets(response):
    """Avoid stale index.html / JSON in the browser while iterating on link logic."""
    p = request.path or ""
    if p in ("/", "/index.html") or p.startswith("/api/"):
        response.headers["Cache-Control"] = "no-store, max-age=0"
    return response

DEFAULT_PGDB = os.path.join(SCRIPT_DIR, "..", "29.6")
_client = None

GRAPH_TIMEOUT = int(os.environ.get("GRAPH_TIMEOUT", "120"))  # UniProt-heavy seeds need more time
# Web BLAST can take several minutes (many ORFs * delay); override with SEQUENCE_ANALYSIS_TIMEOUT
SEQUENCE_ANALYSIS_TIMEOUT = int(os.environ.get("SEQUENCE_ANALYSIS_TIMEOUT", "600"))


def get_client(pgdb_dir):
    global _client
    if _client is None:
        _client = EcoCycClient(pgdb_dir)
    return _client


@app.route("/")
def index():
    return app.send_static_file("index.html")


@app.route("/api/health")
def health():
    return jsonify({"status": "ok"})


@app.route("/api/analyze-sequence", methods=["POST"])
def analyze_sequence():
    """
    Accept raw DNA / FASTA in JSON (`sequence` or `fasta`) or as file upload (`file`).
    ORFs + local blastp when installed, else NCBI web BLASTp; resolves seeds via EcoCyc/UniProt.
    """
    text = ""
    if request.files and request.files.get("file"):
        uf = request.files["file"]
        text = uf.read().decode("utf-8", errors="replace")
    else:
        body = request.get_json(force=True, silent=True) or {}
        text = body.get("sequence") or body.get("fasta") or body.get("text") or ""

    text = (text or "").strip()
    if not text:
        return jsonify(
            {"error": "Send JSON {sequence: ...} or multipart form field `file` (.fa / .fasta)."}
        ), 400

    max_orfs = 10
    body = request.get_json(silent=True) or {}
    raw_max = body.get("max_orfs")
    if raw_max is None and request.form:
        raw_max = request.form.get("max_orfs")
    if raw_max is not None:
        try:
            max_orfs = int(raw_max)
        except (TypeError, ValueError):
            max_orfs = 10
    max_orfs = max(1, min(max_orfs, 20))

    min_id = None
    raw_min = body.get("min_blast_identity")
    if raw_min is None and request.form:
        raw_min = request.form.get("min_blast_identity")
    if raw_min is not None and str(raw_min).strip() != "":
        try:
            min_id = float(raw_min)
        except (TypeError, ValueError):
            min_id = None
    if min_id is not None:
        min_id = max(0.0, min(100.0, min_id))

    out = {}
    exc_holder = []

    def worker():
        try:
            out["result"] = analyze_sequences_for_graph_hints(
                text,
                max_orfs_to_blast=max_orfs,
                client=get_client(app.config["PGDB_DIR"]),
                min_blast_identity_pct=min_id,
            )
        except Exception as e:
            exc_holder.append(e)

    th = threading.Thread(target=worker, daemon=True)
    th.start()
    th.join(timeout=SEQUENCE_ANALYSIS_TIMEOUT)

    if th.is_alive():
        return jsonify(
            {
                "error": f"Timed out after {SEQUENCE_ANALYSIS_TIMEOUT}s (try fewer ORFs, local blastp, or increase SEQUENCE_ANALYSIS_TIMEOUT).",
                "suggested_genes": [],
                "orfs": [],
                "notes": [],
                "resolved_seeds": [],
                "resolution_errors": [],
                "genes_ready_for_graph": [],
            }
        ), 504

    if exc_holder:
        return jsonify({"error": str(exc_holder[0])}), 500

    data = out.get("result") or {}
    return jsonify(data)


@app.route("/api/graph", methods=["POST"])
def graph():
    body = request.get_json(force=True, silent=True) or {}
    gene_tokens = body.get("genes", [])
    include_stringdb = bool(body.get("include_stringdb", False))
    raw_stringdb_min_score = body.get("stringdb_min_score", 700)
    raw_stringdb_species = body.get("stringdb_species", 511145)
    try:
        stringdb_min_score = max(0, min(1000, int(raw_stringdb_min_score)))
    except (TypeError, ValueError):
        stringdb_min_score = 700
    try:
        stringdb_species = int(raw_stringdb_species)
    except (TypeError, ValueError):
        stringdb_species = 511145

    if not gene_tokens or not isinstance(gene_tokens, list):
        return jsonify({"error": "Provide a non-empty 'genes' list."}), 400

    client = get_client(app.config["PGDB_DIR"])

    # ── 1. Resolve gene names (EcoCyc local or UniProt fallback) ───────────
    resolved_genes = []
    resolved_ids = []
    errors = []
    for tok in gene_tokens:
        t = str(tok).strip()
        print(f"[resolve] {t!r} ...", flush=True)
        try:
            r = resolve_gene_input(client, t)
            resolved_genes.append(r)
            resolved_ids.append(seed_graph_node_id(r))
            print(f"[resolve] {t!r} -> {r}", flush=True)
        except (ValueError, RuntimeError) as exc:
            errors.append({"token": t, "error": str(exc)})
            print(f"[resolve] {t!r} -> NOT FOUND: {exc}", flush=True)

    if not resolved_genes:
        return jsonify({
            "orgid": "ECOLI", "nodes": [], "edges": [],
            "resolved": [], "errors": errors,
            "warning": "None of the supplied gene names could be resolved (EcoCyc or UniProt).",
        })

    # ── 2. Build graph in a thread with timeout ────────────────────────────
    result = {}
    exc_holder = []

    def worker():
        try:
            print(f"[graph] building for {resolved_genes} ...", flush=True)
            t0 = time.time()
            G, uniprot_labels = build_graph(
                client,
                resolved_genes,
                include_stringdb=include_stringdb,
                stringdb_min_score=stringdb_min_score,
                stringdb_species=stringdb_species,
            )
            elapsed = time.time() - t0
            print(
                f"[graph] done in {elapsed:.1f}s -- {G.number_of_nodes()} nodes, {G.number_of_edges()} edges",
                flush=True,
            )
            result["G"] = G
            result["uniprot_labels"] = uniprot_labels
        except Exception as e:
            exc_holder.append(e)
            print(f"[graph] ERROR: {e}", flush=True)

    t = threading.Thread(target=worker, daemon=True)
    t.start()
    t.join(timeout=GRAPH_TIMEOUT)

    if t.is_alive():
        print(f"[graph] TIMEOUT after {GRAPH_TIMEOUT}s", flush=True)
        return jsonify({
            "orgid": "ECOLI", "nodes": [], "edges": [],
            "resolved": resolved_ids, "errors": errors,
            "warning": f"Graph build timed out after {GRAPH_TIMEOUT}s. Try fewer genes.",
        })

    if exc_holder:
        return jsonify({"error": str(exc_holder[0])}), 500

    G = result.get("G")
    if G is None:
        return jsonify({"error": "Graph build produced no result."}), 500

    uniprot_labels = result.get("uniprot_labels") or {}

    seed_uniprot_entry_url = {}
    for r in resolved_genes:
        sid = seed_graph_node_id(r)
        acc = getattr(r, "uniprot_accession", None)
        if getattr(r, "source", None) == "uniprot" and sid and acc:
            seed_uniprot_entry_url[sid] = uniprot_entry_url(str(acc).strip())

    # ── 3. Serialise ───────────────────────────────────────────────────────
    name_cache = {}

    def label_for(nid):
        if nid not in name_cache:
            name_cache[nid] = get_node_label(client, nid, uniprot_labels)
        return name_cache[nid]

    nodes = []
    for n in G.nodes():
        row = {"id": n, "label": label_for(n), "kind": _node_kind(n)}
        is_uniprot_style = n.startswith(GENE_NODE_PREFIX) or n.startswith(UNIPROT_NODE_PREFIX)
        if is_uniprot_style:
            row["source"] = "uniprot"
        else:
            row["source"] = "ecocyc"
            row["incoming"] = client.incoming_regulation_summary(n)

        u = seed_uniprot_entry_url.get(n)
        if (
            not u
            and n.startswith(GENE_NODE_PREFIX)
        ):
            u = seed_uniprot_entry_url.get(
                GENE_NODE_PREFIX + n[len(GENE_NODE_PREFIX) :].lower()
            )
        if not u and is_uniprot_style:
            u = uniprot_url_for_graph_node(G, n)
        if not u:
            u = ecocyc_node_uniprot_entry_url(client, n)
        if u:
            row["external_url"] = normalize_legacy_uniprot_gene_search_url(u)
        nodes.append(row)
    edges = []
    for u, v, d in G.edges(data=True):
        row = {"from": u, "to": v, "interaction": d.get("interaction", "")}
        if "source" in d:
            row["source"] = d.get("source")
        if "confidence" in d:
            row["confidence"] = d.get("confidence")
        edges.append(row)

    return jsonify({
        "orgid": "ECOLI",
        "nodes": nodes,
        "edges": edges,
        "resolved": resolved_ids,
        "errors": errors,
        "options": {
            "include_stringdb": include_stringdb,
            "stringdb_min_score": stringdb_min_score,
            "stringdb_species": stringdb_species,
        },
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgdb-dir", default=DEFAULT_PGDB)
    parser.add_argument("--port", type=int, default=5000)
    parser.add_argument("--host", default="127.0.0.1")
    args = parser.parse_args()

    pgdb = os.path.abspath(args.pgdb_dir)
    if not os.path.isdir(pgdb):
        print(f"ERROR: PGDB directory not found: {pgdb}")
        sys.exit(1)

    app.config["PGDB_DIR"] = pgdb
    print(f"Pre-loading EcoCyc PGDB from {pgdb} ...", flush=True)
    get_client(pgdb)
    print(f"Ready - http://{args.host}:{args.port}/", flush=True)
    app.run(host=args.host, port=args.port, debug=False, threaded=True)