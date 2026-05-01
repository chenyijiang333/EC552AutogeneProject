import json
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import urllib.error

from ecocyc_client import ORG_ID, ResolvedGene, get_node_label, parse_gene_ids, resolve_gene_input
from remote_clients import (
    GENE_NODE_PREFIX,
    UNIPROT_NODE_PREFIX,
    STRINGDB_DEFAULT_MIN_SCORE,
    STRINGDB_DEFAULT_SPECIES,
    _fetch_uniprot_entry_json,
    _gene_symbol_from_entry,
    _protein_label_from_entry,
    stringdb_interaction_partners,
    uniprot_interaction_partners,
    uniprot_node_id,
    uniprot_gene_node_id,
)

def uniprot_node_id(accession: str) -> str:
    return f"{UNIPROT_NODE_PREFIX}{accession}"


def uniprot_gene_node_id(gene_symbol: str) -> str:
    """Stable graph id for UniProt-backed gene seeds (case-insensitive symbols share one node)."""
    return f"{GENE_NODE_PREFIX}{gene_symbol.strip().lower()}"


def seed_graph_node_id(r: ResolvedGene) -> str:
    """Primary graph node id for a resolved input (matches build_graph anchors)."""
    if r.source == "ecocyc" and r.ecocyc_id:
        return r.ecocyc_id
    if r.source == "uniprot" and r.uniprot_accession:
        gene_sym = (r.uniprot_gene_symbol or r.uniprot_accession).strip()
        return uniprot_gene_node_id(gene_sym)
    return ""


def resolved_gene_to_preview(token: str, r: ResolvedGene) -> Dict[str, Any]:
    """JSON-friendly summary for API / web (EcoCyc id or UniProt-backed gene seed)."""
    row: Dict[str, Any] = {
        "token": token,
        "source": r.source,
        "graph_seed_id": seed_graph_node_id(r),
    }
    if r.ecocyc_id:
        row["ecocyc_id"] = r.ecocyc_id
    if r.uniprot_accession:
        row["uniprot_accession"] = r.uniprot_accession
    if r.uniprot_gene_symbol:
        row["uniprot_gene_symbol"] = r.uniprot_gene_symbol
    if r.protein_label:
        row["protein_label"] = r.protein_label
    return row


def resolve_blast_tokens_to_graph_seeds(
    client, tokens: List[str]
) -> Tuple[List[Dict[str, Any]], List[Dict[str, str]]]:
    """
    Map BLAST-derived symbols to the same seeds ``build_graph`` uses: local EcoCyc
    gene id when possible, else UniProt (gene + product + PPI path).
    """
    resolved: List[Dict[str, Any]] = []
    errors: List[Dict[str, str]] = []
    for tok in tokens:
        t = str(tok).strip()
        if not t:
            continue
        try:
            r = resolve_gene_input(client, t)
            resolved.append(resolved_gene_to_preview(t, r))
        except (ValueError, RuntimeError) as exc:
            errors.append({"token": t, "error": str(exc)})
    return resolved, errors


def filter_graph_to_paths_between_seeds(G: nx.DiGraph, seed_ids: List[str]) -> nx.DiGraph:
    """
    Keep only seed nodes and vertices that lie on an undirected shortest path between
    some pair of seeds. If fewer than two seeds appear in G, return G unchanged.
    """
    seeds = [s for s in seed_ids if s and s in G]
    if len(seeds) < 2:
        return G

    U = G.to_undirected()
    keep: set = set(seeds)

    for i, s in enumerate(seeds):
        for t in seeds[i + 1 :]:
            dist_s = nx.single_source_shortest_path_length(U, s)
            if t not in dist_s:
                continue
            D = dist_s[t]
            dist_t = nx.single_source_shortest_path_length(U, t)
            for v in G.nodes():
                ds = dist_s.get(v)
                dt = dist_t.get(v)
                if ds is not None and dt is not None and ds + dt == D:
                    keep.add(v)

    return G.subgraph(keep).copy()

def build_graph(
    client,
    resolved: List[ResolvedGene],
    include_stringdb: bool = False,
    stringdb_min_score: int = STRINGDB_DEFAULT_MIN_SCORE,
    stringdb_species: int = STRINGDB_DEFAULT_SPECIES,
):
    G = nx.DiGraph()
    uniprot_labels = {}

    for r in resolved:
        if r.source == "ecocyc":
            gid = r.ecocyc_id
            if not gid:
                continue
            print(f"\nProcessing {gid} (EcoCyc)")

            protein = client.protein_of_gene(gid)
            if protein:
                G.add_edge(gid, protein, interaction="creates")

            query_id = protein if protein else gid

            inhib_xml = client.api("direct-inhibitors", query_id)
            for target in parse_gene_ids(inhib_xml):
                G.add_edge(query_id, target, interaction="inhibits")

            act_xml = client.api("direct-activators", query_id)
            for target in parse_gene_ids(act_xml):
                G.add_edge(query_id, target, interaction="activates")

            if protein:
                has_specific = any(
                    d.get("interaction") in ("inhibits", "activates")
                    for _, _, d in G.out_edges(query_id, data=True)
                )
                if not has_specific:
                    reg_xml = client.api("regulon-of-protein", protein)
                    for target in parse_gene_ids(reg_xml):
                        if not G.has_edge(protein, target):
                            G.add_edge(protein, target, interaction="regulates")
        else:
            acc = r.uniprot_accession
            if not acc:
                continue
            gene_sym = (r.uniprot_gene_symbol or acc).strip()
            gene_nid = uniprot_gene_node_id(gene_sym)
            prot_nid = uniprot_node_id(acc)
            uniprot_labels[gene_nid] = gene_sym
            uniprot_labels[prot_nid] = r.protein_label or acc

            print(
                f"\nProcessing {gene_nid} -> {prot_nid} "
                f"(gene product + UniProt protein-protein interactions)"
            )

            G.add_edge(gene_nid, prot_nid, interaction="creates")

            partners = uniprot_interaction_partners(acc)
            partner_fetch_cache = {}

            def _partner_gene_and_protein_label(pacc: str, hint_gene: Optional[str]):
                if pacc in partner_fetch_cache:
                    return partner_fetch_cache[pacc]
                if hint_gene and hint_gene.strip():
                    try:
                        pe = _fetch_uniprot_entry_json(pacc)
                        plab = _protein_label_from_entry(pe, pacc)
                    except urllib.error.URLError:
                        plab = pacc
                    partner_fetch_cache[pacc] = (hint_gene.strip(), plab)
                    return partner_fetch_cache[pacc]
                try:
                    pe = _fetch_uniprot_entry_json(pacc)
                    gs = _gene_symbol_from_entry(pe)
                    plab = _protein_label_from_entry(pe, pacc)
                    partner_fetch_cache[pacc] = (gs, plab)
                except urllib.error.URLError:
                    partner_fetch_cache[pacc] = (None, pacc)
                return partner_fetch_cache[pacc]

            if not partners:
                print(
                    f"  (no curated binary interactions in UniProt INTERACTION comment for {acc})"
                )

            for partner_acc, partner_gene in partners:
                p_prot = uniprot_node_id(partner_acc)
                g_sym, p_lab = _partner_gene_and_protein_label(partner_acc, partner_gene)
                if p_prot not in uniprot_labels:
                    uniprot_labels[p_prot] = p_lab
                if g_sym:
                    g_n = uniprot_gene_node_id(g_sym)
                    if g_n not in uniprot_labels:
                        uniprot_labels[g_n] = g_sym
                    G.add_edge(g_n, p_prot, interaction="creates")
                G.add_edge(prot_nid, p_prot, interaction="interacts")

    if include_stringdb:
        uniprot_symbol_to_acc: Dict[str, str] = {}
        for r in resolved:
            if r.source != "uniprot":
                continue
            acc = (r.uniprot_accession or "").strip()
            sym = (r.uniprot_gene_symbol or "").strip().lower()
            if acc and sym:
                uniprot_symbol_to_acc[sym] = acc

        def _node_accession(node_id: str) -> Optional[str]:
            if not node_id:
                return None
            if node_id.startswith(UNIPROT_NODE_PREFIX):
                return node_id[len(UNIPROT_NODE_PREFIX) :].strip() or None
            if node_id.startswith(GENE_NODE_PREFIX):
                sym = node_id[len(GENE_NODE_PREFIX) :].strip().lower()
                return uniprot_symbol_to_acc.get(sym)
            return client.uniprot_accession_for_graph_node(node_id)

        def _symbols_from_gene_record(gid: str) -> List[str]:
            rec = client.gene_records.get(gid) or {}
            out: List[str] = []
            for k in ("COMMON-NAME", "SYNONYMS"):
                for s in rec.get(k, []) or []:
                    sym = str(s).strip().lower()
                    if sym:
                        out.append(sym)
            return out

        def _node_symbols(node_id: str) -> List[str]:
            if not node_id:
                return []
            if node_id.startswith(GENE_NODE_PREFIX):
                sym = node_id[len(GENE_NODE_PREFIX) :].strip().lower()
                return [sym] if sym else []
            if node_id in client.gene_records:
                return _symbols_from_gene_record(node_id)
            gid = client.protein_to_gene.get(node_id)
            if gid:
                return _symbols_from_gene_record(gid)
            return []

        partner_score_cache: Dict[str, Dict[str, float]] = {}

        def _partner_scores_for_accession(acc: str) -> Dict[str, float]:
            if acc in partner_score_cache:
                return partner_score_cache[acc]
            score_map: Dict[str, float] = {}
            for partner_gene, conf in stringdb_interaction_partners(
                acc,
                species=stringdb_species,
                min_score=stringdb_min_score,
            ):
                sym = (partner_gene or "").strip().lower()
                if not sym:
                    continue
                cur = score_map.get(sym, 0.0)
                if conf > cur:
                    score_map[sym] = float(conf)
            partner_score_cache[acc] = score_map
            return score_map

        for u, v, d in G.edges(data=True):
            if d.get("interaction") not in ("inhibits", "activates"):
                continue
            u_acc = (_node_accession(u) or "").strip()
            if not u_acc:
                continue
            v_syms = _node_symbols(v)
            if not v_syms:
                continue

            partner_scores = _partner_scores_for_accession(u_acc)
            best_score = 0.0
            for sym in v_syms:
                sc = float(partner_scores.get(sym, 0.0) or 0.0)
                if sc > best_score:
                    best_score = sc
            if best_score <= 0.0:
                continue

            best_score = round(best_score, 3)
            d["stringdb_score"] = best_score
            # Preserve existing viewer compatibility while keeping interaction type unchanged.
            d["confidence"] = best_score

        # Add explicit STRING edges when two existing graph nodes have no direct edge
        # but pass the STRING threshold by accession+symbol matching.
        nodes = list(G.nodes())
        for i, a in enumerate(nodes):
            a_acc = (_node_accession(a) or "").strip()
            if not a_acc:
                continue
            a_syms = _node_symbols(a)
            scores_a = _partner_scores_for_accession(a_acc)
            for b in nodes[i + 1 :]:
                # Respect user expectation: only add STRING edge when no edge already exists.
                if G.has_edge(a, b) or G.has_edge(b, a):
                    continue
                b_syms = _node_symbols(b)
                if not b_syms:
                    continue

                best_score = 0.0
                for sym in b_syms:
                    sc = float(scores_a.get(sym, 0.0) or 0.0)
                    if sc > best_score:
                        best_score = sc

                # Try reverse lookup too (helps when one direction has sparse aliases).
                if best_score <= 0.0:
                    b_acc = (_node_accession(b) or "").strip()
                    if b_acc:
                        scores_b = _partner_scores_for_accession(b_acc)
                        for sym in a_syms:
                            sc = float(scores_b.get(sym, 0.0) or 0.0)
                            if sc > best_score:
                                best_score = sc

                if best_score <= 0.0:
                    continue

                best_score = round(best_score, 3)
                src, dst = (a, b) if str(a) <= str(b) else (b, a)
                G.add_edge(
                    src,
                    dst,
                    interaction="stringdb",
                    source="stringdb",
                    stringdb_score=best_score,
                    confidence=best_score,
                )


    seed_ids = [seed_graph_node_id(r) for r in resolved]
    G = filter_graph_to_paths_between_seeds(G, seed_ids)
    uniprot_labels = {
        k: v for k, v in uniprot_labels.items() if k in G
    }

    return G, uniprot_labels


def _node_kind(frame_id):
    """Rough EcoCyc frame-id classification for the circuit viewer."""
    fid = frame_id or ""
    if fid.startswith(GENE_NODE_PREFIX):
        return "gene"
    if fid.startswith(UNIPROT_NODE_PREFIX):
        return "protein"
    if fid.startswith("EG") and "MONOMER" not in fid:
        return "gene"
    if "MONOMER" in fid:
        return "protein"
    if fid.startswith("G") and "-" in fid:
        return "gene"
    return "other"


def export_graph_json(
    client, graph, path, org_id=ORG_ID, name_cache=None, uniprot_labels=None
):
    """Write nodes/edges + labels for the local circuit viewer (JSON)."""
    cache = name_cache if name_cache is not None else {}

    def label_for(nid):
        if nid in cache:
            return cache[nid]
        name = get_node_label(client, nid, uniprot_labels)
        cache[nid] = name
        return name

    nodes = []
    for nid in graph.nodes():
        nodes.append(
            {
                "id": nid,
                "label": label_for(nid),
                "kind": _node_kind(nid),
            }
        )

    edges = []
    for u, v, data in graph.edges(data=True):
        edges.append(
            {
                "from": u,
                "to": v,
                "interaction": data.get("interaction", ""),
            }
        )

    payload = {"orgid": org_id, "nodes": nodes, "edges": edges}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote circuit JSON for viewer: {path}")
