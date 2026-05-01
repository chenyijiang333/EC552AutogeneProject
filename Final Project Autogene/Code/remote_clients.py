import json
import os
import re
from typing import Any, Dict, List, Optional, Tuple

import urllib.error
import urllib.parse
import urllib.request
import networkx as nx

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb"
UNIPROT_ORGANISM_ID = "83333"
UNIPROT_TAXON_ECOLI = "562"
UNIPROT_NODE_PREFIX = "UniProt:"
GENE_NODE_PREFIX = "Gene:"
STRINGDB_API_BASE = os.environ.get("STRINGDB_API_BASE", "https://version-12-0.string-db.org/api/json")
STRINGDB_DEFAULT_SPECIES = int(os.environ.get("STRINGDB_SPECIES", "511145"))
STRINGDB_DEFAULT_MIN_SCORE = int(os.environ.get("STRINGDB_MIN_SCORE", "700"))
STRINGDB_MAX_PARTNERS = max(1, int(os.environ.get("STRINGDB_MAX_PARTNERS", "40")))

def uniprot_entry_url(accession: str) -> str:
    return f"https://www.uniprot.org/uniprotkb/{urllib.parse.quote(accession.strip())}"


def uniprot_search_url_for_gene_symbol(gene_symbol: str) -> str:
    """
    UniProt website search: **gene name only** + **Swiss-Prot (reviewed)** — no organism filter
    (so symbols like RTEL1 are not incorrectly scoped to *E. coli* K-12).
    """
    sym = gene_symbol.strip()
    q = f"(gene:{sym}) AND (reviewed:true)"
    return f"https://www.uniprot.org/uniprotkb?query={urllib.parse.quote(q)}"


def normalize_legacy_uniprot_gene_search_url(url: Optional[str]) -> Optional[str]:
    """
    Rewrite UniProt **search** links that add organism / taxon filters (including K-12
    ``83333``) to ``(gene:…) AND (reviewed:true)`` only. Skips single-entry ``/uniprotkb/ACCESSION`` URLs.
    """
    if not url or "uniprot.org" not in url.lower():
        return url
    try:
        p = urllib.parse.urlparse(url)
    except Exception:
        return url
    path = (p.path or "").rstrip("/")
    if re.match(r"^/uniprotkb/[^/]+$", path, re.I):
        return url

    dec = urllib.parse.unquote(url.replace("+", " "))
    m = re.search(r"\(gene:([^)]+)\)", dec, re.I)
    if not m:
        return url
    sym = m.group(1).strip()
    low = dec.lower()
    if "organism_id" not in low and "taxonomy_id" not in low and "83333" not in dec:
        return url
    new_q = f"(gene:{sym}) AND (reviewed:true)"
    return f"https://www.uniprot.org/uniprotkb?query={urllib.parse.quote(new_q)}"


def uniprot_url_for_graph_node(graph: nx.DiGraph, node_id: str) -> Optional[str]:
    """Best UniProt URL for a viewer node (protein entry, or gene symbol search)."""
    if node_id.startswith(UNIPROT_NODE_PREFIX):
        return uniprot_entry_url(node_id[len(UNIPROT_NODE_PREFIX) :])
    if node_id.startswith(GENE_NODE_PREFIX):
        for _, tgt, data in graph.out_edges(node_id, data=True):
            if data.get("interaction") == "creates" and tgt.startswith(UNIPROT_NODE_PREFIX):
                return uniprot_entry_url(tgt[len(UNIPROT_NODE_PREFIX) :])
        return uniprot_search_url_for_gene_symbol(node_id[len(GENE_NODE_PREFIX) :])
    return None


def _uniprot_search(query: str, size: int = 25) -> List[dict]:
    params = {"query": query, "format": "json", "size": str(size)}
    url = f"{UNIPROT_SEARCH_URL}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    return payload.get("results") or []


def _protein_label_from_uniprot_hit(hit: dict) -> str:
    """Recommended protein name (product), not gene symbol."""
    pd = hit.get("proteinDescription") or {}
    rn = pd.get("recommendedName") or {}
    fn = rn.get("fullName") or {}
    if fn.get("value"):
        return fn["value"]
    return hit.get("primaryAccession") or ""


def _gene_symbol_from_hit(hit: dict) -> Optional[str]:
    for g in hit.get("genes") or []:
        gn = (g.get("geneName") or {}).get("value")
        if gn:
            return gn.strip()
        for syn in g.get("synonyms") or []:
            sv = (syn.get("value") or "").strip()
            if sv:
                return sv
    return None


def _gene_symbol_from_entry(entry: dict) -> Optional[str]:
    for g in entry.get("genes") or []:
        gn = (g.get("geneName") or {}).get("value")
        if gn:
            return gn.strip()
        for syn in g.get("synonyms") or []:
            sv = (syn.get("value") or "").strip()
            if sv:
                return sv
    return None


def _protein_label_from_entry(entry: dict, accession: str) -> str:
    pd = entry.get("proteinDescription") or {}
    rn = pd.get("recommendedName") or {}
    fn = rn.get("fullName") or {}
    if fn.get("value"):
        return fn["value"]
    return accession


def _looks_like_uniprot_accession_query(token: str) -> bool:
    """True if ``accession:TOKEN`` is valid for UniProt search (avoids HTTP 400 on symbols like ``bla``)."""
    t = token.strip().upper()
    return bool(
        re.match(
            r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:-\d+)?$|^[OPQ][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9](?:-\d+)?$",
            t,
        )
    )


def _uniprot_hit_taxon_id(hit: dict) -> Optional[int]:
    org = hit.get("organism") or {}
    tid = org.get("taxonId")
    if tid is None:
        return None
    try:
        return int(tid)
    except (TypeError, ValueError):
        return None


def _uniprot_hit_gene_tokens(hit: dict) -> List[str]:
    toks: List[str] = []
    for g in hit.get("genes") or []:
        gn = (g.get("geneName") or {}).get("value") or ""
        if gn.strip():
            toks.append(gn.strip().lower())
        for syn in g.get("synonyms") or []:
            sv = (syn.get("value") or "").strip()
            if sv:
                toks.append(sv.lower())
    return toks


def _uniprot_hit_rank_tuple(hit: dict) -> Tuple[int, float]:
    """
    Single numeric score (higher is better) for disambiguating many exact gene-name hits.
    Prefers *E. coli* K-12 / species, then common reference eukaryotes (human > rodent > …),
    then Swiss-Prot and UniProt annotation score.
    """
    tid = _uniprot_hit_taxon_id(hit) or 0
    et = (hit.get("entryType") or "").lower()
    reviewed = 1 if "reviewed" in et else 0
    score = float(hit.get("annotationScore") or 0.0)
    k12 = int(UNIPROT_ORGANISM_ID)
    tax_ecoli = int(UNIPROT_TAXON_ECOLI)
    taxon_tier = 0
    if tid == k12:
        taxon_tier = 1_000_000
    elif tid == tax_ecoli:
        taxon_tier = 500_000
    elif tid == 9606:  # Homo sapiens
        taxon_tier = 200_000
    elif tid in (10090, 10116):  # mouse / rat
        taxon_tier = 120_000
    elif tid in (559292, 4932):  # budding yeast (strain / species)
        taxon_tier = 80_000
    elif tid == 3702:  # Arabidopsis
        taxon_tier = 70_000
    elif tid in (7227, 7955):  # fly / zebrafish
        taxon_tier = 50_000
    return (taxon_tier + reviewed * 10_000 + int(score * 100), score)


def _pick_uniprot_search_hit(results: List[dict], token: str) -> Optional[dict]:
    """Prefer exact gene/synonym match; tie-break toward *E. coli* K-12 / species, then Swiss-Prot."""
    if not results:
        return None
    q = token.strip().lower()
    if not q:
        return results[0]

    exact = [r for r in results if q in _uniprot_hit_gene_tokens(r)]
    pool = exact if exact else list(results)

    acc = token.strip().upper()
    if _looks_like_uniprot_accession_query(acc):
        for r in pool:
            if (r.get("primaryAccession") or "").upper() == acc:
                return r

    if exact:
        return max(exact, key=_uniprot_hit_rank_tuple)

    return max(pool, key=_uniprot_hit_rank_tuple)


def uniprot_search_best_hit(token: str) -> Optional[Tuple[str, str, str]]:
    """
    Return (primaryAccession, gene_symbol, protein_product_name) for token, or None.
    Uses several queries; avoids OR-combines that return HTTP 400 from UniProt.
    """
    t = token.strip()
    if not t:
        return None

    query_attempts = [
        f"gene:{t} AND organism_id:{UNIPROT_ORGANISM_ID}",
        f"gene:{t} AND taxonomy_id:{UNIPROT_TAXON_ECOLI}",
        f"gene:{t} AND reviewed:true",
        f"gene:{t}",
    ]
    if _looks_like_uniprot_accession_query(t):
        query_attempts.insert(1, f"accession:{t}")
    results: List[dict] = []
    for q in query_attempts:
        try:
            results = _uniprot_search(q, size=25)
        except urllib.error.HTTPError:
            continue
        except urllib.error.URLError:
            continue
        if results:
            break
    if not results:
        return None

    hit = _pick_uniprot_search_hit(results, t)
    if not hit:
        return None
    acc = hit.get("primaryAccession")
    if not acc:
        return None
    gene_sym = _gene_symbol_from_hit(hit)
    prot_name = _protein_label_from_uniprot_hit(hit)
    if not gene_sym or not prot_name:
        try:
            entry = _fetch_uniprot_entry_json(acc)
            if not gene_sym:
                gene_sym = _gene_symbol_from_entry(entry) or ""
            if not prot_name or prot_name == acc:
                prot_name = _protein_label_from_entry(entry, acc)
        except urllib.error.URLError:
            pass
    if not gene_sym:
        gene_sym = t
    if not prot_name:
        prot_name = acc
    return acc, gene_sym, prot_name


def _fetch_uniprot_entry_json(accession: str) -> dict:
    url = f"{UNIPROT_ENTRY_URL}/{urllib.parse.quote(accession)}?format=json"
    with urllib.request.urlopen(url, timeout=30) as resp:
        return json.loads(resp.read().decode("utf-8"))


def uniprot_interaction_partners(accession: str) -> List[Tuple[str, Optional[str]]]:
    """
    Protein–protein interactors from UniProt INTERACTION comments (curated binary).
    Returns list of (partner_accession, partner_gene_name_or_None).
    """
    try:
        entry = _fetch_uniprot_entry_json(accession)
    except urllib.error.URLError as e:
        print(f"UniProt entry fetch failed for {accession!r}: {e}")
        return []

    partners: List[Tuple[str, Optional[str]]] = []
    for c in entry.get("comments") or []:
        if c.get("commentType") != "INTERACTION":
            continue
        for it in c.get("interactions") or []:
            one = it.get("interactantOne") or {}
            two = it.get("interactantTwo") or {}
            a1 = one.get("uniProtKBAccession")
            a2 = two.get("uniProtKBAccession")
            if not a1 or not a2:
                continue
            if a1 == accession and a2 != accession:
                partners.append((a2, two.get("geneName")))
            elif a2 == accession and a1 != accession:
                partners.append((a1, one.get("geneName")))
    return partners


def stringdb_interaction_partners(
    accession: str,
    species: int = STRINGDB_DEFAULT_SPECIES,
    min_score: int = STRINGDB_DEFAULT_MIN_SCORE,
) -> List[Tuple[str, float]]:
    """
    STRING functional partners for a UniProt accession.
    Returns list of (partner_gene_symbol, confidence_0_to_1), highest score first.
    """
    def _stringdb_json(method: str, params: Dict[str, str]) -> List[Dict[str, Any]]:
        url = f"{STRINGDB_API_BASE}/{method}?{urllib.parse.urlencode(params)}"
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                return data if isinstance(data, list) else []
        except urllib.error.HTTPError as e:
            # Missing mapping / no partners from STRING often returns 404.
            if e.code != 404:
                print(f"STRING {method} HTTP error for {accession!r}: {e}")
            return []
        except urllib.error.URLError as e:
            print(f"STRING {method} fetch failed for {accession!r}: {e}")
            return []
        except Exception as e:
            print(f"STRING {method} parse failed for {accession!r}: {e}")
            return []

    req_score = max(0, min(1000, int(min_score)))
    id_rows = _stringdb_json(
        "get_string_ids",
        {
            "identifiers": accession.strip(),
            "species": str(int(species)),
            "limit": "1",
            "caller_identity": "autogene-circuit-viewer",
        },
    )
    if not id_rows:
        return []
    sid = str(id_rows[0].get("stringId") or "").strip()
    if not sid:
        return []
    rows = _stringdb_json(
        "interaction_partners",
        {
            "identifiers": sid,
            "species": str(int(species)),
            "required_score": str(req_score),
            "limit": str(STRINGDB_MAX_PARTNERS),
            "caller_identity": "autogene-circuit-viewer",
        },
    )
    if not rows:
        return []

    out: Dict[str, float] = {}
    q_sid = sid.upper()
    for r in rows or []:
        pa = str(r.get("preferredName_A") or "").strip()
        pb = str(r.get("preferredName_B") or "").strip()
        id_a = str(r.get("stringId_A") or "")
        id_b = str(r.get("stringId_B") or "")
        score_raw = r.get("score", 0)
        try:
            score = float(score_raw)
        except (TypeError, ValueError):
            continue
        # STRING API may return score in [0,1] or [0,1000] depending endpoint/version.
        if score > 1.0:
            score = score / 1000.0
        score = max(0.0, min(1.0, score))
        if score * 1000.0 < req_score:
            continue

        choose_b = (
            id_a.upper() == q_sid
            or id_a.upper().endswith("." + q_sid)
        )
        choose_a = (
            id_b.upper() == q_sid
            or id_b.upper().endswith("." + q_sid)
        )
        partner = pb if choose_b else (pa if choose_a else "")
        partner = partner.strip()
        if not partner:
            continue
        key = partner.lower()
        prev = out.get(key)
        if prev is None or score > prev:
            out[key] = score

    ranked = sorted(out.items(), key=lambda kv: kv[1], reverse=True)[:STRINGDB_MAX_PARTNERS]
    return [(sym, score) for sym, score in ranked]


def uniprot_node_id(accession: str) -> str:
    return f"{UNIPROT_NODE_PREFIX}{accession}"


def uniprot_gene_node_id(gene_symbol: str) -> str:
    """Stable graph id for UniProt-backed gene seeds (case-insensitive symbols share one node)."""
    return f"{GENE_NODE_PREFIX}{gene_symbol.strip().lower()}"

