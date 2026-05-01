import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Dict, List, Optional

from remote_clients import GENE_NODE_PREFIX, UNIPROT_NODE_PREFIX, uniprot_entry_url, uniprot_search_best_hit

ORG_ID = "ECOLI"
DEFAULT_PGDB_DIR = "29.6"

@dataclass
class ResolvedGene:
    source: str
    ecocyc_id: Optional[str] = None
    uniprot_accession: Optional[str] = None
    uniprot_gene_symbol: Optional[str] = None
    protein_label: Optional[str] = None

class EcoCycClient:
    def __init__(self, pgdb_dir=DEFAULT_PGDB_DIR):
        self.pgdb_dir = pgdb_dir
        data_dir = os.path.join(pgdb_dir, "data")
        self.genes_dat = os.path.join(data_dir, "genes.dat")
        self.proteins_dat = os.path.join(data_dir, "proteins.dat")
        self.regulation_dat = os.path.join(data_dir, "regulation.dat")
        self.promoters_dat = os.path.join(data_dir, "promoters.dat")
        self.transunits_dat = os.path.join(data_dir, "transunits.dat")

        needed = [
            self.genes_dat,
            self.proteins_dat,
            self.regulation_dat,
            self.promoters_dat,
            self.transunits_dat,
        ]
        missing = [p for p in needed if not os.path.exists(p)]
        if missing:
            raise FileNotFoundError(
                "Missing EcoCyc flat files in local PGDB folder:\n  "
                + "\n  ".join(missing)
            )

        print(f"Using local EcoCyc PGDB: {os.path.abspath(pgdb_dir)}")
        self._load_local_db()

    @staticmethod
    def _iter_dat_records(path):
        record = {}
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                if line == "//":
                    if record:
                        yield record
                        record = {}
                    continue
                if line.startswith("/"):
                    continue
                if " - " not in line:
                    continue
                key, val = line.split(" - ", 1)
                record.setdefault(key, []).append(val.strip())
        if record:
            yield record

    def _load_local_db(self):
        self.gene_records = {}
        self.protein_records = {}
        self.promoter_to_tus = {}
        self.tu_to_genes = {}
        self.regulations_by_regulator = {}
        self.name_index = {}

        self.protein_to_gene = {}
        for rec in self._iter_dat_records(self.genes_dat):
            gid = (rec.get("UNIQUE-ID") or [None])[0]
            if not gid:
                continue
            self.gene_records[gid] = rec
            prod = (rec.get("PRODUCT") or [None])[0]
            if prod:
                self.protein_to_gene[prod] = gid
            for key in ("COMMON-NAME", "SYNONYMS", "ACCESSION-1", "ACCESSION-2"):
                for token in rec.get(key, []):
                    self.name_index.setdefault(token.lower(), set()).add(gid)

        for rec in self._iter_dat_records(self.proteins_dat):
            pid = (rec.get("UNIQUE-ID") or [None])[0]
            if not pid:
                continue
            self.protein_records[pid] = rec

        for rec in self._iter_dat_records(self.transunits_dat):
            tu = (rec.get("UNIQUE-ID") or [None])[0]
            if not tu:
                continue
            genes = [
                c for c in rec.get("COMPONENTS", [])
                if c.startswith("EG") or (c.startswith("G") and "-" in c)
            ]
            self.tu_to_genes[tu] = genes

        for rec in self._iter_dat_records(self.promoters_dat):
            pm = (rec.get("UNIQUE-ID") or [None])[0]
            if not pm:
                continue
            tus = [c for c in rec.get("COMPONENT-OF", []) if c.startswith("TU")]
            self.promoter_to_tus[pm] = tus

        self.gene_incoming: Dict[str, Dict[str, set]] = {}
        for rec in self._iter_dat_records(self.regulation_dat):
            reg = (rec.get("REGULATOR") or [None])[0]
            if not reg:
                continue
            self.regulations_by_regulator.setdefault(reg, []).append(rec)
            entity = (rec.get("REGULATED-ENTITY") or [None])[0]
            if not entity:
                continue
            mode = (rec.get("MODE") or [""])[0]
            effect = (
                "activates"
                if mode == "+"
                else ("inhibits" if mode == "-" else "regulates")
            )
            for g in self._genes_from_regulated_entity(entity):
                bucket = self.gene_incoming.setdefault(
                    g, {"activates": set(), "inhibits": set(), "regulates": set()}
                )
                bucket[effect].add(reg)

        self._incoming_summary_cache = {}

        self.gene_to_uniprot = _parse_gene_links_uniprot(
            os.path.join(self.pgdb_dir, "data", "gene-links.dat")
        )

    def uniprot_accession_for_graph_node(self, node_id: str) -> Optional[str]:
        """Swiss-Prot / UniProt accession for an EcoCyc gene or its polypeptide product, if known."""
        if not node_id:
            return None
        if node_id in self.gene_records:
            return self.gene_to_uniprot.get(node_id)
        gid = self.protein_to_gene.get(node_id)
        if gid:
            return self.gene_to_uniprot.get(gid)
        return None

    def _genes_from_regulated_entity(self, entity: str) -> List[str]:
        """Gene frame IDs whose transcription is touched by this regulated-entity slot."""
        if not entity:
            return []
        if entity.startswith("EG") or (entity.startswith("G") and "-" in entity):
            return [entity]
        if entity.startswith("TU"):
            return list(self.tu_to_genes.get(entity, []))
        if entity.startswith("PM"):
            out = []
            for tu in self.promoter_to_tus.get(entity, []):
                out.extend(self.tu_to_genes.get(tu, []))
            return out
        return []

    def incoming_regulation_summary(
        self, obj_id: str, per_bucket: int = 14
    ) -> Dict[str, List[str]]:
        """EcoCyc: regulators (common name + frame id) grouped by effect on this gene/protein."""
        cache_key = (obj_id, per_bucket)
        if cache_key in self._incoming_summary_cache:
            return self._incoming_summary_cache[cache_key]

        gid = obj_id if obj_id in self.gene_records else self.protein_to_gene.get(obj_id)
        empty = {"activates": [], "inhibits": [], "regulates": []}
        if not gid:
            self._incoming_summary_cache[cache_key] = empty
            return empty
        raw = self.gene_incoming.get(gid)
        if not raw:
            self._incoming_summary_cache[cache_key] = empty
            return empty

        def fmt_regs(regs: set) -> List[str]:
            lines = []
            for r in sorted(regs, key=lambda x: self._record_common_name(x).lower())[
                :per_bucket
            ]:
                cn = self._record_common_name(r)
                lines.append(f"{cn} ({r})")
            return lines

        out = {
            "activates": fmt_regs(raw.get("activates", set())),
            "inhibits": fmt_regs(raw.get("inhibits", set())),
            "regulates": fmt_regs(raw.get("regulates", set())),
        }
        self._incoming_summary_cache[cache_key] = out
        return out

    def _record_common_name(self, obj_id):
        rec = self.gene_records.get(obj_id) or self.protein_records.get(obj_id)
        if not rec:
            return obj_id
        return (rec.get("COMMON-NAME") or [obj_id])[0]

    def _record_product(self, gene_id):
        rec = self.gene_records.get(gene_id)
        if not rec:
            return None
        return (rec.get("PRODUCT") or [None])[0]

    def _candidate_regulator_ids(self, obj_id):
        """Frame IDs to scan as REGULATOR in regulation.dat for this protein/gene product.

        EcoCyc often lists the active form (e.g. CPLX0-226 CRP-cAMP) as regulator, linked from
        the polypeptide via COMPONENT-OF chains (PD00257 → PC00004 → CPLX0-226). A single hop
        misses those edges (e.g. CRP activating malT).
        """
        ids = set()
        frontier = {obj_id}
        while frontier:
            next_frontier = set()
            for oid in frontier:
                if oid in ids:
                    continue
                ids.add(oid)
                rec = self.protein_records.get(oid)
                if not rec:
                    continue
                for parent in rec.get("COMPONENT-OF", []):
                    if parent and parent not in ids:
                        next_frontier.add(parent)
            frontier = next_frontier
        return ids

    def _regulated_genes_for_regulator(self, regulator_id, mode_filter=None):
        out = set()
        for rec in self.regulations_by_regulator.get(regulator_id, []):
            mode = (rec.get("MODE") or [""])[0]
            if mode_filter is not None and mode != mode_filter:
                continue
            entity = (rec.get("REGULATED-ENTITY") or [None])[0]
            if not entity:
                continue

            if entity.startswith("EG") or (entity.startswith("G") and "-" in entity):
                out.add(entity)
                continue
            if entity.startswith("TU"):
                out.update(self.tu_to_genes.get(entity, []))
                continue
            if entity.startswith("PM"):
                for tu in self.promoter_to_tus.get(entity, []):
                    out.update(self.tu_to_genes.get(tu, []))
        return sorted(out)

    @staticmethod
    def _gene_ids_to_xml(gene_ids):
        root = ET.Element("ptools-xml")
        results = ET.SubElement(root, "results")
        for gid in sorted(set(gene_ids)):
            ET.SubElement(results, "Gene", frameid=gid)
        return ET.tostring(root, encoding="unicode")

    def api(self, fn, obj_id):
        try:
            regulators = self._candidate_regulator_ids(obj_id)
            gene_ids = set()
            if fn == "direct-inhibitors":
                for reg in regulators:
                    gene_ids.update(self._regulated_genes_for_regulator(reg, mode_filter="-"))
            elif fn == "direct-activators":
                for reg in regulators:
                    gene_ids.update(self._regulated_genes_for_regulator(reg, mode_filter="+"))
            elif fn == "regulon-of-protein":
                for reg in regulators:
                    gene_ids.update(self._regulated_genes_for_regulator(reg, mode_filter=None))
            else:
                return None
            return self._gene_ids_to_xml(sorted(gene_ids))
        except Exception as e:
            print(f"Local API error ({fn}): {e}")
            return None

    def get_instance(self, obj_id):
        rec = self.gene_records.get(obj_id) or self.protein_records.get(obj_id)
        if not rec:
            return None
        root = ET.Element("ptools-xml")
        entity = ET.SubElement(root, "item", frameid=obj_id)
        ET.SubElement(entity, "common-name").text = (rec.get("COMMON-NAME") or [obj_id])[0]
        product = (rec.get("PRODUCT") or [None])[0]
        if product:
            prod_el = ET.SubElement(entity, "product")
            ET.SubElement(prod_el, "Protein", frameid=product)
        return ET.tostring(root, encoding="unicode")

    def protein_of_gene(self, gene_id):
        return self._record_product(gene_id)

    def name_search_genes(self, name, org_id=ORG_ID):
        """Local name-search: map a gene name/accession to gene frame IDs (JSON-like dict)."""
        _ = org_id
        hits = []
        for gid in sorted(self.name_index.get(name.lower(), set())):
            rec = self.gene_records.get(gid, {})
            hits.append({
                "OBJECT-ID": gid,
                "COMMON-NAME": (rec.get("COMMON-NAME") or [""])[0],
            })
        return {"RESULTS": hits}


def _parse_gene_links_uniprot(path: str) -> Dict[str, str]:
    """Map EcoCyc gene UNIQUE-ID -> UniProt accession from ``gene-links.dat`` (optional)."""
    out: Dict[str, str] = {}
    if not os.path.exists(path):
        return out
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            gid = (parts[0] or "").strip()
            acc = (parts[5] or "").strip()
            if gid and acc:
                out[gid] = acc
    return out


def ecocyc_node_uniprot_entry_url(client: EcoCycClient, node_id: str) -> Optional[str]:
    """Direct UniProt entry URL for EcoCyc gene / polypeptide nodes when ``gene-links.dat`` has an ID."""
    if node_id.startswith(GENE_NODE_PREFIX) or node_id.startswith(UNIPROT_NODE_PREFIX):
        return None
    acc = client.uniprot_accession_for_graph_node(node_id)
    if acc:
        return uniprot_entry_url(acc)
    return None


# =========================
# XML Parsing
# =========================
def parse_gene_ids(xml):
    """Extract gene frame IDs (EGxxxx)"""
    ids = set()
    if not xml:
        return []

    root = ET.fromstring(xml)

    for el in root.findall(".//*[@frameid]"):
        fid = el.get("frameid")
        if fid and fid.startswith("EG"):
            ids.add(fid)

    return list(ids)


def get_common_name(client, obj_id):
    xml = client.get_instance(obj_id)
    if not xml:
        return obj_id
    try:
        root = ET.fromstring(xml)
        name = root.findtext(".//common-name")
        return name if name else obj_id
    except:
        return obj_id


def get_node_label(client, obj_id, uniprot_labels=None):
    if uniprot_labels and obj_id in uniprot_labels:
        return uniprot_labels[obj_id]
    return get_common_name(client, obj_id)


def _looks_like_biocyc_gene_frame_id(token):
    t = token.strip()
    if re.match(r"^EG\d+$", t, re.I):
        return True
    if re.match(r"^G\d+-\d+$", t, re.I):
        return True
    return False


def _pick_gene_search_result(results, query):
    """Prefer case-insensitive exact match on common name or frame id."""
    q = query.strip().lower()
    for r in results:
        oid = (r.get("OBJECT-ID") or "").strip()
        cn = (r.get("COMMON-NAME") or "").strip()
        if oid.lower() == q:
            return oid
        if cn.lower() == q:
            return oid
    for r in results:
        cn = (r.get("COMMON-NAME") or "").strip()
        if cn and q in cn.lower():
            return r.get("OBJECT-ID")
    return (results[0].get("OBJECT-ID") or "").strip()


def resolve_gene_input(client, token, org_id=ORG_ID) -> ResolvedGene:
    """
    Resolve to EcoCyc gene id, or to UniProt (gene symbol + encoded protein + PPI; no EcoCyc mapping).
    Accepts frame ids (e.g. EG10525, G0-10441) or names (e.g. LacI, trpA).
    """
    t = token.strip()
    if not t:
        raise ValueError("empty gene token")
    if _looks_like_biocyc_gene_frame_id(t):
        gid = t.upper() if re.match(r"^EG\d+$", t, re.I) else t
        return ResolvedGene(source="ecocyc", ecocyc_id=gid)

    data = client.name_search_genes(t, org_id)
    if not data:
        raise RuntimeError(f"name search failed for {t!r}")
    results = data.get("RESULTS") or []
    if not results:
        print(f"No local EcoCyc hit for {t!r}; resolving gene via UniProt, then encoded protein + interactions...")
        hit = uniprot_search_best_hit(t)
        if not hit:
            raise ValueError(f"No EcoCyc gene and no UniProt hit for {t!r}.")
        acc, gene_sym, prot_label = hit
        print(
            f"Resolved {t!r} -> gene {gene_sym} encodes {acc} ({prot_label})"
        )
        return ResolvedGene(
            source="uniprot",
            uniprot_accession=acc,
            uniprot_gene_symbol=gene_sym,
            protein_label=prot_label,
        )

    if len(results) > 1:
        print(f"\nMultiple matches for {t!r} ({len(results)} results):")
        for r in results:
            oid = r.get("OBJECT-ID", "?")
            cn = r.get("COMMON-NAME", "")
            print(f"  {oid}  ({cn})")

    chosen = _pick_gene_search_result(results, t)
    if not chosen:
        raise ValueError(f"Could not pick a gene id for {t!r}")
    if len(results) > 1:
        print(f"Using gene id: {chosen} (disambiguate by passing EG… id if wrong)")
    return ResolvedGene(source="ecocyc", ecocyc_id=chosen)

