import argparse
import io
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout, as_completed
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import networkx as nx

try:
    from Bio.Blast import NCBIXML
except ImportError:
    NCBIXML = None  # type: ignore
try:
    from Bio.Blast import NCBIWWW
except ImportError:
    NCBIWWW = None  # type: ignore

ORG_ID = "ECOLI"
DEFAULT_PGDB_DIR = "29.6"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb"
UNIPROT_ORGANISM_ID = "83333"  # E. coli K-12
UNIPROT_TAXON_ECOLI = "562"
UNIPROT_NODE_PREFIX = "UniProt:"
GENE_NODE_PREFIX = "Gene:"
STRINGDB_API_BASE = os.environ.get("STRINGDB_API_BASE", "https://version-12-0.string-db.org/api/json")
STRINGDB_DEFAULT_SPECIES = int(os.environ.get("STRINGDB_SPECIES", "511145"))
STRINGDB_DEFAULT_MIN_SCORE = int(os.environ.get("STRINGDB_MIN_SCORE", "700"))
STRINGDB_MAX_PARTNERS = max(1, int(os.environ.get("STRINGDB_MAX_PARTNERS", "40")))

# ── DNA / ORF / local BLAST (in-repo; replaces former Protein.py) ─────────
SEQ_MIN_ORF_AA = 50
SEQ_MAX_BLAST_WORKERS = 3
# NCBI web BLAST (when local blastp is missing): delay between jobs, max ORFs per request
SEQ_REMOTE_BLAST_DELAY_SEC = float(os.environ.get("SEQ_REMOTE_BLAST_DELAY_SEC", "4"))
SEQ_REMOTE_ORF_CAP = max(1, int(os.environ.get("SEQ_REMOTE_ORF_CAP", "5")))
SEQ_REMOTE_BLAST_TIMEOUT_SEC = max(30, int(os.environ.get("SEQ_REMOTE_BLAST_TIMEOUT_SEC", "180")))
# Minimum aligned identity (%) on the best BLAST HSP to promote a hit to suggested gene names
# (same idea as ``Protein.py`` high-confidence filter; default 75).
SEQ_BLAST_MIN_IDENTITY_DEFAULT = float(os.environ.get("SEQ_BLAST_MIN_IDENTITY_PCT", "75"))

_DNA_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

# NCBI genetic code table 11 (bacterial), codons in NCBI nucleotide order
_GENETIC_CODE_11: Dict[str, str] = {}
_gc11_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_rna_order = [
    "TTT",
    "TTC",
    "TTA",
    "TTG",
    "TCT",
    "TCC",
    "TCA",
    "TCG",
    "TAT",
    "TAC",
    "TAA",
    "TAG",
    "TGT",
    "TGC",
    "TGA",
    "TGG",
    "CTT",
    "CTC",
    "CTA",
    "CTG",
    "CCT",
    "CCC",
    "CCA",
    "CCG",
    "CAT",
    "CAC",
    "CAA",
    "CAG",
    "CGT",
    "CGC",
    "CGA",
    "CGG",
    "ATT",
    "ATC",
    "ATA",
    "ATG",
    "ACT",
    "ACC",
    "ACA",
    "ACG",
    "AAT",
    "AAC",
    "AAA",
    "AAG",
    "AGT",
    "AGC",
    "AGA",
    "AGG",
    "GTT",
    "GTC",
    "GTA",
    "GTG",
    "GCT",
    "GCC",
    "GCA",
    "GCG",
    "GAT",
    "GAC",
    "GAA",
    "GAG",
    "GGT",
    "GGC",
    "GGA",
    "GGG",
]
for _i, _codon in enumerate(_rna_order):
    _GENETIC_CODE_11[_codon] = _gc11_aa[_i]


def clean_dna_sequence(raw: str) -> str:
    return re.sub(r"[^ACGTacgt]", "", raw or "").upper()


def normalize_dna_clipboard_input(text: str) -> str:
    """Strip optional Python ``RAW_DNA = \"\"\" ... \"\"\"`` wrapper from pasted snippets."""
    t = (text or "").strip()
    m = re.search(r'RAW_DNA\s*=\s*"""(.*?)"""', t, re.DOTALL | re.IGNORECASE)
    if m:
        return m.group(1).strip()
    return t


def get_reverse_complement(dna: str) -> str:
    return clean_dna_sequence(dna).translate(_DNA_COMPLEMENT)[::-1]


def parse_fasta_text(text: str) -> List[Tuple[str, str]]:
    """Return list of (record_id, ACGT DNA) from FASTA or raw sequence."""
    text = normalize_dna_clipboard_input(text)
    if not text.strip():
        return []
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    records: List[Tuple[str, str]] = []
    cur_id = "sequence"
    buf: List[str] = []
    for ln in lines:
        if ln.startswith(">"):
            if buf:
                seq = clean_dna_sequence("".join(buf))
                if seq:
                    records.append((cur_id, seq))
            cur_id = ln[1:].strip().split()[0] or "sequence"
            buf = []
        else:
            buf.append(ln)
    if buf:
        seq = clean_dna_sequence("".join(buf))
        if seq:
            records.append((cur_id, seq))
    if not records and not any(ln.startswith(">") for ln in lines):
        seq = clean_dna_sequence("".join(lines))
        if seq:
            records.append(("sequence", seq))
    return records


def _pbr322_plasmid_gene_hints(record_ids: List[str], longest_seq_len: int) -> Tuple[List[str], List[str]]:
    """
    Curated pBR322 **protein** symbols only: ``bla`` (AmpR / beta-lactamase) and ``tetA`` (TcR).

    ORF-level BLAST often returns spurious cross-database gene names on plasmid ORFs; when this
    heuristic fires, those hits are **not** promoted into ``suggested_genes`` (see per-ORF
    ``blast_gene`` fields if you need them). **ColE1 ori**, **bom**, and the **AmpR promoter** are
    non-coding; add them manually if you want them beside the protein list.
    """
    hints: List[str] = []
    notes: List[str] = []
    joined = " ".join(record_ids).lower()
    name_hit = bool(re.search(r"pbr322|pbr-322|pbr_322", joined))
    len_hit = bool(longest_seq_len and 4310 <= longest_seq_len <= 4435)
    if not (name_hit or len_hit):
        return hints, notes
    hints = ["bla", "tetA"]
    notes.append(
        "pBR322-like sequence (FASTA id or length ~4.3 kb): suggested genes are **bla** (AmpR) and "
        "**tetA** (TcR) only. ColE1 **ori** (replication origin) is not a resolvable protein symbol "
        "here - add the label \"ori\" or \"ColE1 ori\" yourself if you track it next to the list. "
        "Per-ORF BLAST names remain under each ORF below for debugging; they are not auto-imported."
    )
    return hints, notes


def _translate_orf_dna(dna_nt: str) -> str:
    """Translate DNA from first base through in-frame codons; stops before stop codon (excluded)."""
    dna_nt = clean_dna_sequence(dna_nt)
    aa: List[str] = []
    for i in range(0, len(dna_nt) - 2, 3):
        codon = dna_nt[i : i + 3]
        if codon in _STOP_CODONS:
            break
        aa.append(_GENETIC_CODE_11.get(codon, "X"))
    return "".join(aa)


def find_orfs(dna: str, strand: str, seq_len: int) -> List[Dict[str, Any]]:
    """
    ORFs with in-frame ATG and downstream in-frame stop; one longest ORF per stop (same strand).
    ``strand`` is '+' (dna is forward) or '-' (dna is reverse complement of record).
    """
    dna = clean_dna_sequence(dna)
    n = len(dna)
    best_per_stop: Dict[int, Dict[str, Any]] = {}
    for frame_idx in range(3):
        frame_label = frame_idx + 1 if strand == "+" else -(frame_idx + 1)
        i = frame_idx
        while i + 3 <= n:
            if dna[i : i + 3] != "ATG":
                i += 3
                continue
            start = i
            j = i + 3
            while j + 3 <= n:
                codon = dna[j : j + 3]
                if codon in _STOP_CODONS:
                    stop_pos = j
                    cds = dna[i:j]
                    protein = _translate_orf_dna(cds)
                    if len(protein) < SEQ_MIN_ORF_AA:
                        break
                    cand = {
                        "frame": frame_label,
                        "start": start + 1,
                        "length": len(protein),
                        "one_letter": protein,
                    }
                    prev = best_per_stop.get(stop_pos)
                    if prev is None or len(protein) > prev["length"]:
                        best_per_stop[stop_pos] = cand
                    break
                j += 3
            i += 3
    return list(best_per_stop.values())


def parse_swissprot_gene_from_defline(defline: str) -> str:
    """
    Prefer GN= from UniProt-style text; else compact ``GENE_ORG`` after second ``|``;
    else primary accession (e.g. ``P42212``) so ``resolve_gene_input`` can hit UniProt.
    """
    if not defline:
        return ""
    m = re.search(r"\bGN=([A-Za-z0-9_.-]+)\b", defline)
    if m:
        g = m.group(1).strip()
        if g and g.upper() != "N/A":
            return g
    parts = [p.strip() for p in defline.split("|")]
    if len(parts) >= 3 and parts[2]:
        third = parts[2]
        if "_" in third and len(third) < 48 and not third.lower().startswith("recname"):
            return third.split("_")[0].strip()
    # sp|ACCESSION.version|... or tr|...
    m2 = re.match(r"^[a-z]{2}\|([^|]+)\|", defline.strip(), re.I)
    if m2:
        acc = m2.group(1).split(".")[0].strip()
        if re.match(r"^[A-Z0-9]{4,12}$", acc, re.I):
            return acc
    return ""


def resolve_blastp_executable() -> Optional[str]:
    override = os.environ.get("BLASTP_EXE", "").strip()
    if override and os.path.isfile(override):
        return override
    return shutil.which("blastp") or shutil.which("blastp.exe")


def resolve_swissprot_db_prefix() -> str:
    for key in ("SWISSPROT_DB", "BLAST_DB"):
        v = (os.environ.get(key) or "").strip().strip('"')
        if v:
            return v
    root = os.path.dirname(os.path.abspath(__file__))
    local = os.path.join(root, "swissprot", "swissprot")
    if os.path.isfile(local + ".pin") or os.path.isfile(local + ".pal"):
        return local
    return "swissprot"


def _blast_db_argv_and_cwd(db_path: str) -> Tuple[List[str], Optional[str]]:
    db_path = os.path.normpath(db_path)
    if not db_path or db_path == "swissprot":
        return ["-db", "swissprot"], None
    parent = os.path.dirname(db_path)
    base = os.path.basename(db_path)
    if parent and base and (
        " " in db_path
        and (
            os.path.isfile(os.path.join(parent, base + ".pin"))
            or os.path.isfile(os.path.join(parent, base + ".pal"))
        )
    ):
        return ["-db", base], parent
    return ["-db", db_path], None


_blastp_missing_logged = False


def protein_blast_env_notes() -> List[str]:
    notes: List[str] = []
    root = os.path.dirname(os.path.abspath(__file__))
    sp_dir = os.path.join(root, "swissprot")
    fasta = os.path.join(sp_dir, "uniprot_sprot.fasta")
    prefix_default = os.path.join(sp_dir, "swissprot")
    if os.path.isfile(fasta) and not (
        os.path.isfile(prefix_default + ".pin") or os.path.isfile(prefix_default + ".pal")
    ):
        notes.append(
            "Found swissprot/uniprot_sprot.fasta but no BLAST index (*.pin). "
            "cd into swissprot, then: makeblastdb -in uniprot_sprot.fasta -dbtype prot -out swissprot"
        )
    if not resolve_blastp_executable():
        if NCBIWWW is not None and os.environ.get("SEQ_REMOTE_BLAST", "1").strip().lower() not in (
            "0",
            "false",
            "no",
        ):
            notes.append(
                "Local blastp.exe not found (optional for speed/offline). "
                "The app will try NCBI web BLASTp when you import a sequence (needs internet). "
                "To use local BLAST instead: install NCBI BLAST+, then set BLASTP_EXE to the full path, "
                "e.g. C:\\Program Files\\NCBI\\blast-2.16.0+\\bin\\blastp.exe"
            )
        else:
            notes.append(
                "blastp not found: install NCBI BLAST+ and add bin to PATH, "
                "or set BLASTP_EXE to the full path of blastp.exe. "
                "Or enable web BLAST (default on): ensure biopython includes Bio.Blast.NCBIWWW "
                "and do not set SEQ_REMOTE_BLAST=0."
            )
    db_path = resolve_swissprot_db_prefix()
    if db_path and db_path != "swissprot":
        db_norm = os.path.normpath(db_path)
        parent = os.path.dirname(db_norm)
        base = os.path.basename(db_norm)
        candidates = [db_norm + ".pin", db_norm + ".pal"]
        if parent and base:
            candidates.extend(
                [os.path.join(parent, base + ".pin"), os.path.join(parent, base + ".pal")]
            )
        if not any(os.path.isfile(p) for p in candidates):
            notes.append(
                f"SWISSPROT_DB/BLAST_DB is {db_path!r} but no .pin/.pal at that prefix."
            )
    if NCBIXML is None:
        notes.append("Install biopython for BLAST XML parsing: pip install biopython")
    return notes


def blast_protein_local(one_letter_seq: str, orf_id: int) -> Optional[Dict[str, Any]]:
    """Run local blastp vs SwissProt-formatted DB; return best hit metadata or None."""
    global _blastp_missing_logged
    if NCBIXML is None:
        return None
    exe = resolve_blastp_executable()
    if not exe:
        if not _blastp_missing_logged:
            print(
                "Local BLAST: blastp not found. Install NCBI BLAST+ or set BLASTP_EXE.",
                flush=True,
            )
            _blastp_missing_logged = True
        return None

    db_path = resolve_swissprot_db_prefix()
    db_args, blast_cwd = _blast_db_argv_and_cwd(db_path)
    fd, temp_query = tempfile.mkstemp(prefix="orf_", suffix=".fasta", text=True)
    try:
        with os.fdopen(fd, "w") as f:
            f.write(f">ORF_{orf_id}\n{one_letter_seq}\n")
        cmd = [exe, "-query", temp_query] + db_args + ["-outfmt", "5", "-max_target_seqs", "1"]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            cwd=blast_cwd,
        )
        if result.returncode != 0:
            err = (result.stderr or result.stdout or "").strip()[:400]
            print(f"Local BLAST: blastp exit {result.returncode} (db={db_path!r}). {err}", flush=True)
            return None
        for record in NCBIXML.parse(io.StringIO(result.stdout)):
            if record.alignments:
                best_algn = record.alignments[0]
                best_hsp = best_algn.hsps[0]
                gn = parse_swissprot_gene_from_defline(best_algn.title)
                return {
                    "orf_id": orf_id,
                    "gene_name": gn,
                    "name": best_algn.title,
                    "identity": round(
                        (best_hsp.identities / best_hsp.align_length) * 100, 1
                    ),
                    "e_value": f"{best_hsp.expect:.2e}",
                    "score": best_hsp.score,
                }
    except OSError as e:
        if getattr(e, "winerror", None) == 2 or e.errno == 2:
            if not _blastp_missing_logged:
                print(
                    "Local BLAST: could not run blastp. Set BLASTP_EXE to blastp.exe if needed.",
                    flush=True,
                )
                _blastp_missing_logged = True
        else:
            print(f"Local BLAST Error: {e}", flush=True)
    except Exception as e:
        print(f"Local BLAST Error: {e}", flush=True)
    finally:
        try:
            if os.path.exists(temp_query):
                os.remove(temp_query)
        except OSError:
            pass
    return None


def _remote_ncbi_blast_enabled() -> bool:
    if NCBIXML is None or NCBIWWW is None:
        return False
    v = (os.environ.get("SEQ_REMOTE_BLAST") or "1").strip().lower()
    return v not in ("0", "false", "no")


def _ncbi_entrez_email() -> str:
    return (
        (os.environ.get("NCBI_EMAIL") or os.environ.get("ENTREZ_EMAIL") or "").strip()
        or "autogene-sequence@users.noreply.github.com"
    )


def blast_protein_remote_ncbi(one_letter_seq: str, orf_id: int) -> Optional[Dict[str, Any]]:
    """
    Run blastp vs SwissProt on NCBI's servers (no local blastp). Slower; needs network.
    Set NCBI_EMAIL or ENTREZ_EMAIL to your address for NCBI's usage guidelines.
    """
    if not _remote_ncbi_blast_enabled():
        return None
    seq = (one_letter_seq or "").strip()
    if len(seq) < 12:
        return None

    def _qblast():
        return NCBIWWW.qblast(
            "blastp",
            "swissprot",
            seq,
            hitlist_size=1,
            expect=10.0,
        )

    try:
        from Bio import Entrez

        Entrez.email = _ncbi_entrez_email()
        with ThreadPoolExecutor(max_workers=1) as pool:
            fut = pool.submit(_qblast)
            handle = fut.result(timeout=SEQ_REMOTE_BLAST_TIMEOUT_SEC)
    except FuturesTimeout:
        print(f"Remote BLAST: timeout after {SEQ_REMOTE_BLAST_TIMEOUT_SEC}s (ORF {orf_id})", flush=True)
        return None
    except Exception as e:
        print(f"Remote BLAST error (ORF {orf_id}): {e}", flush=True)
        return None

    try:
        raw = handle.read()
        if hasattr(handle, "close"):
            handle.close()
        if isinstance(raw, bytes):
            xml = raw.decode("utf-8", errors="replace")
        else:
            xml = str(raw)
        record = NCBIXML.read(io.StringIO(xml))
    except Exception as e:
        print(f"Remote BLAST XML parse error (ORF {orf_id}): {e}", flush=True)
        return None

    if not record.alignments:
        return None
    best_algn = record.alignments[0]
    if not best_algn.hsps:
        return None
    best_hsp = best_algn.hsps[0]
    gn = parse_swissprot_gene_from_defline(best_algn.title)
    return {
        "orf_id": orf_id,
        "gene_name": gn,
        "name": best_algn.title,
        "identity": round((best_hsp.identities / best_hsp.align_length) * 100, 1),
        "e_value": f"{best_hsp.expect:.2e}",
        "score": best_hsp.score,
        "blast_source": "ncbi_web",
    }


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


@dataclass
class ResolvedGene:
    """EcoCyc gene id, or UniProt gene + encoded protein (PPI path; no EcoCyc mapping)."""

    source: str  # "ecocyc" | "uniprot"
    ecocyc_id: Optional[str] = None
    uniprot_accession: Optional[str] = None
    uniprot_gene_symbol: Optional[str] = None
    protein_label: Optional[str] = None


# =========================
# EcoCyc Client
# =========================
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


def seed_graph_node_id(r: ResolvedGene) -> str:
    """Primary graph node id for a resolved input (matches build_graph anchors)."""
    if r.source == "ecocyc" and r.ecocyc_id:
        return r.ecocyc_id
    if r.source == "uniprot" and r.uniprot_accession:
        gene_sym = (r.uniprot_gene_symbol or r.uniprot_accession).strip()
        return uniprot_gene_node_id(gene_sym)
    return ""


def _swissprot_defline_gene_symbols(title: str) -> List[str]:
    """Gene symbols from Swiss-Prot / UniProt style deflines (e.g. GN=lacZ) when BLAST omits gene_name."""
    if not title:
        return []
    seen: set = set()
    out: List[str] = []
    for m in re.finditer(r"\bGN=([A-Za-z0-9_.-]+)\b", title):
        g = m.group(1).strip()
        if not g or g.upper() == "N/A":
            continue
        key = g.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(g)
    return out


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
    client: EcoCycClient, tokens: List[str]
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


def _sequence_hints_empty(*note_lines: str) -> Dict[str, Any]:
    """Uniform JSON shape for sequence analysis API (including errors / empty input)."""
    return {
        "suggested_genes": [],
        "orfs": [],
        "notes": [n for n in note_lines if n],
        "resolved_seeds": [],
        "resolution_errors": [],
        "genes_ready_for_graph": [],
    }


def analyze_sequences_for_graph_hints(
    dna_fasta_text: str,
    max_orfs_to_blast: int = 10,
    client: Optional[EcoCycClient] = None,
    resolve_blast_genes: bool = True,
    min_blast_identity_pct: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Parse FASTA or raw DNA, predict ORFs, run local BLASTp (when ``blastp`` and a
    SwissProt-style DB are available) to propose gene symbols, then optionally
    resolve each symbol through EcoCyc (local PGDB) or UniProt for
    ``build_graph`` / ``resolve_gene_input``.

    Requires optional ``blastp`` + ``makeblastdb`` SwissProt indices and biopython
    (``Bio.Blast.NCBIXML``) for BLAST parsing. ORFs are always returned when DNA
    parses successfully.

    When ``client`` is set and ``resolve_blast_genes`` is True, the response includes
    ``resolved_seeds``, ``resolution_errors``, and ``genes_ready_for_graph``.

    Best BLAST alignments must meet ``min_blast_identity_pct`` (default from
    ``SEQ_BLAST_MIN_IDENTITY_PCT`` / 75, matching ``Protein.py``-style high-confidence
    filtering) before their gene symbols are added to ``suggested_genes``.

    pBR322-like inserts (FASTA id or length ~4.3 kb) restrict ``suggested_genes`` to ``bla`` and
    ``tetA`` (AmpR and TcR); ORF BLAST hits stay on each ORF row only. **ori** / **bom** / promoters
    are non-coding and must be added manually if needed.
    """
    text = (dna_fasta_text or "").strip()
    if not text:
        return _sequence_hints_empty("Empty sequence.")

    records = parse_fasta_text(text)
    if not records:
        return _sequence_hints_empty("No DNA records parsed.")

    eff_min_identity = (
        float(min_blast_identity_pct)
        if min_blast_identity_pct is not None
        else float(SEQ_BLAST_MIN_IDENTITY_DEFAULT)
    )
    eff_min_identity = max(0.0, min(100.0, eff_min_identity))

    notes: List[str] = []
    notes.append(
        f"BLAST: best-hit aligned identity must be >= {eff_min_identity:g}% to add a gene symbol "
        "to the suggested list (see ``Protein.py`` high-confidence pattern)."
    )
    pooled: List[Dict[str, Any]] = []
    blast_counter = 0

    for rec_id, dna in records:
        if len(dna) < SEQ_MIN_ORF_AA * 3:
            notes.append(f"Skipped {rec_id!r}: sequence shorter than minimum ORF length.")
            continue
        seq_len = len(dna)
        orfs: List[Dict] = []
        orfs.extend(find_orfs(dna, "+", seq_len))
        rev = get_reverse_complement(dna)
        orfs.extend(find_orfs(rev, "-", seq_len))
        orfs.sort(key=lambda x: x["length"], reverse=True)
        for o in orfs:
            blast_counter += 1
            pooled.append({**o, "record_id": rec_id, "blast_idx": blast_counter})

    pooled.sort(key=lambda x: x["length"], reverse=True)
    if not pooled:
        return _sequence_hints_empty(
            *notes,
            "No ORFs found (need in-frame ATG through stop, min length per "
            f"{SEQ_MIN_ORF_AA} aa).",
        )

    longest_seq = max((len(dna) for _, dna in records), default=0)
    plasmid_hints, plasmid_notes = _pbr322_plasmid_gene_hints(
        [rid for rid, _ in records], longest_seq
    )
    notes.extend(plasmid_notes)

    exe_ok = bool(resolve_blastp_executable())
    remote_ok = (not exe_ok) and _remote_ncbi_blast_enabled()
    eff_max = max(1, max_orfs_to_blast)
    # Short replicons (e.g. pBR322 ~4.3 kb) carry several CDS; include more ORFs so *rop* is not dropped.
    if longest_seq and longest_seq <= 8000:
        eff_max = min(max(eff_max, 18), len(pooled), 25)
    if remote_ok:
        remote_cap = (
            max(SEQ_REMOTE_ORF_CAP, 15)
            if (longest_seq and longest_seq <= 8000)
            else SEQ_REMOTE_ORF_CAP
        )
        eff_max = min(eff_max, remote_cap)
        if eff_max < max_orfs_to_blast:
            notes.append(
                f"Web BLAST: only the top {eff_max} longest ORFs are queried (cap={remote_cap}; "
                "install local blastp to BLAST more ORFs per request)."
            )
    to_blast = pooled[:eff_max]
    blast_results: Dict[int, Tuple[Optional[Dict[str, Any]], Dict[str, Any]]] = {}
    workers = max(1, min(SEQ_MAX_BLAST_WORKERS, len(to_blast)))

    if exe_ok:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            future_to_idx = {}
            for i, o in enumerate(to_blast):
                fut = pool.submit(
                    blast_protein_local, o["one_letter"], int(o["blast_idx"])
                )
                future_to_idx[fut] = (i, o)
            for fut in as_completed(future_to_idx):
                i, o = future_to_idx[fut]
                try:
                    hit = fut.result()
                except Exception:
                    hit = None
                blast_results[i] = (hit, o)
    elif remote_ok:
        notes.append(
            "Using NCBI web BLASTp (no local blastp). Requires internet; can take a minute. "
            "Optional: set NCBI_EMAIL to your email for NCBI policy. "
            "For faster/offline use, install BLAST+ and set BLASTP_EXE."
        )
        for i, o in enumerate(to_blast):
            if i > 0:
                time.sleep(SEQ_REMOTE_BLAST_DELAY_SEC)
            try:
                hit = blast_protein_remote_ncbi(o["one_letter"], int(o["blast_idx"]))
            except Exception:
                hit = None
            blast_results[i] = (hit, o)
    else:
        for i, o in enumerate(to_blast):
            blast_results[i] = (None, o)

    orf_out: List[Dict[str, Any]] = []
    genes_order: List[str] = []
    seen: set = set()

    for i in range(len(to_blast)):
        hit, o = blast_results[i]
        entry: Dict[str, Any] = {
            "record_id": o["record_id"],
            "frame": o["frame"],
            "start": o["start"],
            "length_aa": o["length"],
            "protein_prefix": (o.get("one_letter") or "")[:48],
        }
        if hit:
            title = (hit.get("name") or "").strip()
            gn = (hit.get("gene_name") or "").strip()
            if not gn or gn.upper() == "N/A":
                for alt in _swissprot_defline_gene_symbols(title):
                    gn = alt
                    break
            id_pct = hit.get("identity")
            try:
                id_f = float(id_pct) if id_pct is not None else None
            except (TypeError, ValueError):
                id_f = None
            pass_identity = id_f is None or id_f >= eff_min_identity - 1e-9
            entry["blast_hit_title"] = title[:160]
            entry["identity_pct"] = id_pct
            entry["e_value"] = hit.get("e_value")
            entry["blast_gene"] = (
                gn if gn and gn.upper() != "N/A" and pass_identity else None
            )
            if gn and gn.upper() != "N/A" and gn not in seen and pass_identity:
                seen.add(gn)
                genes_order.append(gn)
        else:
            entry["blast_gene"] = None
        orf_out.append(entry)

    if plasmid_hints:
        # Do not merge weak / spurious BLAST gene names (e.g. RTEL1, galK) into the import list.
        genes_order = list(plasmid_hints)

    resolved_seeds: List[Dict[str, Any]] = []
    resolution_errors: List[Dict[str, str]] = []
    genes_ready_for_graph: List[str] = []

    if client is not None and resolve_blast_genes and genes_order:
        resolved_seeds, resolution_errors = resolve_blast_tokens_to_graph_seeds(
            client, genes_order
        )
        genes_ready_for_graph = [row["token"] for row in resolved_seeds]
        src = "pBR322 backbone" if plasmid_hints else "BLAST-derived"
        notes.append(
            f"EcoCyc/UniProt: resolved {len(resolved_seeds)}/{len(genes_order)} "
            f"{src} gene name(s) to graph seeds."
        )

    if not genes_order:
        env_notes = protein_blast_env_notes()
        if env_notes:
            notes.extend(env_notes)
        else:
            notes.append(
                "No BLAST gene names for these ORFs: BLAST may have run but returned no "
                "alignments, or the database is not SwissProt protein. For a typical setup, "
                "install NCBI BLAST+, set BLASTP_EXE on Windows if needed, run makeblastdb on "
                "uniprot_sprot.fasta, and set SWISSPROT_DB to the DB prefix (path to "
                "swissprot/swissprot with no extension)."
            )

    out: Dict[str, Any] = {
        "suggested_genes": genes_order,
        "orfs": orf_out,
        "notes": notes,
        "resolved_seeds": resolved_seeds,
        "resolution_errors": resolution_errors,
        "genes_ready_for_graph": genes_ready_for_graph,
        "blast_min_identity_pct": eff_min_identity,
    }
    if plasmid_hints:
        out["plasmid_hint_genes"] = plasmid_hints
        out["plasmid_backbone_labels"] = [
            "AmpR (gene bla)",
            "TcR (gene tetA)",
            "ColE1 ori (non-coding; add manually)",
        ]
    return out


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


# =========================
# Graph Builder
# =========================
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
        seed_symbol_to_id: Dict[str, str] = {}
        seed_to_acc: Dict[str, str] = {}

        for r in resolved:
            sid = seed_graph_node_id(r)
            if not sid:
                continue
            if r.source == "uniprot":
                acc = (r.uniprot_accession or "").strip()
                sym = (r.uniprot_gene_symbol or "").strip().lower()
                if sym:
                    seed_symbol_to_id.setdefault(sym, sid)
            else:
                gid = (r.ecocyc_id or "").strip()
                acc = (client.uniprot_accession_for_graph_node(gid) or "").strip()
                rec = client.gene_records.get(gid) or {}
                for k in ("COMMON-NAME", "SYNONYMS"):
                    for s in rec.get(k, []) or []:
                        sym = str(s).strip().lower()
                        if sym:
                            seed_symbol_to_id.setdefault(sym, sid)
            if acc:
                seed_to_acc[sid] = acc

        for sid, acc in seed_to_acc.items():
            for partner_gene, conf in stringdb_interaction_partners(
                acc,
                species=stringdb_species,
                min_score=stringdb_min_score,
            ):
                tid = seed_symbol_to_id.get((partner_gene or "").strip().lower())
                if not tid or tid == sid:
                    continue
                prev = G.get_edge_data(sid, tid)
                if prev and prev.get("interaction") != "stringdb":
                    continue
                new_conf = round(float(conf), 3)
                if prev and prev.get("interaction") == "stringdb":
                    old_conf = float(prev.get("confidence", 0.0) or 0.0)
                    if old_conf >= new_conf:
                        continue
                G.add_edge(
                    sid,
                    tid,
                    interaction="stringdb",
                    source="stringdb",
                    confidence=new_conf,
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


# =========================
# MAIN
# =========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build EcoCyc regulation graph and export JSON for the circuit viewer.")
    parser.add_argument(
        "--pgdb-dir",
        default=DEFAULT_PGDB_DIR,
        help="Path to local EcoCyc PGDB directory (default: 29.6).",
    )
    parser.add_argument(
        "-o",
        "--output-json",
        default="circuit_data.json",
        help="Path to write graph JSON consumed by circuit_viewer/index.html (default: circuit_data.json).",
    )
    parser.add_argument(
        "genes",
        nargs="*",
        metavar="GENE",
        help='Gene names (e.g. LacI lacZ) or EcoCyc gene ids (e.g. EG10525). '
        "If not in EcoCyc, resolves the gene via UniProt, adds Gene->protein (creates), "
        "then protein-protein edges from UniProt INTERACTION data. "
        "Default: lacI lacZ.",
    )
    args = parser.parse_args()

    client = EcoCycClient(args.pgdb_dir)

    tokens = list(args.genes) if args.genes else ["LacI", "lacZ"]
    resolved: List[ResolvedGene] = []
    for tok in tokens:
        r = resolve_gene_input(client, tok)
        if r.source == "ecocyc":
            print(f"Resolved {tok!r} -> {r.ecocyc_id}")
        resolved.append(r)

    G, uniprot_labels = build_graph(client, resolved)

    print("\n=== Graph Summary ===")
    print("Nodes:", G.number_of_nodes())
    print("Edges:", G.number_of_edges())

    print("\n=== Edges (Readable) ===")
    name_cache = {}

    def cached_name(nid):
        if nid not in name_cache:
            name_cache[nid] = get_node_label(client, nid, uniprot_labels)
        return name_cache[nid]

    for u, v, d in G.edges(data=True):
        print(f"{cached_name(u)} --{d['interaction']}--> {cached_name(v)}")

    export_graph_json(
        client,
        G,
        args.output_json,
        name_cache=name_cache,
        uniprot_labels=uniprot_labels,
    )

