import io
import os
import re
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout, as_completed
from typing import Any, Dict, List, Optional, Tuple

import urllib.error
import urllib.parse
import urllib.request

try:
    from Bio.Blast import NCBIXML
except ImportError:
    NCBIXML = None  # type: ignore
try:
    from Bio.Blast import NCBIWWW
except ImportError:
    NCBIWWW = None  # type: ignore

SEQ_MIN_ORF_AA = 50
SEQ_MAX_BLAST_WORKERS = 3
SEQ_REMOTE_BLAST_DELAY_SEC = float(os.environ.get("SEQ_REMOTE_BLAST_DELAY_SEC", "4"))
SEQ_REMOTE_ORF_CAP = max(1, int(os.environ.get("SEQ_REMOTE_ORF_CAP", "5")))
SEQ_REMOTE_BLAST_TIMEOUT_SEC = max(30, int(os.environ.get("SEQ_REMOTE_BLAST_TIMEOUT_SEC", "180")))
SEQ_BLAST_MIN_IDENTITY_DEFAULT = float(os.environ.get("SEQ_BLAST_MIN_IDENTITY_PCT", "75"))

_DNA_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

_GENETIC_CODE_11: Dict[str, str] = {}
_gc11_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_rna_order = [
    "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
    "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
    "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
    "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
    "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
    "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
    "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
    "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG",
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
    client: Optional["EcoCycClient"] = None,
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


# Late import to avoid circular dependency at import time.
def resolve_blast_tokens_to_graph_seeds(client, tokens):
    from graph_builder import resolve_blast_tokens_to_graph_seeds as _impl
    return _impl(client, tokens)
