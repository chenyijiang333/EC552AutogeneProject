#!/usr/bin/env python3
"""
DNA to Amino Acid Translator
Reads a DNA sequence, finds all ORFs (Open Reading Frames),
and lists the protein chains encoded in the sequence.
Automatically collects the top BLAST hit for each ORF into an array.

Multithreaded BLAST edition — concurrent requests with NCBI rate-limit compliance.
"""
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
from typing import Optional, List, Dict, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import time
import pprint

# ─────────────────────────────────────────────
# NCBI CONFIGURATION ← required for remote access
# ─────────────────────────────────────────────
Entrez.email = "your_email@example.com"   # ← replace with your real email

# ─────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────
MIN_ORF_LENGTH    = 50   # minimum amino acids to keep an ORF
MAX_ORFS_TO_BLAST = 20   # raised from 10 → 20 so ROP isn't cut off
NCBI_DELAY_SEC    = 3    # seconds between BLAST calls (NCBI policy)

# NCBI allows ~3 requests/sec for registered users, ~1/sec otherwise.
# MAX_WORKERS=3 + the semaphore delay keeps us well within limits.
# If you hit HTTP 429s, lower MAX_WORKERS to 2.
MAX_WORKERS       = 3

# ─────────────────────────────────────────────
# RATE-LIMIT SEMAPHORE
# Controls concurrent access so we never fire more than MAX_WORKERS
# requests at once, and each thread waits NCBI_DELAY_SEC before calling.
# ─────────────────────────────────────────────
_ncbi_semaphore = threading.Semaphore(MAX_WORKERS)

# ─────────────────────────────────────────────
# CODON TABLE (DNA codons → Amino acid)
# ─────────────────────────────────────────────
CODON_TABLE = {
    # Phenylalanine
    'TTT': 'Phe (F)', 'TTC': 'Phe (F)',
    # Leucine
    'TTA': 'Leu (L)', 'TTG': 'Leu (L)',
    'CTT': 'Leu (L)', 'CTC': 'Leu (L)', 'CTA': 'Leu (L)', 'CTG': 'Leu (L)',
    # Isoleucine
    'ATT': 'Ile (I)', 'ATC': 'Ile (I)', 'ATA': 'Ile (I)',
    # Methionine (Start)
    'ATG': 'Met (M)',
    # Valine
    'GTT': 'Val (V)', 'GTC': 'Val (V)', 'GTA': 'Val (V)', 'GTG': 'Val (V)',
    # Serine
    'TCT': 'Ser (S)', 'TCC': 'Ser (S)', 'TCA': 'Ser (S)', 'TCG': 'Ser (S)',
    'AGT': 'Ser (S)', 'AGC': 'Ser (S)',
    # Proline
    'CCT': 'Pro (P)', 'CCC': 'Pro (P)', 'CCA': 'Pro (P)', 'CCG': 'Pro (P)',
    # Threonine
    'ACT': 'Thr (T)', 'ACC': 'Thr (T)', 'ACA': 'Thr (T)', 'ACG': 'Thr (T)',
    # Alanine
    'GCT': 'Ala (A)', 'GCC': 'Ala (A)', 'GCA': 'Ala (A)', 'GCG': 'Ala (A)',
    # Tyrosine
    'TAT': 'Tyr (Y)', 'TAC': 'Tyr (Y)',
    # STOP codons
    'TAA': '*** STOP ', 'TAG': ' STOP ', 'TGA': ' STOP ***',
    # Histidine
    'CAT': 'His (H)', 'CAC': 'His (H)',
    # Glutamine
    'CAA': 'Gln (Q)', 'CAG': 'Gln (Q)',
    # Asparagine
    'AAT': 'Asn (N)', 'AAC': 'Asn (N)',
    # Lysine
    'AAA': 'Lys (K)', 'AAG': 'Lys (K)',
    # Aspartic acid
    'GAT': 'Asp (D)', 'GAC': 'Asp (D)',
    # Glutamic acid
    'GAA': 'Glu (E)', 'GAG': 'Glu (E)',
    # Cysteine
    'TGT': 'Cys (C)', 'TGC': 'Cys (C)',
    # Tryptophan
    'TGG': 'Trp (W)',
    # Arginine
    'CGT': 'Arg (R)', 'CGC': 'Arg (R)', 'CGA': 'Arg (R)', 'CGG': 'Arg (R)',
    'AGA': 'Arg (R)', 'AGG': 'Arg (R)',
    # Glycine
    'GGT': 'Gly (G)', 'GGC': 'Gly (G)', 'GGA': 'Gly (G)', 'GGG': 'Gly (G)',
}

SHORT_NAME = {
    'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
    'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'Stop', 'TAG': 'Stop',
    'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'Stop', 'TGG': 'Trp',
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}

ONE_LETTER = {
    'Phe': 'F', 'Leu': 'L', 'Ile': 'I', 'Met': 'M', 'Val': 'V',
    'Ser': 'S', 'Pro': 'P', 'Thr': 'T', 'Ala': 'A', 'Tyr': 'Y',
    'His': 'H', 'Gln': 'Q', 'Asn': 'N', 'Lys': 'K', 'Asp': 'D',
    'Glu': 'E', 'Cys': 'C', 'Trp': 'W', 'Arg': 'R', 'Gly': 'G',
    'Stop': '*',
}

# ─────────────────────────────────────────────
# DNA SEQUENCE
# ─────────────────────────────────────────────
RAW_DNA = """
TTCTCATGTTTGACAGCTTATCATCGATAAGCTTTAATGCGGTAGTTTATCACAGTTAAATTGCTAACGC
AGTCAGGCACCGTGTATGAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTCACCCTGGATGCTG
TAGGCATAGGCTTGGTTATGCCGGTACTGCCGGGCCTCTTGCGGGATATCGTCCATTCCGACAGCATCGC
CAGTCACTATGGCGTGCTGCTAGCGCTATATGCGTTGATGCAATTTCTATGCGCACCCGTTCTCGGAGCA
CTGTCCGACCGCTTTGGCCGCCGCCCAGTCCTGCTCGCTTCGCTACTTGGAGCCACTATCGACTACGCGA
TCATGGCGACCACACCCGTCCTGTGGATCCTCTACGCCGGACGCATCGTGGCCGGCATCACCGGCGCCAC
AGGTGCGGTTGCTGGCGCCTATATCGCCGACATCACCGATGGGGAAGATCGGGCTCGCCACTTCGGGCTC
ATGAGCGCTTGTTTCGGCGTGGGTATGGTGGCAGGCCCCGTGGCCGGGGGACTGTTGGGCGCCATCTCCT
TGCATGCACCATTCCTTGCGGCGGCGGTGCTCAACGGCCTCAACCTACTACTGGGCTGCTTCCTAATGCA
GGAGTCGCATAAGGGAGAGCGTCGACCGATGCCCTTGAGAGCCTTCAACCCAGTCAGCTCCTTCCGGTGG
GCGCGGGGCATGACTATCGTCGCCGCACTTATGACTGTCTTCTTTATCATGCAACTCGTAGGACAGGTGC
CGGCAGCGCTCTGGGTCATTTTCGGCGAGGACCGCTTTCGCTGGAGCGCGACGATGATCGGCCTGTCGCT
TGCGGTATTCGGAATCTTGCACGCCCTCGCTCAAGCCTTCGTCACTGGTCCCGCCACCAAACGTTTCGGC
GAGAAGCAGGCCATTATCGCCGGCATGGCGGCCGACGCGCTGGGCTACGTCTTGCTGGCGTTCGCGACGC
GAGGCTGGATGGCCTTCCCCATTATGATTCTTCTCGCTTCCGGCGGCATCGGGATGCCCGCGTTGCAGGC
CATGCTGTCCAGGCAGGTAGATGACGACCATCAGGGACAGCTTCAAGGATCGCTCGCGGCTCTTACCAGC
CTAACTTCGATCACTGGACCGCTGATCGTCACGGCGATTTATGCCGCCTCGGCGAGCACATGGAACGGGT
TGGCATGGATTGTAGGCGCCGCCCTATACCTTGTCTGCCTCCCCGCGTTGCGTCGCGGTGCATGGAGCCG
GGCCACCTCGACCTGAATGGAAGCCGGCGGCACCTCGCTAACGGATTCACCACTCCAAGAATTGGAGCCA
ATCAATTCTTGCGGAGAACTGTGAATGCGCAAACCAACCCTTGGCAGAACATATCCATCGCGTCCGCCAT
CTCCAGCAGCCGCACGCGGCGCATCTCGGGCAGCGTTGGGTCCTGGCCACGGGTGCGCATGATCGTGCTC
CTGTCGTTGAGGACCCGGCTAGGCTGGCGGGGTTGCCTTACTGGTTAGCAGAATGAATCACCGATACGCG
AGCGAACGTGAAGCGACTGCTGCTGCAAAACGTCTGCGACCTGAGCAACAACATGAATGGTCTTCGGTTT
CCGTGTTTCGTAAAGTCTGGAAACGCGGAAGTCAGCGCCCTGCACCATTATGTTCCGGATCTGCATCGCA
GGATGCTGCTGGCTACCCTGTGGAACACCTACATCTGTATTAACGAAGCGCTGGCATTGACCCTGAGTGA
TTTTTCTCTGGTCCCGCCGCATCCATACCGCCAGTTGTTTACCCTCACAACGTTCCAGTAACCGGGCATG
TTCATCATCAGTAACCCGTATCGTGAGCATCCTCTCTCGTTTCATCGGTATCATTACCCCCATGAACAGA
AATCCCCCTTACACGGAGGCATCAGTGACCAAACAGGAAAAAACCGCCCTTAACATGGCCCGCTTTATCA
GAAGCCAGACATTAACGCTTCTGGAGAAACTCAACGAGCTGGACGCGGATGAACAGGCAGACATCTGTGA
ATCGCTTCACGACCACGCTGATGAGCTTTACCGCAGCTGCCTCGCGCGTTTCGGTGATGACGGTGAAAAC
CTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCC
GTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCGCAGCCATGACCCAGTCACGTAGCGATAGCGG
AGTGTATACTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAA
TACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCT
GCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAA
TCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCG
CGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAG
GTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCT
GTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATA
GCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCC
CGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTA
TCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCT
TGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGT
TACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTT
GTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGT
CTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCAC
CTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGAC
AGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCT
GACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACC
GCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGA
AGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTT
CGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTGCAGGCATCGTGGTGTCACGCTCGTCGTTTGG
TATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAA
GCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTA
TGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTC
AACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAACACGGGATAAT
ACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAA
GGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTT
TACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCG
ACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTC
TCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCG
AAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACG
AGGCCCTTTCGTCTTCAAGAA
"""

# ─────────────────────────────────────────────
# FUNCTIONS
# ─────────────────────────────────────────────

def clean_sequence(raw: str) -> str:
    """Strip whitespace and uppercase the DNA string."""
    return raw.replace('\n', '').replace(' ', '').upper()


def get_reverse_complement(dna: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(dna))


def translate_orf(dna: str, start: int) -> Tuple[List[str], bool]:
    """
    Translate a DNA sequence from 'start' into a list of amino acid short names.

    Returns
    -------
    protein   : list of three-letter amino acid names (e.g. ['Met', 'Ala', ...])
    has_stop  : True if translation ended with a stop codon (complete ORF),
                False if it ran off the end of the sequence (incomplete ORF)
    """
    protein: List[str] = []
    has_stop = False

    for i in range(start, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa = SHORT_NAME.get(codon, None)

        if aa is None:
            continue

        if aa == 'Stop':
            has_stop = True
            break

        protein.append(aa)

    return protein, has_stop


def find_orfs(
    dna: str,
    strand: str,
    original_length: int,
    min_length: int = MIN_ORF_LENGTH,
) -> List[Dict]:
    """
    Scan all three reading frames of dna and return ORFs that are at least
    min_length amino acids long AND have a confirmed stop codon.
    Deduplicates by stop-codon position, keeping only the longest ORF per stop.
    """
    orfs: List[Dict] = []

    for frame_idx in range(3):
        best_per_stop: Dict[int, Dict] = {}

        i = frame_idx
        while i <= len(dna) - 3:
            if dna[i:i + 3] == 'ATG':
                protein, has_stop = translate_orf(dna, i)

                if not has_stop:
                    i += 3
                    continue

                if len(protein) >= min_length:
                    stop_pos = i + len(protein) * 3

                    one_letter    = ''.join(ONE_LETTER.get(aa, '') for aa in protein)
                    display_pos   = (i + 1) if strand == '+' else (original_length - i)
                    display_frame = (frame_idx + 1) if strand == '+' else -(frame_idx + 1)

                    candidate = {
                        'frame'     : display_frame,
                        'start'     : display_pos,
                        'length'    : len(protein),
                        'protein'   : protein,
                        'one_letter': one_letter,
                        'stop_pos'  : stop_pos,
                    }

                    if (stop_pos not in best_per_stop or
                            len(protein) > best_per_stop[stop_pos]['length']):
                        best_per_stop[stop_pos] = candidate

            i += 3

        orfs.extend(best_per_stop.values())

    return orfs


def blast_protein(one_letter_seq: str, orf_id: int) -> Optional[Dict]:
    """
    Submit a protein sequence to BLASTp (SwissProt) and return the best hit.

    Thread-safe: uses _ncbi_semaphore to cap concurrent requests at MAX_WORKERS.
    Each thread sleeps NCBI_DELAY_SEC *while holding the semaphore* so the next
    thread can only acquire it after the delay, keeping requests spread out.
    """
    with _ncbi_semaphore:
        print(f"  🔍 BLASTing ORF {orf_id} ({len(one_letter_seq)} aa) "
              f"— waiting {NCBI_DELAY_SEC}s (NCBI policy)...")
        time.sleep(NCBI_DELAY_SEC)

        result_handle = None
        try:
            result_handle = NCBIWWW.qblast(
                "blastp", "swissprot", one_letter_seq,
                hitlist_size=1,
            )

            records = list(NCBIXML.parse(result_handle))

            if not records:
                print(f"     ⚠️  ORF {orf_id}: BLAST returned an empty response.")
                return None

            blast_record = records[0]

            if blast_record.alignments:
                best_algn = blast_record.alignments[0]
                best_hsp  = best_algn.hsps[0]
                return {
                    'orf_id'  : orf_id,
                    'name'    : best_algn.title.split(' >')[0],
                    'identity': round((best_hsp.identities / best_hsp.align_length) * 100, 1),
                    'e_value' : f"{best_hsp.expect:.2e}",
                    'score'   : best_hsp.score,
                }

        except Exception as exc:
            print(f"     ❌ ORF {orf_id} BLAST error: {exc}")

        finally:
            if result_handle is not None:
                result_handle.close()

    return None


# ─────────────────────────────────────────────
# MAIN EXECUTION
# ─────────────────────────────────────────────

def main() -> None:
    dna_fwd = clean_sequence(RAW_DNA)
    dna_rev = get_reverse_complement(dna_fwd)
    seq_len = len(dna_fwd)

    print(f"--- Plasmid Analysis Started (Length: {seq_len} bp) ---")

    # ── Identify ORFs on both strands ─────────────────────────────────────
    all_orfs: List[Dict] = []
    all_orfs.extend(find_orfs(dna_fwd, '+', seq_len))
    all_orfs.extend(find_orfs(dna_rev, '-', seq_len))

    # Sort longest → shortest
    all_orfs.sort(key=lambda x: x['length'], reverse=True)

    print(f"Detected {len(all_orfs)} unique protein-coding regions "
          f"(deduplicated by stop codon, confirmed stop required).")

    # ── Debug: print every ORF found before the BLAST cap ─────────────────
    print(f"\nAll ORFs (pre-BLAST):")
    for i, orf in enumerate(all_orfs, 1):
        print(f"  {i:>3}. Frame {orf['frame']:>3}  "
              f"Start {orf['start']:>6}  "
              f"Length {orf['length']:>4} aa  "
              f"Seq: {orf['one_letter'][:25]}...")
    print()

    # Apply the cap AFTER sorting so we always BLAST the longest ORFs
    orfs_to_blast = all_orfs[:MAX_ORFS_TO_BLAST]
    skipped = len(all_orfs) - len(orfs_to_blast)
    if skipped:
        print(f"⚠️  Only the top {MAX_ORFS_TO_BLAST} ORFs will be BLASTed "
              f"({skipped} shorter ORFs skipped).\n")

    # ── BLAST all ORFs concurrently ───────────────────────────────────────
    print(f"🚀 Launching {len(orfs_to_blast)} BLAST jobs "
          f"with up to {MAX_WORKERS} concurrent workers...\n")

    # future → (idx, orf) so we can reconstruct order after completion
    future_to_meta: Dict = {}
    best_results_summary: List[Optional[Dict]] = [None] * len(orfs_to_blast)

    start_time = time.time()

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for idx, orf in enumerate(orfs_to_blast, 1):
            future = executor.submit(blast_protein, orf['one_letter'], idx)
            future_to_meta[future] = (idx, orf)

        for future in as_completed(future_to_meta):
            idx, orf = future_to_meta[future]
            try:
                hit = future.result()
            except Exception as exc:
                print(f"     ❌ ORF #{idx} raised an unexpected exception: {exc}")
                hit = None

            if hit:
                hit['frame'] = orf['frame']
                hit['start'] = orf['start']
                best_results_summary[idx - 1] = hit
                print(f"  ✅ ORF #{idx} MATCH: {hit['name'][:60]}...")
                print(f"     Identity: {hit['identity']}% | E-value: {hit['e_value']}")
            else:
                print(f"  ❌ ORF #{idx}: No significant matches found.")

    elapsed = time.time() - start_time
    print(f"\n⏱  Total BLAST time: {elapsed:.1f}s  "
          f"(vs ~{len(orfs_to_blast) * NCBI_DELAY_SEC}s sequential minimum)")

    # Remove None entries (ORFs with no hit) and preserve sorted order
    best_results_summary = [r for r in best_results_summary if r is not None]

    # ── Final summary table ───────────────────────────────────────────────
    print("\n" + "█" * 80)
    print("  BLAST BEST RESULTS ARRAY (SUMMARY)")
    print("█" * 80)

    if not best_results_summary:
        print("No matches were found to store in the array.")
    else:
        header = f"{'ORF':<4} {'Frame':<6} {'Start':<8} {'Identity':<10} {'Protein Name'}"
        print(header)
        print("-" * 80)
        for entry in best_results_summary:
            print(
                f"{entry['orf_id']:<4} {entry['frame']:<6} {entry['start']:<8} "
                f"{entry['identity']:<9}% {entry['name'][:50]}"
            )

        print("\n--- Raw Data Array (list of dicts) ---")
        pprint.pprint(best_results_summary)


if __name__ == "__main__":
    main()