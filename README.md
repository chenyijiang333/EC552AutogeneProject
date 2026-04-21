# EcoCyc + UniProt + STRING Circuit Viewer

This project builds gene/protein interaction graphs from:

- local EcoCyc PGDB data (`29.6`)
- UniProt resolution for genes not found in EcoCyc
- BLAST-based gene hints from pasted DNA / FASTA
- STRING DB interaction scoring for UniProt-resolved proteins

The web app is served by Flask from `circuit_viewer/server.py`.

## 1) Requirements

- Python 3.10+
- Local EcoCyc PGDB folder in this repo as `29.6` (already expected by default)

Python packages:

```bash
pip install flask flask-cors networkx biopython
```

- Install Blast+ Software:

Mac: brew install blast
Ubuntu/Linux: sudo apt-get install ncbi-blast+
Windows: Download and run the .exe from NCBI's FTP site --> https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

- Download the Swiss-Prot file from Uniprot --> https://www.uniprot.org/help/downloads
- Unzip the file in the working directory, it should have created swissprot.fasta
- While in the same directory as swissprot.fasta, run the command

```bash
makeblastdb -in swissprot.fasta -dbtype prot -out swissprot
````

## 2) Start the Web Server

From the project root:

```bash
python circuit_viewer/server.py
```

Open:

- [http://127.0.0.1:5000](http://127.0.0.1:5000)

Optional flags:

```bash
python circuit_viewer/server.py --pgdb-dir 29.6 --host 127.0.0.1 --port 5000
```

## 3) Using the Web Page

### A) Build graph from gene list

1. In the left sidebar, enter genes/IDs one per line (e.g. `lacI`, `EG10525`).
2. Set **STRING interaction score threshold** slider (0.00 to 1.00).
3. Click **Build graph**.

What happens:

- each token is resolved in EcoCyc first
- if not found, it resolves through UniProt
- graph is built and shown in the viewer
- for UniProt seeds, STRING edges are added only when score >= threshold

### B) Build from DNA / FASTA

1. Paste DNA / FASTA text or upload a FASTA file.
2. Click **Import sequence -> gene list**.
3. The app runs ORF detection and BLAST to suggest genes, then auto-builds the graph.

Notes are shown under the button (ORFs found, resolved seeds, BLAST notes).

## 4) BLAST Modes

The sequence pipeline in `Autogene.py` supports:

- **Local BLASTp** (faster/offline, preferred if installed)
- **NCBI web BLASTp** fallback (slower, requires internet)

### Local BLAST setup (optional but recommended)

Install NCBI BLAST+ and point to `blastp`:

- Add BLAST `bin` folder to `PATH`, or
- set `BLASTP_EXE` to full path of `blastp.exe`

Prepare SwissProt DB prefix:

```bash
makeblastdb -in swissprot/uniprot_sprot.fasta -dbtype prot -out swissprot/swissprot
```

Then set:

- `SWISSPROT_DB` (or `BLAST_DB`) to DB prefix path (no extension), for example:
  - `.../swissprot/swissprot`

## 5) Environment Variables

### Timeouts / performance

- `GRAPH_TIMEOUT` (default: `120`)  
  max seconds for `/api/graph`
- `SEQUENCE_ANALYSIS_TIMEOUT` (default: `600`)  
  max seconds for `/api/analyze-sequence`

### Remote BLAST controls

- `SEQ_REMOTE_BLAST` (default: enabled)  
  set `0`/`false`/`no` to disable web BLAST fallback
- `SEQ_REMOTE_ORF_CAP` (default: `5`)  
  max ORFs sent to web BLAST
- `SEQ_REMOTE_BLAST_DELAY_SEC` (default: `4`)  
  delay between web BLAST ORFs
- `SEQ_REMOTE_BLAST_TIMEOUT_SEC` (default: `180`)  
  timeout per web BLAST request
- `NCBI_EMAIL` (or `ENTREZ_EMAIL`)  
  optional email for NCBI requests

### STRING behavior

- threshold is sent from UI slider as `string_score_threshold` in `/api/graph`
- default if missing is `0.7`

## 6) Troubleshooting

- **No edges for UniProt-only seeds**  
  Lower STRING threshold slider (try `0.4`).

- **`blastp not found` note**  
  Install BLAST+ and set `BLASTP_EXE`, or rely on web BLAST fallback.

- **`Graph build timed out`**  
  Increase `GRAPH_TIMEOUT`, reduce input genes, or lower graph complexity.

- **`POST /api/graph` returns 200 but empty graph**  
  Means genes resolved but no qualifying interactions at current filters/threshold.

## 7) API Endpoints (used by UI)

- `GET /api/health`
- `POST /api/analyze-sequence`
  - input: JSON `{"sequence": "...", "max_orfs": 10}` or multipart `file`
  - output includes: `suggested_genes`, `orfs`, `notes`, `resolved_seeds`
- `POST /api/graph`
  - input: `{"genes": [...], "string_score_threshold": 0.7}`
  - output: `nodes`, `edges`, `resolved`, `errors`

