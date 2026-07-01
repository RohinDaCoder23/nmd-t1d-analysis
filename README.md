# NMD Analysis Pipeline — NOD Mouse T1D Study

Python pipeline for classifying mouse transcripts as NMD-sensitive or NMD-insensitive using the 50-nucleotide rule. Built to investigate whether Nonsense-Mediated Decay is inhibited during Type 1 Diabetes development in the NOD mouse model.

Developed for the Sussel Lab at Columbia University.

---

## Background

Nonsense-Mediated Decay (NMD) is a cellular quality control system that detects and destroys faulty messenger RNA transcripts before they can be translated into harmful proteins. The primary rule governing NMD sensitivity is the **50-nucleotide rule**: if a stop codon appears more than 50 nucleotides upstream of the last exon-exon junction, the transcript is flagged for destruction.

This project tests the hypothesis that pro-inflammatory cytokines associated with Type 1 Diabetes suppress NMD in pancreatic beta cells, allowing faulty transcripts to accumulate. To test this, the pipeline compares NMD-sensitive transcript prevalence between:

- **NOD mice** — develop Type 1 Diabetes naturally (disease group)
- **NOD scid mice** — similar genetic background but lack a working immune system (control group)

Data type: Long-read single-cell sequencing of pancreatic islets.

---

## What the Pipeline Does

For every transcript in a provided GTF + FASTA dataset, the pipeline performs four steps:

1. Locates the start codon and converts it from genomic to transcript coordinates
2. Performs in silico translation using BioPython
3. Identifies the first stop codon in transcript coordinates
4. Applies the 50nt NMD rule by comparing the stop codon position to the last exon-exon junction

Each transcript is classified as **NMD-sensitive** or **NMD-insensitive**, results are saved to CSV, and the classifications are validated against Ensembl's own NMD annotations.

---

## Validation

Validated against Ensembl NMD annotations on the mouse reference genome (GRCm38.96) across 47,534 protein-coding transcripts:

| Metric | Result |
|---|---|
| Overall agreement with Ensembl | 96.2% |
| Recall of known NMD targets | 98.5% |
| Precision | 74.2% |
| NMD-sensitive transcripts found | 6,767 (14.2%) |
| NMD-insensitive transcripts found | 40,767 (85.8%) |

---

## Requirements

### Dependencies
Install required Python libraries:

```bash
pip install -r requirements.txt
```

Or manually:

```bash
pip install gffutils biopython pandas
```

### Input Files
You must provide a matching GTF + FASTA pair and place both in the same folder as the script:

| File | Description |
|---|---|
| GTF file | Genome/transcript annotation (exons, start codons, biotypes) |
| FASTA file | Matching transcript sequences |

These files are not included in this repository due to their size. Mouse reference genome files can be downloaded from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html). The GTF and FASTA must come from the same source/assembly so their transcript IDs match.

---

## Usage

**1. Clone the repository**
```bash
git clone https://github.com/YOUR_USERNAME/nmd-t1d-analysis.git
cd nmd-t1d-analysis
```

**2. Install dependencies**
```bash
pip install -r requirements.txt
```

**3. Place your GTF and FASTA files in the project folder**

**4. Open `nmd_analysis.py` and edit the SETTINGS block at the top**

```python
# SETTINGS
GTF_FILE          = "Mus_musculus.GRCm38.96.gtf"   # input GTF annotation
FASTA_FILE        = "transcrips.transcripts.fa"     # input transcript FASTA
DB_FILE           = "mouse_genome.db"               # gffutils database (auto-built)

OUTPUT_FULL       = "nmd_results_full.csv"          # all classifications
OUTPUT_VALIDATED  = "nmd_results_validated.csv"     # classifications + Ensembl comparison

NMD_THRESHOLD     = 50    # nt upstream of last junction defining NMD-sensitive
ATG_SEARCH_WINDOW = 10    # nt either side of the estimated start to search for ATG
```

To run on a new dataset, change `GTF_FILE`, `FASTA_FILE`, the output names, and — importantly — `DB_FILE` to a new name. If `DB_FILE` already exists, the script loads that cached database instead of rebuilding from the new GTF, so give each dataset its own database filename.

**5. Run the pipeline**
```bash
python nmd_analysis.py
```

Runtime is approximately 1–2 minutes for ~140,000 transcripts (plus a one-time 5–10 minute database build the first time a new GTF is used).

---

## How It Works (Code Overview)

The script is organized into clearly separated sections:

| Section | Purpose |
|---|---|
| **Settings** | All filenames and parameters in one place — the only part you edit per dataset |
| **Data loading** | `load_database()` builds/loads the gffutils database from the GTF; `load_fasta()` loads transcript sequences |
| **Core functions** | Coordinate conversion, true-ATG finder, exon mapping, junction detection, and the NMD classification rule |
| **Pre-flight checks** | `run_accuracy_checks()` tests a small sample and blocks the full run if critical failures are found |
| **Full analysis** | `run_full_analysis()` classifies every transcript; `print_summary()` reports results and sanity checks |
| **Ensembl validation** | `validate_against_ensembl()` cross-checks calls against Ensembl biotype labels and reports precision/recall |

---

## Output

Results are saved to two CSV files.

`nmd_results_full.csv` — all classifications:

| Column | Description |
|---|---|
| `transcript_id` | Transcript ID |
| `chromosome` | Chromosome |
| `strand` | + or - (reading direction) |
| `transcript_length_nt` | Total transcript length in nucleotides |
| `num_exons` | Number of exons |
| `num_junctions` | Number of exon-exon junctions |
| `start_codon_transcript` | Start codon position in transcript coordinates |
| `stop_codon_transcript` | Stop codon position in transcript coordinates |
| `stop_codon_exon` | Which exon the stop codon falls in |
| `last_junction_pos` | Position of the last exon-exon junction |
| `distance_stop_to_junction` | Distance from stop codon to last junction (negative = after junction = normal) |
| `protein_length_aa` | Translated protein length in amino acids |
| `nmd_status` | NMD-sensitive or NMD-insensitive |

`nmd_results_validated.csv` — the above plus Ensembl comparison columns (`ensembl_biotype`, `ensembl_nmd_label`, `pipeline_vs_ensembl`).

---

## Limitations

The pipeline implements the primary 50nt NMD rule only. The following NMD mechanisms are not currently implemented:

- Upstream open reading frames (uORFs)
- Retained introns
- Long 3' UTR rule

The NMD-sensitive transcripts identified are therefore a conservative undercount of true NMD targets.

---

## Project Status

- [x] Pipeline built, tested, and validated on mouse reference genome
- [ ] NOD scid mouse islet data — partially received
- [ ] NOD mouse islet data — pending
- [ ] Statistical comparison of NMD-sensitive transcript prevalence between NOD and NOD scid

---

## References

- Katsioudi et al. 2023 — *Science Advances* — doi: 10.1126/sciadv.ade2828
- Karousis et al. 2021 — *Genome Biology* — PMC8361881
