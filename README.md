# NMD Analysis Pipeline — NOD Mouse T1D Study

Python pipeline for classifying mouse transcripts as NMD-sensitive or NMD-insensitive using the 50-nucleotide rule. Built to investigate whether Nonsense-Mediated Decay is inhibited during Type 1 Diabetes development in the NOD mouse model.


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

1. Locates the start codon in transcript coordinates using the GTF annotation
2. Performs in silico translation using BioPython
3. Identifies the first stop codon in transcript coordinates
4. Applies the 50nt NMD rule by comparing the stop codon position to the last exon-exon junction

Each transcript is classified as **NMD-sensitive** or **NMD-insensitive** and results are saved to a CSV file.

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
Install required Python libraries using:

```bash
pip install -r requirements.txt
```

Or manually:

```bash
pip install gffutils biopython pandas
```

### Input Files
You must provide two files and place them in the same folder as the script:

| File | Description |
|---|---|
| `your_file.gtf` | GTF genome annotation file for your organism/assembly |
| `your_file.fa` | Matching FASTA file of transcript sequences |

These files are not included in this repository due to their size. Mouse reference genome files can be downloaded from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html).

---

## Usage

**1. Clone the repository**
```bash
git clone https://github.com/YOUR_USERNAME/nmd-analysis-pipeline.git
cd nmd-analysis-pipeline
```

**2. Install dependencies**
```bash
pip install -r requirements.txt
```

**3. Place your GTF and FASTA files in the project folder**

**4. Open `nmd_analysis.py` and edit the settings block at the top of the file**
```python
# ── USER SETTINGS ──────────────────────────────
GTF_FILE      = "your_file.gtf"       # replace with your GTF filename
FASTA_FILE    = "your_file.fa"        # replace with your FASTA filename
OUTPUT_FILE   = "results/nmd_results_full.csv"
# ───────────────────────────────────────────────
```

**5. Run the pipeline**
```bash
python nmd_analysis.py
```

Runtime is approximately 1–2 minutes for ~140,000 transcripts.

---

## Output

Results are saved to `results/nmd_results_full.csv` with the following columns:

| Column | Description |
|---|---|
| `transcript_id` | Ensembl transcript ID |
| `chromosome` | Chromosome number |
| `strand` | + or - indicating reading direction |
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

---

## Limitations

The pipeline implements the primary 50nt NMD rule only. The following NMD mechanisms are not currently implemented:

- Upstream open reading frames (uORFs)
- Retained introns
- Long 3' UTR rule

---

## Project Status

- [x] Pipeline built, tested, and validated on mouse reference genome
- [ ] NOD mouse islet sequencing data — pending
- [ ] NOD scid mouse islet sequencing data — pending
- [ ] Statistical comparison of NMD-sensitive transcript prevalence between NOD and NOD scid

---

## References

- Katsioudi et al. 2023 — *Science Advances* — doi: 10.1126/sciadv.ade2828
- Karousis et al. 2021 — *Genome Biology* — PMC8361881
