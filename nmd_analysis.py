"""
NMD Analysis Pipeline

Classifies transcripts as NMD-sensitive or NMD-insensitive using the 50-nucleotide
rule, then validates the classifications against Ensembl's own NMD annotations.

For every transcript the pipeline:
  1. Finds the start codon and converts it to transcript coordinates
  2. Performs in silico translation with BioPython
  3. Finds the first stop codon in transcript coordinates
  4. Compares the stop codon to the last exon-exon junction and applies the 50nt rule

Edit the SETTINGS block below to point at your own GTF and FASTA files, then run:
    python nmd_analysis.py
"""

import os
import time

import gffutils
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

# ──────────────────────────────────────────────────────────────────────────────
# SETTINGS
# ──────────────────────────────────────────────────────────────────────────────

GTF_FILE          = "Mus_musculus.GRCm38.96.gtf"   # input GTF annotation
FASTA_FILE        = "transcrips.transcripts.fa"     # input transcript FASTA
DB_FILE           = "mouse_genome.db"               # gffutils database (auto-built)

OUTPUT_FULL       = "nmd_results_full.csv"          # all classifications
OUTPUT_VALIDATED  = "nmd_results_validated.csv"     # classifications + Ensembl comparison

NMD_THRESHOLD     = 50    # nt upstream of last junction defining NMD-sensitive
ATG_SEARCH_WINDOW = 10    # nt either side of the estimated start to search for ATG
PROGRESS_INTERVAL = 5000  # print progress every N transcripts
PREFLIGHT_SAMPLE  = 10    # transcripts to test before the full run

STOP_CODONS = ("TAA", "TAG", "TGA")


# ──────────────────────────────────────────────────────────────────────────────
# DATA LOADING
# ──────────────────────────────────────────────────────────────────────────────

def load_database(gtf_file, db_file):
    """Load the gffutils database, building it from the GTF on first run."""
    if os.path.exists(db_file):
        print(f"Loading database from {db_file}")
        return gffutils.FeatureDB(db_file)

    print(f"Building database from {gtf_file} (first run, 5-10 minutes)...")
    db = gffutils.create_db(
        gtf_file,
        dbfn=db_file,
        force=False,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
    )
    print("Database built and saved.")
    return db


def load_fasta(fasta_file):
    """Load all transcript sequences into a dict keyed by transcript ID."""
    print(f"Loading sequences from {fasta_file}")
    return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))


# ──────────────────────────────────────────────────────────────────────────────
# CORE FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def get_exons(db, transcript_id):
    """Return all exons for a transcript, sorted by genomic start."""
    try:
        return list(db.children(transcript_id, featuretype="exon", order_by="start"))
    except Exception:
        return []


def get_start_codon_genomic(db, transcript_id):
    """Return the genomic start position of the start codon, or None if absent."""
    for sc in db.children(transcript_id, featuretype="start_codon"):
        return sc.start
    return None


def genomic_to_transcript(genomic_pos, exons, strand):
    """
    Convert a genomic coordinate to a transcript coordinate.

    Genomic coordinates run along the whole chromosome; transcript coordinates run
    from 0 within the spliced transcript. Negative-strand transcripts are read in
    reverse, so the exon order is flipped before walking through them.
    """
    transcript_pos = 0

    if strand == "-":
        exons = sorted(exons, key=lambda x: x.start, reverse=True)
        for exon in exons:
            if exon.start <= genomic_pos <= exon.end:
                return transcript_pos + (exon.end - genomic_pos)
            transcript_pos += (exon.end - exon.start)
    else:
        for exon in exons:
            if exon.start <= genomic_pos <= exon.end:
                return transcript_pos + (genomic_pos - exon.start)
            transcript_pos += (exon.end - exon.start)

    return None


def find_true_atg(sequence, estimated_pos, search_window=ATG_SEARCH_WINDOW):
    """
    Find the actual ATG near an estimated start position.

    GTF files are 1-based, Python is 0-based, and strand orientation adds small
    shifts, so the converted position is usually within a few nt of the true ATG
    but not exact. Searching a small window avoids relying on a hardcoded offset.
    """
    search_start = max(0, estimated_pos - search_window)
    search_end = min(len(sequence) - 3, estimated_pos + search_window)

    for pos in range(search_start, search_end + 1):
        if str(sequence[pos:pos + 3]).upper() == "ATG":
            return pos
    return None


def get_exon_coordinate_map(exons, strand):
    """Map every exon onto transcript coordinates (transcript_start/transcript_end)."""
    exon_coords = []
    cumulative = 0

    if strand == "-":
        exons = sorted(exons, key=lambda x: x.start, reverse=True)

    for i, exon in enumerate(exons):
        exon_length = exon.end - exon.start
        exon_coords.append({
            "exon_number": i + 1,
            "genomic_start": exon.start,
            "genomic_end": exon.end,
            "transcript_start": cumulative,
            "transcript_end": cumulative + exon_length,
            "length": exon_length,
        })
        cumulative += exon_length

    return exon_coords


def get_junction_positions(exon_coords):
    """Return exon-exon junction positions (always len(exons) - 1 of them)."""
    return [e["transcript_end"] for e in exon_coords[:-1]]


def classify_nmd(stop_pos, junctions, threshold=NMD_THRESHOLD):
    """
    Apply the 50nt NMD rule.

    distance = last_junction - stop_pos
      positive  -> stop codon is BEFORE the last junction (NMD risk)
      negative  -> stop codon is AFTER the last junction (normal)
    A transcript is NMD-sensitive when distance exceeds the threshold.
    """
    if not junctions:
        return "NMD-insensitive", None, None

    last_junction = junctions[-1]
    distance = last_junction - stop_pos
    status = "NMD-sensitive" if distance > threshold else "NMD-insensitive"
    return status, last_junction, distance


def analyze_transcript(db, fasta_sequences, transcript):
    """
    Run all four steps on a single transcript.

    Returns a result dict, or a string skip-reason if the transcript can't be
    classified (no FASTA entry, no start codon, too few exons, no ATG, no stop).
    """
    transcript_id = transcript.id
    strand = transcript.strand

    if transcript_id not in fasta_sequences:
        return "no_fasta"

    start_codon_genomic = get_start_codon_genomic(db, transcript_id)
    if start_codon_genomic is None:
        return "no_start"

    exons = get_exons(db, transcript_id)
    if len(exons) < 2:
        return "no_exons"

    # Step 1: start codon in transcript coordinates
    estimated_start = genomic_to_transcript(start_codon_genomic, exons, strand)
    if estimated_start is None:
        return "no_coord"

    sequence = fasta_sequences[transcript_id].seq
    true_start = find_true_atg(sequence, estimated_start)
    if true_start is None:
        return "no_atg"

    # Step 2: in silico translation
    coding_sequence = Seq(str(sequence[true_start:]).upper())
    translated = coding_sequence.translate(to_stop=False)

    # Step 3: first stop codon in transcript coordinates
    stop_in_protein = str(translated).find("*")
    if stop_in_protein == -1:
        return "no_stop"
    stop_pos = true_start + (stop_in_protein * 3)

    # Step 4: exon map, junctions, NMD classification
    exon_coords = get_exon_coordinate_map(exons, strand)
    junctions = get_junction_positions(exon_coords)
    nmd_status, last_junction, distance = classify_nmd(stop_pos, junctions)

    stop_exon = None
    for e in exon_coords:
        if e["transcript_start"] <= stop_pos < e["transcript_end"]:
            stop_exon = e["exon_number"]
            break

    return {
        "transcript_id": transcript_id,
        "chromosome": transcript.chrom,
        "strand": strand,
        "transcript_length_nt": len(sequence),
        "num_exons": len(exons),
        "num_junctions": len(junctions),
        "start_codon_transcript": true_start,
        "stop_codon_transcript": stop_pos,
        "stop_codon_exon": stop_exon,
        "last_junction_pos": last_junction,
        "distance_stop_to_junction": distance,
        "protein_length_aa": stop_in_protein,
        "nmd_status": nmd_status,
    }


# ──────────────────────────────────────────────────────────────────────────────
# PRE-FLIGHT ACCURACY CHECKS
# ──────────────────────────────────────────────────────────────────────────────

def run_accuracy_checks(db, fasta_sequences, sample_size=PREFLIGHT_SAMPLE):
    """
    Verify the pipeline on a small sample before the full run.

    Critical failures (ATG not found, invalid stop codon, wrong junction count)
    indicate a code or coordinate problem and block the full analysis. Transcripts
    with no stop codon are incomplete sequences, not errors, and are reported but
    do not block the run.
    """
    print(f"Running accuracy checks on {sample_size} transcripts before full analysis\n")

    tested = 0
    atg_confirmed = 0
    valid_stops = 0
    incomplete = 0
    junction_ok = 0
    critical_failures = 0

    for transcript in db.features_of_type("transcript"):
        if tested >= sample_size:
            break

        tid = transcript.id
        strand = transcript.strand

        if tid not in fasta_sequences:
            continue

        start_codon_genomic = get_start_codon_genomic(db, tid)
        if start_codon_genomic is None:
            continue

        exons = get_exons(db, tid)
        if len(exons) < 2:
            continue

        estimated_start = genomic_to_transcript(start_codon_genomic, exons, strand)
        if estimated_start is None:
            continue

        sequence = fasta_sequences[tid].seq

        # Check 1: ATG found within the search window
        true_start = find_true_atg(sequence, estimated_start)
        if true_start is None:
            critical_failures += 1
            print(f"  CRITICAL: ATG not found for {tid} (estimated {estimated_start})")
            tested += 1
            continue
        atg_confirmed += 1

        # Check 2: valid stop codon
        translated = Seq(str(sequence[true_start:]).upper()).translate(to_stop=False)
        stop_in_protein = str(translated).find("*")
        if stop_in_protein == -1:
            incomplete += 1
            print(f"  Incomplete transcript {tid}: no stop codon (will be skipped)")
            tested += 1
            continue

        stop_pos = true_start + (stop_in_protein * 3)
        stop_letters = str(sequence[stop_pos:stop_pos + 3]).upper()
        if stop_letters in STOP_CODONS:
            valid_stops += 1
        else:
            critical_failures += 1
            print(f"  CRITICAL: invalid stop codon for {tid} (found {stop_letters})")

        # Check 3: junction count equals exons - 1
        exon_coords = get_exon_coordinate_map(exons, strand)
        junctions = get_junction_positions(exon_coords)
        if len(junctions) == len(exons) - 1:
            junction_ok += 1
        else:
            critical_failures += 1
            print(f"  CRITICAL: junction count wrong for {tid} "
                  f"(expected {len(exons) - 1}, got {len(junctions)})")

        tested += 1

    complete_tested = tested - incomplete

    print("\nAccuracy check summary")
    print(f"  Transcripts tested:         {tested}")
    print(f"  Complete transcripts:       {complete_tested}")
    print(f"  Incomplete (no stop codon): {incomplete}")
    print(f"  ATG found within window:    {atg_confirmed}/{tested}")
    print(f"  Valid stop codons:          {valid_stops}/{complete_tested}")
    print(f"  Junction counts correct:    {junction_ok}/{complete_tested}")
    print(f"  Critical failures:          {critical_failures}")

    if critical_failures == 0:
        print("\nAll critical checks passed. Starting full analysis.\n")
        return True

    print("\nCritical failures found. Full analysis will not run.")
    print("If ATG checks are failing, try increasing ATG_SEARCH_WINDOW.\n")
    return False


# ──────────────────────────────────────────────────────────────────────────────
# FULL ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────

def run_full_analysis(db, fasta_sequences):
    """Classify every transcript in the database. Returns (DataFrame, counts dict)."""
    results = []
    counts = {
        "total": 0, "no_fasta": 0, "no_start": 0, "no_exons": 0,
        "no_coord": 0, "no_atg": 0, "no_stop": 0, "errors": 0,
    }

    total_transcripts = sum(1 for _ in db.features_of_type("transcript"))
    start_time = time.time()

    print("Starting full analysis")
    print(f"  Transcripts to process: {total_transcripts:,}")
    print(f"  NMD threshold:          {NMD_THRESHOLD} nt")
    print(f"  ATG search window:      +/- {ATG_SEARCH_WINDOW} nt\n")

    for transcript in db.features_of_type("transcript"):
        counts["total"] += 1

        if counts["total"] % PROGRESS_INTERVAL == 0:
            pct = 100 * counts["total"] / total_transcripts
            print(f"  {pct:.1f}% complete ({counts['total']:,} / {total_transcripts:,})")

        try:
            result = analyze_transcript(db, fasta_sequences, transcript)
            if isinstance(result, dict):
                results.append(result)
            else:
                counts[result] += 1
        except Exception:
            counts["errors"] += 1

    counts["analyzed"] = len(results)
    counts["elapsed_minutes"] = (time.time() - start_time) / 60
    return pd.DataFrame(results), counts


def print_summary(results_df, counts):
    """Print classification counts, accuracy indicators, and sanity checks."""
    total_classified = len(results_df)
    sensitive = int((results_df["nmd_status"] == "NMD-sensitive").sum())
    insensitive = int((results_df["nmd_status"] == "NMD-insensitive").sum())

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)

    print("\nTranscript disposition")
    print(f"  Total in GTF:                   {counts['total']:,}")
    print(f"  Successfully analyzed:          {counts['analyzed']:,} "
          f"({100 * counts['analyzed'] / counts['total']:.1f}%)")
    print(f"  Skipped - not in FASTA:         {counts['no_fasta']:,}")
    print(f"  Skipped - no start codon:       {counts['no_start']:,}")
    print(f"  Skipped - fewer than 2 exons:   {counts['no_exons']:,}")
    print(f"  Skipped - no ATG in window:     {counts['no_atg']:,}")
    print(f"  Skipped - no stop codon:        {counts['no_stop']:,}")
    print(f"  Skipped - coordinate issues:    {counts['no_coord']:,}")
    print(f"  Skipped - unexpected errors:    {counts['errors']:,}")

    print("\nNMD classification results")
    print(f"  NMD-sensitive:    {sensitive:,} ({100 * sensitive / total_classified:.1f}%)")
    print(f"  NMD-insensitive:  {insensitive:,} ({100 * insensitive / total_classified:.1f}%)")
    print(f"  Total classified: {total_classified:,}")

    print("\nAccuracy indicators")
    atg_fail_rate = 100 * counts["no_atg"] / counts["total"]
    error_rate = 100 * counts["errors"] / counts["total"]
    print(f"  ATG not found rate: {atg_fail_rate:.2f}% "
          f"({'acceptable' if atg_fail_rate < 10 else 'high - raise ATG_SEARCH_WINDOW'})")
    print(f"  Error rate:         {error_rate:.2f}% "
          f"({'acceptable' if error_rate < 5 else 'high - review results'})")

    print("\nBiological sanity checks")
    avg_protein = results_df["protein_length_aa"].mean()
    avg_exons = results_df["num_exons"].mean()
    print(f"  Average protein length: {avg_protein:.0f} aa "
          f"({'reasonable' if 50 < avg_protein < 2000 else 'unusual - investigate'})")
    print(f"  Average exon count:     {avg_exons:.1f} per transcript")

    print(f"\nRuntime: {counts['elapsed_minutes']:.1f} minutes")


# ──────────────────────────────────────────────────────────────────────────────
# ENSEMBL VALIDATION
# ──────────────────────────────────────────────────────────────────────────────

def get_ensembl_labels(db):
    """
    Read Ensembl's own NMD calls from the GTF transcript biotype attribute.

    Returns (ensembl_nmd set, ensembl_biotypes dict). Transcripts whose biotype is
    'nonsense_mediated_decay' are Ensembl's NMD-sensitive calls.
    """
    ensembl_nmd = set()
    ensembl_biotypes = {}

    for transcript in db.features_of_type("transcript"):
        tid = transcript.id
        biotype = None
        for attr_name in ("transcript_biotype", "transcript_type", "biotype"):
            try:
                biotype = transcript.attributes[attr_name][0]
                break
            except (KeyError, IndexError):
                continue

        if biotype:
            ensembl_biotypes[tid] = biotype
            if biotype == "nonsense_mediated_decay":
                ensembl_nmd.add(tid)

    return ensembl_nmd, ensembl_biotypes


def validate_against_ensembl(results_df, db):
    """Compare pipeline calls to Ensembl labels and report agreement metrics."""
    print("\n" + "=" * 60)
    print("VALIDATION AGAINST ENSEMBL NMD ANNOTATIONS")
    print("=" * 60)

    ensembl_nmd, ensembl_biotypes = get_ensembl_labels(db)
    print(f"\n  Transcripts with biotype labels: {len(ensembl_biotypes):,}")
    print(f"  Ensembl NMD-labeled transcripts: {len(ensembl_nmd):,}")

    if not ensembl_biotypes:
        print("\n  No biotype labels found in GTF. Skipping validation.")
        return results_df

    true_pos = true_neg = false_pos = false_neg = no_label = 0

    for _, row in results_df.iterrows():
        tid = row["transcript_id"]
        our_call = row["nmd_status"]

        if tid not in ensembl_biotypes:
            no_label += 1
            continue

        ensembl_sensitive = tid in ensembl_nmd
        our_sensitive = our_call == "NMD-sensitive"

        if our_sensitive and ensembl_sensitive:
            true_pos += 1
        elif not our_sensitive and not ensembl_sensitive:
            true_neg += 1
        elif our_sensitive and not ensembl_sensitive:
            false_pos += 1
        else:
            false_neg += 1

    total_comparable = true_pos + true_neg + false_pos + false_neg
    if total_comparable == 0:
        print("\n  No transcripts had comparable Ensembl labels. Skipping metrics.")
        return results_df

    agreement = 100 * (true_pos + true_neg) / total_comparable
    precision = 100 * true_pos / (true_pos + false_pos) if (true_pos + false_pos) else 0
    recall = 100 * true_pos / (true_pos + false_neg) if (true_pos + false_neg) else 0

    print("\n  Confusion matrix (pipeline vs Ensembl)")
    print(f"    True Positive  (both NMD-sensitive):   {true_pos:,}")
    print(f"    True Negative  (both insensitive):     {true_neg:,}")
    print(f"    False Positive (we NMD, Ensembl no):   {false_pos:,}")
    print(f"    False Negative (Ensembl NMD, we no):   {false_neg:,}")
    print(f"    No Ensembl label:                      {no_label:,}")

    print("\n  Metrics")
    print(f"    Overall agreement: {agreement:.1f}%")
    print(f"    Precision:         {precision:.1f}%")
    print(f"    Recall:            {recall:.1f}%")

    # Add comparison columns and save
    results_df = results_df.copy()
    results_df["ensembl_biotype"] = results_df["transcript_id"].map(
        lambda x: ensembl_biotypes.get(x, "unknown")
    )
    results_df["ensembl_nmd_label"] = results_df["transcript_id"].map(
        lambda x: "NMD-sensitive" if x in ensembl_nmd
        else ("NMD-insensitive" if x in ensembl_biotypes else "no_label")
    )

    def label_outcome(row):
        ours = row["nmd_status"]
        ens = row["ensembl_nmd_label"]
        if ours == "NMD-sensitive" and ens == "NMD-sensitive":
            return "True_Positive"
        if ours == "NMD-insensitive" and ens == "NMD-insensitive":
            return "True_Negative"
        if ours == "NMD-sensitive" and ens == "NMD-insensitive":
            return "False_Positive"
        if ours == "NMD-insensitive" and ens == "NMD-sensitive":
            return "False_Negative"
        return "No_Ensembl_Label"

    results_df["pipeline_vs_ensembl"] = results_df.apply(label_outcome, axis=1)
    return results_df


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("NMD ANALYSIS PIPELINE")
    print("=" * 60 + "\n")

    db = load_database(GTF_FILE, DB_FILE)
    fasta_sequences = load_fasta(FASTA_FILE)

    transcript_count = sum(1 for _ in db.features_of_type("transcript"))
    print(f"  Transcripts in GTF:   {transcript_count:,}")
    print(f"  Sequences in FASTA:   {len(fasta_sequences):,}\n")

    if not run_accuracy_checks(db, fasta_sequences):
        return

    results_df, counts = run_full_analysis(db, fasta_sequences)
    print_summary(results_df, counts)

    results_df.to_csv(OUTPUT_FULL, index=False)
    print(f"\nSaved classifications to {OUTPUT_FULL}")

    validated_df = validate_against_ensembl(results_df, db)
    validated_df.to_csv(OUTPUT_VALIDATED, index=False)
    print(f"Saved validated results to {OUTPUT_VALIDATED}")


if __name__ == "__main__":
    main()
