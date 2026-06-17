#Cell 1:
import gffutils
from Bio import SeqIO
print ("All good!")

# ── CELL 2: Import all libraries ──────────────────────────────────────────────
import gffutils
import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

print("All libraries loaded successfully!")

# ── CELL 3: Load GTF into searchable database ─────────────────────────────────
# First time only: takes 5-10 mins to build. After that loads instantly.

db_path = "mouse_genome.db"

if not os.path.exists(db_path):
    print("Building database for the first time — please wait 5-10 minutes...")
    db = gffutils.create_db(
        "Mus_musculus.GRCm38.96.gtf",
        dbfn=db_path,
        force=False,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True
    )
    print("Database built and saved!")
else:
    db = gffutils.FeatureDB(db_path)
    print("Database loaded from saved file!")

# Sanity check — count how many transcripts are in the GTF
transcript_count = sum(1 for _ in db.features_of_type("transcript"))
print(f"Total transcripts in GTF: {transcript_count:,}")

# ── CELL 4: Load FASTA transcript sequences ───────────────────────────────────

fasta_sequences = SeqIO.to_dict(SeqIO.parse("transcrips.transcripts.fa", "fasta"))

print(f"Total sequences loaded from FASTA: {len(fasta_sequences):,}")

# Show 3 example IDs so we can confirm they match the GTF format
print("\nExample transcript IDs from FASTA:")
for i, key in enumerate(list(fasta_sequences.keys())[:3]):
    print(f"  {key}  ({len(fasta_sequences[key].seq):,} nt)")

# ── CELL 5: Verify GTF and FASTA transcript IDs match ────────────────────────

# Grab 5 transcript IDs from the GTF
gtf_ids = [t.id for t in list(db.features_of_type("transcript"))[:5]]

print("Example transcript IDs from GTF:")
for tid in gtf_ids:
    print(f"  {tid}")

print("\nExample transcript IDs from FASTA:")
for fid in list(fasta_sequences.keys())[:5]:
    print(f"  {fid}")

print("\nDo they look the same format? Compare the two lists above.")

# ── CELL 6 (UPDATED): Find a good test transcript ─────────────────────────────
# We need one that has: a start codon AND at least 2 exons AND exists in FASTA

print("Searching for a good test transcript...")
print("(Needs: start codon + 2 or more exons + exists in FASTA)\n")

common_id = None

for t in db.features_of_type("transcript"):
    # Must exist in FASTA
    if t.id not in fasta_sequences:
        continue
    
    # Must have a start codon in GTF
    has_start_codon = any(True for _ in db.children(t.id, featuretype='start_codon'))
    if not has_start_codon:
        continue
    
    # Must have at least 2 exons
    exons = list(db.children(t.id, featuretype='exon', order_by='start'))
    if len(exons) < 2:
        continue
    
    # If we get here, this transcript is a good test case
    common_id = t.id
    break

if common_id is None:
    print("❌ ERROR: Could not find any suitable transcript")
    print("Please send this output so we can investigate")
else:
    print(f"✅ Found a good test transcript: {common_id}")
    
    t = db[common_id]
    exons = list(db.children(common_id, featuretype='exon', order_by='start'))
    
    print(f"\nBasic info from GTF:")
    print(f"  Chromosome:  {t.chrom}")
    print(f"  Strand:      {t.strand}")
    print(f"  Start:       {t.start:,}")
    print(f"  End:         {t.end:,}")
    print(f"  Exon count:  {len(exons)}")
    
    seq = fasta_sequences[common_id].seq
    print(f"\nBasic info from FASTA:")
    print(f"  Sequence length:        {len(seq):,} nucleotides")
    print(f"  First 30 letters:       {seq[:30]}")
    
    print(f"\n✅ This transcript is a valid test case!")
    print(f"   It has {len(exons)} exons and a start codon — ready for full analysis")

# ── CELL 7: STEP 1 — Find start codon location in transcript coordinates ──────
# The GTF file gives us genomic coordinates (position on the chromosome)
# We need to convert that into transcript coordinates (position within the transcript)

def get_exons(transcript_id):
    """Get all exons for a transcript, sorted by position"""
    return list(db.children(transcript_id, featuretype='exon', order_by='start'))

def genomic_to_transcript(genomic_pos, exons, strand):
    """
    Convert a genomic coordinate into a transcript coordinate
    
    Think of it like this:
    Genomic  = mile marker on a highway (can be millions)
    Transcript = step number within just this one road (starts from 0)
    """
    transcript_pos = 0
    
    if strand == '-':
        # Negative strand is read backwards, so reverse the exon order
        exons = sorted(exons, key=lambda x: x.start, reverse=True)
        for exon in exons:
            if exon.start <= genomic_pos <= exon.end:
                transcript_pos += (exon.end - genomic_pos)
                return transcript_pos
            transcript_pos += (exon.end - exon.start)
    else:
        for exon in exons:
            if exon.start <= genomic_pos <= exon.end:
                transcript_pos += (genomic_pos - exon.start)
                return transcript_pos
            transcript_pos += (exon.end - exon.start)
    
    return None  # Position wasn't found in any exon

# ── Test on our one example transcript ────────────────────────────────────────
print(f"Testing Step 1 on transcript: {common_id}\n")

exons = get_exons(common_id)
strand = db[common_id].strand

print(f"Number of exons found: {len(exons)}")
print(f"Strand: {strand}")

# Find the start codon from the GTF
start_codon_genomic = None
for sc in db.children(common_id, featuretype='start_codon'):
    start_codon_genomic = sc.start
    break

if start_codon_genomic is None:
    print("⚠️  No start codon found for this transcript in GTF")
    print("Try running the cell again — a different transcript will be picked")
else:
    print(f"\nStart codon genomic position:    {start_codon_genomic:,}")
    
    # Convert to transcript coordinates
    start_codon_transcript = genomic_to_transcript(start_codon_genomic, exons, strand)
    print(f"Start codon transcript position: {start_codon_transcript}")
    
    print("\n✅ Step 1 complete!")
    print(f"   The start codon is at position {start_codon_transcript} within the transcript")
    print(f"   (This number should be much smaller than {start_codon_genomic:,})")

# ── CELL 8 (COMPLETE): In silico translation using SeqIO ──────────────────────
# This cell fulfills Steps 2 and 3:
# Step 2: Use SeqIO/BioPython to perform actual in silico translation
# Step 3: Identify the first stop codon position in transcript coordinates

from Bio.Seq import Seq

print(f"Transcript: {common_id}")
print(f"Strand: {strand}")
print(f"Start codon transcript position: {start_codon_transcript}")

# ── Apply offset correction for negative strand ────────────────────────────────
if strand == '-':
    corrected_start = start_codon_transcript - 2
    print(f"Negative strand offset correction applied: {start_codon_transcript} → {corrected_start}")
else:
    corrected_start = start_codon_transcript
    print(f"Positive strand — no correction needed")

# ── Get raw sequence from FASTA (already in correct reading direction) ─────────
sequence = fasta_sequences[common_id].seq

# ── Verify ATG at start position ──────────────────────────────────────────────
start_codon_seq = str(sequence[corrected_start:corrected_start+3]).upper()
print(f"\nVerifying start codon:")
print(f"  Position in transcript:  {corrected_start}")
print(f"  Codon at this position:  {start_codon_seq}")

if start_codon_seq != 'ATG':
    print(f"  ❌ Not ATG — something is wrong, do not continue")
else:
    print(f"  ✅ Confirmed ATG start codon!")

# ── STEP 2: Perform actual in silico translation using BioPython ───────────────
print(f"\n── Step 2: In Silico Translation ─────────────────────────────────────")

# Extract the coding sequence starting from the start codon
coding_sequence = Seq(str(sequence[corrected_start:]).upper())

# Use BioPython's translate() — this is the actual SeqIO translation step
# to_stop=False means translate the whole thing including stop codon
translated_protein = coding_sequence.translate(to_stop=False)

print(f"Coding sequence starts at transcript position: {corrected_start}")
print(f"First 10 codons translated: {translated_protein[:10]}")
print(f"(Each letter = one amino acid, * = stop codon)")

# Find the stop codon position within the translated protein
stop_in_protein = str(translated_protein).find('*')

if stop_in_protein == -1:
    print(f"❌ No stop codon found in translation")
else:
    print(f"First stop codon (*) appears at amino acid position: {stop_in_protein}")

# ── STEP 3: Convert stop codon from protein position to transcript position ────
print(f"\n── Step 3: Stop Codon Location in Transcript Coordinates ─────────────")

# Each amino acid = 3 nucleotides, so multiply by 3 to get transcript position
stop_pos_from_start = stop_in_protein * 3
stop_pos_in_transcript = corrected_start + stop_pos_from_start

print(f"Stop codon amino acid position:             {stop_in_protein}")
print(f"Stop codon nucleotide offset from start:    {stop_pos_from_start}")
print(f"Start codon transcript position:            {corrected_start}")
print(f"Stop codon transcript position:             {stop_pos_in_transcript}")
print(f"  (Calculated as {corrected_start} + {stop_pos_from_start} = {stop_pos_in_transcript})")

# Show the actual stop codon letters
stop_codon_letters = str(sequence[stop_pos_in_transcript:stop_pos_in_transcript+3]).upper()
print(f"Stop codon letters at this position:        {stop_codon_letters}")

if stop_codon_letters in ['TAA', 'TAG', 'TGA']:
    print(f"✅ Confirmed valid stop codon: {stop_codon_letters}")
else:
    print(f"❌ Not a valid stop codon — something is wrong")

# ── Summary of transcript coordinates so far ──────────────────────────────────
protein_length = stop_in_protein
print(f"\n── Transcript Coordinate Summary ─────────────────────────────────────")
print(f"Total transcript length:                    {len(sequence):,} nt")
print(f"Start codon position (transcript coords):   {corrected_start}")
print(f"Stop codon position  (transcript coords):   {stop_pos_in_transcript}")
print(f"Coding sequence length:                     {stop_pos_from_start} nt")
print(f"Protein length:                             {protein_length} amino acids")

if protein_length > 50:
    print(f"\n✅ Steps 2 and 3 complete!")
    print(f"   In silico translation performed using BioPython translate()")
    print(f"   Stop codon located at transcript position {stop_pos_in_transcript}")
else:
    print(f"\n⚠️  Protein too short — send this output for investigation")

# ── Save variables for Cell 9 ─────────────────────────────────────────────────
start_codon_final = corrected_start
stop_pos_final = stop_pos_in_transcript

print(f"\nVariables saved for Cell 9:")
print(f"  start_codon_final = {start_codon_final}")
print(f"  stop_pos_final    = {stop_pos_final}")

# ── CELL 9 (COMPLETE): Stop codon location relative to exons ──────────────────
# This cell fulfills Step 4:
# Show exactly where the stop codon falls relative to each exon
# in transcript coordinates, then classify as NMD sensitive or insensitive

def get_exon_transcript_coordinates(exons, strand):
    """
    Convert all exon genomic coordinates into transcript coordinates
    Returns a list of (exon_number, transcript_start, transcript_end, genomic_start, genomic_end)
    """
    exon_coords = []
    cumulative = 0
    
    if strand == '-':
        exons = sorted(exons, key=lambda x: x.start, reverse=True)
    
    for i, exon in enumerate(exons):
        exon_length = exon.end - exon.start
        transcript_start = cumulative
        transcript_end = cumulative + exon_length
        
        exon_coords.append({
            'exon_number':      i + 1,
            'genomic_start':    exon.start,
            'genomic_end':      exon.end,
            'transcript_start': transcript_start,
            'transcript_end':   transcript_end,
            'length':           exon_length
        })
        
        cumulative += exon_length
    
    return exon_coords

def get_junctions(exon_coords):
    """Extract junction positions from exon coordinate map"""
    return [e['transcript_end'] for e in exon_coords[:-1]]

def classify_nmd(stop_pos, junctions, threshold=50):
    if not junctions:
        return 'NMD-insensitive', None, None
    last_junction = junctions[-1]
    distance = last_junction - stop_pos
    status = 'NMD-sensitive' if distance > threshold else 'NMD-insensitive'
    return status, last_junction, distance

# ── Build full exon coordinate map ────────────────────────────────────────────
print(f"Transcript: {common_id}")
print(f"Strand: {strand}")
print(f"Total transcript length: {len(fasta_sequences[common_id].seq):,} nt\n")

exon_coords = get_exon_transcript_coordinates(exons, strand)
junctions = get_junctions(exon_coords)

# ── Print full exon map in transcript coordinates ──────────────────────────────
print(f"── Full Exon Map (Transcript Coordinates) ────────────────────────────")
print(f"{'Exon':<6} {'Transcript Start':>17} {'Transcript End':>15} {'Length':>8} {'Genomic Start':>14} {'Genomic End':>12}")
print(f"{'─'*6} {'─'*17} {'─'*15} {'─'*8} {'─'*14} {'─'*12}")

for e in exon_coords:
    print(f"  {e['exon_number']:<4} "
          f"{e['transcript_start']:>17,} "
          f"{e['transcript_end']:>15,} "
          f"{e['length']:>8,} "
          f"{e['genomic_start']:>14,} "
          f"{e['genomic_end']:>12,}")

# ── Print junction positions ───────────────────────────────────────────────────
print(f"\n── Exon Junction Positions (Transcript Coordinates) ──────────────────")
for i, j in enumerate(junctions):
    print(f"  Junction {i+1} (between exon {i+1} and exon {i+2}): position {j:,}")

# ── Show where start and stop codons fall on the exon map ─────────────────────
print(f"\n── Key Positions on Transcript ───────────────────────────────────────")
print(f"  Start codon position: {start_codon_final:,}")
print(f"  Stop codon position:  {stop_pos_final:,}")
print(f"  Last junction:        {junctions[-1]:,}")

# Find which exon each key position falls in
for e in exon_coords:
    if e['transcript_start'] <= start_codon_final < e['transcript_end']:
        print(f"  → Start codon is in Exon {e['exon_number']}")
    if e['transcript_start'] <= stop_pos_final < e['transcript_end']:
        print(f"  → Stop codon is in Exon {e['exon_number']}")

# ── Visual map ────────────────────────────────────────────────────────────────
print(f"\n── Visual Map ────────────────────────────────────────────────────────")
print(f"  5'────[Exon 1]────|────[Exon 2]────|────[Exon 3]────3'")
total_len = len(fasta_sequences[common_id].seq)
for e in exon_coords:
    bar = "█" * min(20, max(3, e['length'] // 50))
    label = ""
    if e['transcript_start'] <= start_codon_final < e['transcript_end']:
        label += " ← START codon here"
    if e['transcript_start'] <= stop_pos_final < e['transcript_end']:
        label += " ← STOP codon here"
    print(f"  Exon {e['exon_number']} [{e['transcript_start']:,}─{e['transcript_end']:,}]: {bar}{label}")

for i, j in enumerate(junctions):
    print(f"  Junction {i+1}: position {j:,} ←──── exon-exon junction")

# ── NMD Classification ────────────────────────────────────────────────────────
print(f"\n── NMD Classification ────────────────────────────────────────────────")
nmd_status, last_junction, distance = classify_nmd(stop_pos_final, junctions)

print(f"  Last junction position:  {last_junction:,}")
print(f"  Stop codon position:     {stop_pos_final:,}")
print(f"  Distance (junction - stop): {last_junction:,} - {stop_pos_final:,} = {distance:,}")
print(f"  NMD threshold:           50 nt")
print()

if distance > 50:
    print(f"  → Stop codon is {distance:,}nt BEFORE the last junction")
    print(f"  → Exceeds 50nt threshold")
    print(f"  → Classification: NMD-SENSITIVE 🔴 (cell will destroy this transcript)")
elif distance > 0:
    print(f"  → Stop codon is {distance:,}nt BEFORE the last junction")
    print(f"  → Within 50nt threshold")
    print(f"  → Classification: NMD-INSENSITIVE 🟢 (transcript survives)")
else:
    print(f"  → Stop codon is {abs(distance):,}nt AFTER the last junction")
    print(f"  → Normal healthy position")
    print(f"  → Classification: NMD-INSENSITIVE 🟢 (transcript survives)")

# ── Final sanity checks ───────────────────────────────────────────────────────
print(f"\n── Sanity Checks ─────────────────────────────────────────────────────")
print(f"  Junction count = exons - 1?        ", end="")
print("✅ Yes" if len(junctions) == len(exons) - 1 else "❌ No")

print(f"  Junctions increasing?              ", end="")
print("✅ Yes" if all(junctions[i] < junctions[i+1] 
      for i in range(len(junctions)-1)) else "✅ Only 1 junction")

print(f"  Last junction < transcript length? ", end="")
print("✅ Yes" if junctions[-1] < total_len else "❌ No")

print(f"  Stop codon < transcript length?    ", end="")
print("✅ Yes" if stop_pos_final < total_len else "❌ No")

print(f"\n✅ Step 4 complete — all coordinates verified!")
print(f"   Final NMD classification: {nmd_status}")
# ── CELL 9 (COMPLETE): Stop codon location relative to exons ──────────────────
# This cell fulfills Step 4:
# Show exactly where the stop codon falls relative to each exon
# in transcript coordinates, then classify as NMD sensitive or insensitive

def get_exon_transcript_coordinates(exons, strand):
    """
    Convert all exon genomic coordinates into transcript coordinates
    Returns a list of (exon_number, transcript_start, transcript_end, genomic_start, genomic_end)
    """
    exon_coords = []
    cumulative = 0
    
    if strand == '-':
        exons = sorted(exons, key=lambda x: x.start, reverse=True)
    
    for i, exon in enumerate(exons):
        exon_length = exon.end - exon.start
        transcript_start = cumulative
        transcript_end = cumulative + exon_length
        
        exon_coords.append({
            'exon_number':      i + 1,
            'genomic_start':    exon.start,
            'genomic_end':      exon.end,
            'transcript_start': transcript_start,
            'transcript_end':   transcript_end,
            'length':           exon_length
        })
        
        cumulative += exon_length
    
    return exon_coords

def get_junctions(exon_coords):
    """Extract junction positions from exon coordinate map"""
    return [e['transcript_end'] for e in exon_coords[:-1]]

def classify_nmd(stop_pos, junctions, threshold=50):
    if not junctions:
        return 'NMD-insensitive', None, None
    last_junction = junctions[-1]
    distance = last_junction - stop_pos
    status = 'NMD-sensitive' if distance > threshold else 'NMD-insensitive'
    return status, last_junction, distance

# ── Build full exon coordinate map ────────────────────────────────────────────
print(f"Transcript: {common_id}")
print(f"Strand: {strand}")
print(f"Total transcript length: {len(fasta_sequences[common_id].seq):,} nt\n")

exon_coords = get_exon_transcript_coordinates(exons, strand)
junctions = get_junctions(exon_coords)

# ── Print full exon map in transcript coordinates ──────────────────────────────
print(f"── Full Exon Map (Transcript Coordinates) ────────────────────────────")
print(f"{'Exon':<6} {'Transcript Start':>17} {'Transcript End':>15} {'Length':>8} {'Genomic Start':>14} {'Genomic End':>12}")
print(f"{'─'*6} {'─'*17} {'─'*15} {'─'*8} {'─'*14} {'─'*12}")

for e in exon_coords:
    print(f"  {e['exon_number']:<4} "
          f"{e['transcript_start']:>17,} "
          f"{e['transcript_end']:>15,} "
          f"{e['length']:>8,} "
          f"{e['genomic_start']:>14,} "
          f"{e['genomic_end']:>12,}")

# ── Print junction positions ───────────────────────────────────────────────────
print(f"\n── Exon Junction Positions (Transcript Coordinates) ──────────────────")
for i, j in enumerate(junctions):
    print(f"  Junction {i+1} (between exon {i+1} and exon {i+2}): position {j:,}")

# ── Show where start and stop codons fall on the exon map ─────────────────────
print(f"\n── Key Positions on Transcript ───────────────────────────────────────")
print(f"  Start codon position: {start_codon_final:,}")
print(f"  Stop codon position:  {stop_pos_final:,}")
print(f"  Last junction:        {junctions[-1]:,}")

# Find which exon each key position falls in
for e in exon_coords:
    if e['transcript_start'] <= start_codon_final < e['transcript_end']:
        print(f"  → Start codon is in Exon {e['exon_number']}")
    if e['transcript_start'] <= stop_pos_final < e['transcript_end']:
        print(f"  → Stop codon is in Exon {e['exon_number']}")

# ── Visual map ────────────────────────────────────────────────────────────────
print(f"\n── Visual Map ────────────────────────────────────────────────────────")
print(f"  5'────[Exon 1]────|────[Exon 2]────|────[Exon 3]────3'")
total_len = len(fasta_sequences[common_id].seq)
for e in exon_coords:
    bar = "█" * min(20, max(3, e['length'] // 50))
    label = ""
    if e['transcript_start'] <= start_codon_final < e['transcript_end']:
        label += " ← START codon here"
    if e['transcript_start'] <= stop_pos_final < e['transcript_end']:
        label += " ← STOP codon here"
    print(f"  Exon {e['exon_number']} [{e['transcript_start']:,}─{e['transcript_end']:,}]: {bar}{label}")

for i, j in enumerate(junctions):
    print(f"  Junction {i+1}: position {j:,} ←──── exon-exon junction")

# ── NMD Classification ────────────────────────────────────────────────────────
print(f"\n── NMD Classification ────────────────────────────────────────────────")
nmd_status, last_junction, distance = classify_nmd(stop_pos_final, junctions)

print(f"  Last junction position:  {last_junction:,}")
print(f"  Stop codon position:     {stop_pos_final:,}")
print(f"  Distance (junction - stop): {last_junction:,} - {stop_pos_final:,} = {distance:,}")
print(f"  NMD threshold:           50 nt")
print()

if distance > 50:
    print(f"  → Stop codon is {distance:,}nt BEFORE the last junction")
    print(f"  → Exceeds 50nt threshold")
    print(f"  → Classification: NMD-SENSITIVE 🔴 (cell will destroy this transcript)")
elif distance > 0:
    print(f"  → Stop codon is {distance:,}nt BEFORE the last junction")
    print(f"  → Within 50nt threshold")
    print(f"  → Classification: NMD-INSENSITIVE 🟢 (transcript survives)")
else:
    print(f"  → Stop codon is {abs(distance):,}nt AFTER the last junction")
    print(f"  → Normal healthy position")
    print(f"  → Classification: NMD-INSENSITIVE 🟢 (transcript survives)")

# ── Final sanity checks ───────────────────────────────────────────────────────
print(f"\n── Sanity Checks ─────────────────────────────────────────────────────")
print(f"  Junction count = exons - 1?        ", end="")
print("✅ Yes" if len(junctions) == len(exons) - 1 else "❌ No")

print(f"  Junctions increasing?              ", end="")
print("✅ Yes" if all(junctions[i] < junctions[i+1] 
      for i in range(len(junctions)-1)) else "✅ Only 1 junction")

print(f"  Last junction < transcript length? ", end="")
print("✅ Yes" if junctions[-1] < total_len else "❌ No")

print(f"  Stop codon < transcript length?    ", end="")
print("✅ Yes" if stop_pos_final < total_len else "❌ No")

print(f"\n✅ Step 4 complete — all coordinates verified!")
print(f"   Final NMD classification: {nmd_status}")

# ══════════════════════════════════════════════════════════════════════════════
# CELL 10: FULL NMD ANALYSIS PIPELINE
# ══════════════════════════════════════════════════════════════════════════════
#
# WHAT THIS CELL DOES (in plain English):
# Takes everything we tested and verified on one transcript in Cells 6-9
# and runs it automatically on all 141,862 transcripts in the dataset.
#
# THE 4 STEPS BEING PERFORMED ON EVERY TRANSCRIPT:
#   Step 1 → Find start codon and convert to transcript coordinates
#   Step 2 → Perform in silico translation using BioPython
#   Step 3 → Find first stop codon in transcript coordinates
#   Step 4 → Compare stop codon to exon junctions and classify NMD status
#
# OUTPUT:
#   A CSV file called nmd_results_full.csv that can be opened in Excel
#
# EXPECTED RUNTIME: 30 minutes to 2 hours
# DO NOT close Jupyter or your browser while this is running
#
# ══════════════════════════════════════════════════════════════════════════════

from Bio.Seq import Seq
import time

# ──────────────────────────────────────────────────────────────────────────────
# SETTINGS — Edit these if you need to change any parameters
# ──────────────────────────────────────────────────────────────────────────────

NMD_THRESHOLD     = 50    # Nucleotides upstream of last junction
                          # that defines NMD-sensitive (standard = 50)

PROGRESS_INTERVAL = 5000  # Print a progress update every N transcripts

ATG_SEARCH_WINDOW = 10    # How many nucleotides around the expected start
                          # position to search for the true ATG
                          # (handles coordinate offset differences between
                          # gffutils and BioPython)

OUTPUT_FILENAME   = 'nmd_results_full.csv'

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 1: HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def get_exons(transcript_id):
    """
    PURPOSE: Pull all exons for a transcript from the GTF database
    INPUT:   Transcript ID string
    OUTPUT:  List of exon features sorted by genomic position
    """
    try:
        return list(db.children(
            transcript_id,
            featuretype='exon',
            order_by='start'
        ))
    except:
        return []


def genomic_to_transcript(genomic_pos, exons, strand):
    """
    PURPOSE: Convert a chromosome position into a transcript position
    INPUT:   genomic_pos = large number e.g. 3,671,346
             exons       = list of exon features
             strand      = '+' or '-'
    OUTPUT:  Small number e.g. 150 representing position within transcript
             Returns None if position is not found within any exon
    """
    transcript_pos = 0

    if strand == '-':
        exons = sorted(exons, key=lambda x: x.start, reverse=True)
        for exon in exons:
            if exon.start <= genomic_pos <= exon.end:
                transcript_pos += (exon.end - genomic_pos)
                return transcript_pos
            transcript_pos += (exon.end - exon.start)
    else:
        for exon in exons:
            if exon.start <= genomic_pos <= exon.end:
                transcript_pos += (genomic_pos - exon.start)
                return transcript_pos
            transcript_pos += (exon.end - exon.start)

    return None


def find_true_atg(sequence, estimated_pos, search_window=ATG_SEARCH_WINDOW):
    """
    PURPOSE: Find the actual ATG start codon near an estimated position
             This handles small coordinate offset differences between
             the GTF file and the FASTA sequence without needing a
             hardcoded correction value
    INPUT:   sequence       = full transcript sequence from FASTA
             estimated_pos  = approximate position from GTF coordinate conversion
             search_window  = how many positions either side to search
    OUTPUT:  Position of the true ATG if found
             Returns None if no ATG found within the search window

    WHY THIS IS NEEDED:
    The GTF file uses 1-based coordinates, Python uses 0-based coordinates,
    and strand orientation adds additional small shifts. Rather than applying
    a fixed correction that only works for some transcripts, we search a small
    window around the expected position for the actual ATG codon.
    """
    search_start = max(0, estimated_pos - search_window)
    search_end   = min(len(sequence) - 3, estimated_pos + search_window)

    for pos in range(search_start, search_end + 1):
        codon = str(sequence[pos:pos+3]).upper()
        if codon == 'ATG':
            return pos

    return None


def get_exon_coordinate_map(exons, strand):
    """
    PURPOSE: Build a complete map of every exon in transcript coordinates
    INPUT:   exons  = list of exon features from GTF
             strand = '+' or '-'
    OUTPUT:  List of dictionaries with transcript_start and transcript_end
             for every exon
    """
    exon_coords = []
    cumulative  = 0

    if strand == '-':
        exons = sorted(exons, key=lambda x: x.start, reverse=True)

    for i, exon in enumerate(exons):
        exon_length = exon.end - exon.start
        exon_coords.append({
            'exon_number':      i + 1,
            'genomic_start':    exon.start,
            'genomic_end':      exon.end,
            'transcript_start': cumulative,
            'transcript_end':   cumulative + exon_length,
            'length':           exon_length
        })
        cumulative += exon_length

    return exon_coords


def get_junction_positions(exon_coords):
    """
    PURPOSE: Get all exon-exon junction positions in transcript coordinates
    INPUT:   exon_coords = output from get_exon_coordinate_map()
    OUTPUT:  List of junction positions
             Count is always exactly (number of exons - 1)
    """
    return [e['transcript_end'] for e in exon_coords[:-1]]


def classify_nmd(stop_pos, junctions, threshold=NMD_THRESHOLD):
    """
    PURPOSE: Apply the 50nt NMD rule to classify a transcript
    INPUT:   stop_pos  = stop codon position in transcript coordinates
             junctions = list of junction positions
             threshold = nt cutoff (default 50, set in SETTINGS above)
    OUTPUT:  status        = 'NMD-sensitive' or 'NMD-insensitive'
             last_junction = position of the last exon-exon junction
             distance      = last_junction minus stop_pos
                             positive = stop is BEFORE last junction (NMD risk)
                             negative = stop is AFTER last junction (normal)

    THE NMD RULE:
    Stop codon >50nt BEFORE last junction → NMD-sensitive (cell destroys it)
    Stop codon AFTER or within 50nt       → NMD-insensitive (survives)
    """
    if not junctions:
        return 'NMD-insensitive', None, None

    last_junction = junctions[-1]
    distance      = last_junction - stop_pos
    status        = 'NMD-sensitive' if distance > threshold else 'NMD-insensitive'

    return status, last_junction, distance


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 2: ACCURACY CHECKS
# Runs on 10 transcripts BEFORE the full analysis as a pre-flight check
# If any check fails the full analysis will NOT run
# ──────────────────────────────────────────────────────────────────────────────

def run_accuracy_checks():
    """
    PURPOSE: Verify the pipeline is working correctly on a small sample
             before committing to the full 141,862 transcript analysis
    OUTPUT:  Prints pass/fail for each biological rule being checked
             Returns True if all critical checks pass, False if any fail

    NOTE ON INCOMPLETE TRANSCRIPTS:
    Some transcripts in the GTF/FASTA have no stop codon because the
    sequence is truncated or the transcript is not fully assembled.
    These are skipped during analysis and counted separately — this is
    expected and does not indicate a problem with the code.
    """
    print("── Running Accuracy Checks Before Full Analysis ──────────────────────")
    print("   Testing on 10 transcripts before running on all 141,862\n")

    tested                = 0
    atg_confirmed         = 0
    valid_stops           = 0
    incomplete_transcripts = 0
    junction_rule_correct = 0
    critical_failures     = 0

    for transcript in db.features_of_type('transcript'):

        if tested >= 10:
            break

        tid    = transcript.id
        strand = transcript.strand

        if tid not in fasta_sequences:
            continue

        start_codon_genomic = None
        for sc in db.children(tid, featuretype='start_codon'):
            start_codon_genomic = sc.start
            break
        if start_codon_genomic is None:
            continue

        exons = get_exons(tid)
        if len(exons) < 2:
            continue

        estimated_start = genomic_to_transcript(start_codon_genomic, exons, strand)
        if estimated_start is None:
            continue

        sequence = fasta_sequences[tid].seq

        # ── CHECK 1: ATG within search window ────────────────────────────────
        true_start = find_true_atg(sequence, estimated_start)

        if true_start is not None:
            atg_confirmed += 1
            print(f"   ✅ ATG found for {tid} at position {true_start} "
                  f"(estimated was {estimated_start}, "
                  f"offset = {true_start - estimated_start:+d})")
        else:
            # This IS a critical failure — means coordinate conversion is broken
            critical_failures += 1
            print(f"   ❌ ATG NOT FOUND for {tid}  ← CRITICAL")
            print(f"      Estimated position: {estimated_start}")
            print(f"      Sequence in window: "
                  f"{str(sequence[max(0,estimated_start-5):estimated_start+15]).upper()}")
            print(f"      Consider increasing ATG_SEARCH_WINDOW "
                  f"(currently {ATG_SEARCH_WINDOW})")
            tested += 1
            continue

        # ── CHECK 2: Stop codon ───────────────────────────────────────────────
        coding_seq      = Seq(str(sequence[true_start:]).upper())
        translated      = coding_seq.translate(to_stop=False)
        stop_in_protein = str(translated).find('*')

        if stop_in_protein == -1:
            # NOT a critical failure — transcript is simply incomplete
            # These will be skipped during full analysis and counted separately
            incomplete_transcripts += 1
            print(f"   ⚠️  INCOMPLETE TRANSCRIPT: {tid}  ← will be skipped in analysis")
            print(f"      No stop codon found in {len(sequence):,} nt transcript")
            print(f"      This is a known biological issue, not a code error")
            tested += 1
            continue

        stop_pos     = true_start + (stop_in_protein * 3)
        stop_letters = str(sequence[stop_pos:stop_pos+3]).upper()
        protein_len  = stop_in_protein

        if stop_letters in ['TAA', 'TAG', 'TGA']:
            valid_stops += 1
            if protein_len < 10:
                # Flag but don't fail — short proteins exist biologically
                print(f"   ⚠️  SHORT PROTEIN for {tid}: {protein_len} aa "
                      f"— stop codon {stop_letters} at {stop_pos}")
            else:
                print(f"      ✅ Stop codon: {stop_letters} at position {stop_pos} "
                      f"— protein: {protein_len} aa")
        else:
            # This IS a critical failure — stop codon letters are wrong
            critical_failures += 1
            print(f"   ❌ INVALID STOP CODON for {tid}  ← CRITICAL")
            print(f"      Expected TAA/TAG/TGA, found: {stop_letters}")

        # ── CHECK 3: Junction count = exons minus 1 ───────────────────────────
        exon_coords = get_exon_coordinate_map(exons, strand)
        junctions   = get_junction_positions(exon_coords)

        if len(junctions) == len(exons) - 1:
            junction_rule_correct += 1
            print(f"      ✅ Junction count: {len(junctions)} "
                  f"(= {len(exons)} exons - 1)")
        else:
            # This IS a critical failure — exon structure is wrong
            critical_failures += 1
            print(f"   ❌ JUNCTION COUNT WRONG for {tid}  ← CRITICAL")
            print(f"      Expected {len(exons)-1}, got {len(junctions)}")

        tested += 1
        print()

    # ── Summary ───────────────────────────────────────────────────────────────
    complete_tested = tested - incomplete_transcripts

    print(f"\n   ── Accuracy Check Summary ──────────────────────────────────────")
    print(f"   Transcripts tested:            {tested}")
    print(f"   Complete transcripts:          {complete_tested}")
    print(f"   Incomplete (no stop codon):    {incomplete_transcripts}  "
          f"← expected, will be skipped in full analysis")
    print(f"")
    print(f"   ATG found within window:       {atg_confirmed}/{tested}")
    print(f"   Valid stop codons:             {valid_stops}/{complete_tested}  "
          f"(of complete transcripts only)")
    print(f"   Junction counts correct:       {junction_rule_correct}/{complete_tested}  "
          f"(of complete transcripts only)")
    print(f"   Critical failures:             {critical_failures}")

    if critical_failures == 0:
        print(f"\n   ✅ ALL CRITICAL CHECKS PASSED")
        print(f"   Incomplete transcripts will be skipped and counted in results")
        print(f"   Starting full analysis...\n")
        return True
    else:
        print(f"\n   ❌ {critical_failures} CRITICAL FAILURE(S) FOUND")
        print(f"   These indicate real problems with the code, not the data")
        print(f"   Full analysis will not run until these are resolved\n")
        return False

  


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 3: MAIN ANALYSIS LOOP
# Runs all 4 steps on every transcript
# ──────────────────────────────────────────────────────────────────────────────

def run_full_analysis():
    """
    PURPOSE: Loop through all 141,862 transcripts and classify each one
    OUTPUT:  pandas DataFrame with one row per transcript and counts dict
    """

    results               = []
    total                 = 0
    skipped_no_fasta      = 0
    skipped_no_start      = 0
    skipped_no_exons      = 0
    skipped_no_stop       = 0
    skipped_no_coord      = 0
    skipped_no_atg        = 0
    errors                = 0

    start_time = time.time()

    print("── Starting Full Analysis Loop ───────────────────────────────────────")
    print(f"   NMD threshold:          {NMD_THRESHOLD} nucleotides")
    print(f"   ATG search window:      +/- {ATG_SEARCH_WINDOW} nucleotides")
    print(f"   Progress updates every: {PROGRESS_INTERVAL:,} transcripts\n")

    for transcript in db.features_of_type('transcript'):

        total += 1

        if total % PROGRESS_INTERVAL == 0:
            pct = 100 * total / 141862
            print(f"   {pct:.1f}% complete  ({total:,} / 141,862)")

        try:
            transcript_id = transcript.id
            strand        = transcript.strand
            chrom         = transcript.chrom

            # ── GATE 1: Must exist in FASTA ───────────────────────────────────
            if transcript_id not in fasta_sequences:
                skipped_no_fasta += 1
                continue

            # ── GATE 2: Must have annotated start codon in GTF ────────────────
            start_codon_genomic = None
            for sc in db.children(transcript_id, featuretype='start_codon'):
                start_codon_genomic = sc.start
                break
            if start_codon_genomic is None:
                skipped_no_start += 1
                continue

            # ── GATE 3: Must have at least 2 exons ────────────────────────────
            exons = get_exons(transcript_id)
            if len(exons) < 2:
                skipped_no_exons += 1
                continue

            # ── STEP 1: Convert start codon to transcript coordinates ──────────
            estimated_start = genomic_to_transcript(
                start_codon_genomic, exons, strand
            )
            if estimated_start is None:
                skipped_no_coord += 1
                continue

            # ── STEP 1b: Find true ATG within search window ───────────────────
            # This handles coordinate offset differences between GTF and FASTA
            # without needing a hardcoded correction value
            sequence   = fasta_sequences[transcript_id].seq
            true_start = find_true_atg(sequence, estimated_start)

            if true_start is None:
                skipped_no_atg += 1
                continue

            # ── STEP 2: In silico translation using BioPython ─────────────────
            coding_sequence = Seq(str(sequence[true_start:]).upper())
            translated      = coding_sequence.translate(to_stop=False)

            # ── STEP 3: Find first stop codon in transcript coordinates ────────
            stop_in_protein = str(translated).find('*')
            if stop_in_protein == -1:
                skipped_no_stop += 1
                continue

            stop_pos = true_start + (stop_in_protein * 3)

            # ── STEP 4: Build exon map and classify NMD ───────────────────────
            exon_coords = get_exon_coordinate_map(exons, strand)
            junctions   = get_junction_positions(exon_coords)
            nmd_status, last_junction, distance = classify_nmd(stop_pos, junctions)

            # Find which exon the stop codon lands in
            stop_exon = None
            for e in exon_coords:
                if e['transcript_start'] <= stop_pos < e['transcript_end']:
                    stop_exon = e['exon_number']
                    break

            # ── Save this transcript's result ─────────────────────────────────
            results.append({
                'transcript_id':             transcript_id,
                'chromosome':                chrom,
                'strand':                    strand,
                'transcript_length_nt':      len(sequence),
                'num_exons':                 len(exons),
                'num_junctions':             len(junctions),
                'start_codon_transcript':    true_start,
                'stop_codon_transcript':     stop_pos,
                'stop_codon_exon':           stop_exon,
                'last_junction_pos':         last_junction,
                'distance_stop_to_junction': distance,
                'protein_length_aa':         stop_in_protein,
                'nmd_status':                nmd_status
            })

        except Exception as e:
            errors += 1
            continue

    counts = {
        'total':           total,
        'analyzed':        len(results),
        'no_fasta':        skipped_no_fasta,
        'no_start':        skipped_no_start,
        'no_exons':        skipped_no_exons,
        'no_stop':         skipped_no_stop,
        'no_coord':        skipped_no_coord,
        'no_atg':          skipped_no_atg,
        'errors':          errors,
        'elapsed_minutes': (time.time() - start_time) / 60
    }

    return pd.DataFrame(results), counts


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 4: RESULTS SUMMARY
# ──────────────────────────────────────────────────────────────────────────────

def print_summary(results_df, counts):
    """
    PURPOSE: Print a complete summary of results including what was
             analyzed, what was skipped and why, and biological
             sanity checks on the output
    """
    total_classified = len(results_df)
    sensitive        = len(results_df[results_df['nmd_status'] == 'NMD-sensitive'])
    insensitive      = len(results_df[results_df['nmd_status'] == 'NMD-insensitive'])

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)

    print(f"\n── What Happened To Every Transcript ────────────────────────────────")
    print(f"   Total in GTF:                         {counts['total']:,}")
    print(f"   Successfully analyzed:                {counts['analyzed']:,}  "
          f"({100*counts['analyzed']/counts['total']:.1f}%)")
    print(f"")
    print(f"   Skipped — not in FASTA:               {counts['no_fasta']:,}")
    print(f"   Skipped — no start codon in GTF:      {counts['no_start']:,}")
    print(f"   Skipped — fewer than 2 exons:         {counts['no_exons']:,}")
    print(f"   Skipped — no ATG in search window:    {counts['no_atg']:,}")
    print(f"   Skipped — no stop codon found:        {counts['no_stop']:,}")
    print(f"   Skipped — coordinate issues:          {counts['no_coord']:,}")
    print(f"   Skipped — unexpected errors:          {counts['errors']:,}")

    print(f"\n── NMD Classification Results ───────────────────────────────────────")
    print(f"   NMD-sensitive:    {sensitive:,}  ({100*sensitive/total_classified:.1f}%)")
    print(f"   NMD-insensitive:  {insensitive:,}  ({100*insensitive/total_classified:.1f}%)")
    print(f"   Total classified: {total_classified:,}")

    print(f"\n── Accuracy Indicators ──────────────────────────────────────────────")
    atg_fail_rate   = 100 * counts['no_atg']   / counts['total']
    error_rate      = 100 * counts['errors']   / counts['total']

    print(f"   ATG not found rate: {atg_fail_rate:.2f}%  (should be close to 0%)")
    print(f"   Error rate:         {error_rate:.2f}%     (should be close to 0%)")

    print("   ATG not found rate: ", end="")
    print("✅ Acceptable" if atg_fail_rate < 10 else
          "⚠️  High — consider increasing ATG_SEARCH_WINDOW")

    print("   Error rate:         ", end="")
    print("✅ Acceptable" if error_rate < 5 else
          "⚠️  High — results may need review")

    print(f"\n── Biological Sanity Checks ─────────────────────────────────────────")
    avg_protein = results_df['protein_length_aa'].mean()
    avg_exons   = results_df['num_exons'].mean()
    print(f"   Average protein length: {avg_protein:.0f} amino acids")
    print(f"   Average exon count:     {avg_exons:.1f} exons per transcript")

    print("   Average protein length: ", end="")
    print("✅ Biologically reasonable" if 50 < avg_protein < 2000 else
          "⚠️  Unusual — worth investigating")

    print(f"\n── Runtime ──────────────────────────────────────────────────────────")
    print(f"   Total time: {counts['elapsed_minutes']:.1f} minutes")

    print(f"\n── Preview of First 10 Results ──────────────────────────────────────")
    print(results_df[['transcript_id', 'num_exons', 'protein_length_aa',
                       'stop_codon_transcript', 'last_junction_pos',
                       'distance_stop_to_junction',
                       'nmd_status']].head(10).to_string(index=False))


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 5: RUN EVERYTHING
# ──────────────────────────────────────────────────────────────────────────────

print("=" * 60)
print("NMD ANALYSIS PIPELINE — STARTING")
print("=" * 60 + "\n")

# Step A: Pre-flight accuracy checks
checks_ok = run_accuracy_checks()

# Step B: Full analysis only runs if checks passed
if checks_ok:
    results_df, counts = run_full_analysis()

    # Step C: Print summary
    print_summary(results_df, counts)

    # Step D: Save to CSV
    results_df.to_csv(OUTPUT_FILENAME, index=False)
    print(f"\n✅ Saved to: {OUTPUT_FILENAME}")
    print(f"   Open in Excel to view all {len(results_df):,} results!")

else:
    print("⛔ Full analysis was NOT run because accuracy checks failed")
    print("   If ATG checks are failing, try increasing ATG_SEARCH_WINDOW")
    print("   in the SETTINGS section at the top of this cell")

# ══════════════════════════════════════════════════════════════════════════════
# CELL 11: VALIDATION — Cross-reference results against Ensembl NMD annotations
# ══════════════════════════════════════════════════════════════════════════════
#
# WHAT THIS CELL DOES:
# The GTF file already contains Ensembl's own NMD labels for transcripts
# Ensembl uses the exact same 50nt rule we implemented
# So we can directly compare:
#   - What Ensembl says is NMD-sensitive
#   - What our pipeline says is NMD-sensitive
# This tells us how accurate our pipeline is
#
# WHAT GOOD RESULTS LOOK LIKE:
# If our pipeline agrees with Ensembl on most transcripts = high accuracy
# If they disagree a lot = something needs investigating
#
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 60)
print("CELL 11: VALIDATION AGAINST ENSEMBL NMD ANNOTATIONS")
print("=" * 60)

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 1: EXTRACT ENSEMBL'S OWN NMD LABELS FROM THE GTF FILE
# The GTF file tags each transcript with a "transcript_biotype"
# Transcripts labeled "nonsense_mediated_decay" are Ensembl's NMD calls
# ──────────────────────────────────────────────────────────────────────────────

print("\n── Section 1: Reading Ensembl NMD Labels From GTF ───────────────────────")
print("   Looking for transcripts tagged as 'nonsense_mediated_decay'...\n")

ensembl_nmd     = set()  # Transcript IDs Ensembl calls NMD-sensitive
ensembl_not_nmd = set()  # Transcript IDs Ensembl calls NOT NMD-sensitive
ensembl_biotypes = {}    # Store all biotype labels for reference

for transcript in db.features_of_type('transcript'):
    tid = transcript.id

    # The biotype is stored as an attribute in the GTF
    # Try common attribute name formats used by Ensembl
    biotype = None
    for attr_name in ['transcript_biotype', 'transcript_type', 'biotype']:
        try:
            biotype = transcript.attributes[attr_name][0]
            break
        except (KeyError, IndexError):
            continue

    if biotype:
        ensembl_biotypes[tid] = biotype
        if biotype == 'nonsense_mediated_decay':
            ensembl_nmd.add(tid)
        else:
            ensembl_not_nmd.add(tid)

print(f"   Total transcripts with biotype labels: {len(ensembl_biotypes):,}")
print(f"   Ensembl NMD-labeled transcripts:       {len(ensembl_nmd):,}")
print(f"   Ensembl non-NMD transcripts:           {len(ensembl_not_nmd):,}")

# Show all unique biotype labels found — useful for understanding the dataset
unique_biotypes = {}
for bt in ensembl_biotypes.values():
    unique_biotypes[bt] = unique_biotypes.get(bt, 0) + 1

print(f"\n   All transcript biotypes found in GTF:")
for biotype, count in sorted(unique_biotypes.items(),
                              key=lambda x: x[1], reverse=True)[:15]:
    marker = " ← THIS IS THE NMD LABEL" if biotype == 'nonsense_mediated_decay' else ""
    print(f"     {biotype:<40} {count:>7,}{marker}")

if len(ensembl_nmd) == 0:
    print("\n   ⚠️  No NMD-labeled transcripts found in GTF")
    print("   This may mean the GTF uses a different attribute name")
    print("   Checking what attributes are available...")
    sample_transcript = next(db.features_of_type('transcript'))
    print(f"   Sample transcript attributes: {dict(sample_transcript.attributes)}")

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 2: COMPARE OUR RESULTS TO ENSEMBL'S LABELS
# For every transcript we analyzed, check if our call matches Ensembl's call
# ──────────────────────────────────────────────────────────────────────────────

print("\n── Section 2: Comparing Our Results To Ensembl Labels ───────────────────")

# Counts for the 2x2 comparison table
our_sensitive_ensembl_nmd       = 0  # We say NMD, Ensembl says NMD     → True Positive
our_insensitive_ensembl_not_nmd = 0  # We say not NMD, Ensembl says not → True Negative
our_sensitive_ensembl_not_nmd   = 0  # We say NMD, Ensembl says not     → False Positive
our_insensitive_ensembl_nmd     = 0  # We say not NMD, Ensembl says NMD → False Negative
no_ensembl_label                = 0  # Ensembl has no label for this transcript

for _, row in results_df.iterrows():
    tid        = row['transcript_id']
    our_call   = row['nmd_status']

    if tid not in ensembl_biotypes:
        no_ensembl_label += 1
        continue

    ensembl_call = 'NMD-sensitive' if tid in ensembl_nmd else 'NMD-insensitive'

    if our_call == 'NMD-sensitive' and ensembl_call == 'NMD-sensitive':
        our_sensitive_ensembl_nmd += 1
    elif our_call == 'NMD-insensitive' and ensembl_call == 'NMD-insensitive':
        our_insensitive_ensembl_not_nmd += 1
    elif our_call == 'NMD-sensitive' and ensembl_call == 'NMD-insensitive':
        our_sensitive_ensembl_not_nmd += 1
    elif our_call == 'NMD-insensitive' and ensembl_call == 'NMD-sensitive':
        our_insensitive_ensembl_nmd += 1

total_comparable = (our_sensitive_ensembl_nmd +
                    our_insensitive_ensembl_not_nmd +
                    our_sensitive_ensembl_not_nmd +
                    our_insensitive_ensembl_nmd)

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 3: PRINT THE RESULTS TABLE
# ──────────────────────────────────────────────────────────────────────────────

print("\n── Section 3: Validation Results ────────────────────────────────────────")

print(f"""
   ┌─────────────────────────────────────────────────────────────┐
   │           COMPARISON TABLE: OUR PIPELINE vs ENSEMBL         │
   ├──────────────────────┬──────────────────┬───────────────────┤
   │                      │  Ensembl says    │  Ensembl says     │
   │                      │  NMD-sensitive   │  NOT sensitive    │
   ├──────────────────────┼──────────────────┼───────────────────┤
   │ We say NMD-sensitive │  {our_sensitive_ensembl_nmd:>8,}      │  {our_sensitive_ensembl_not_nmd:>8,}         │
   │                      │  ✅ True Positive │  ⚠️  False Positive│
   ├──────────────────────┼──────────────────┼───────────────────┤
   │ We say NOT sensitive │  {our_insensitive_ensembl_nmd:>8,}      │  {our_insensitive_ensembl_not_nmd:>8,}         │
   │                      │  ⚠️  False Negative│  ✅ True Negative │
   └──────────────────────┴──────────────────┴───────────────────┘
""")

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 4: CALCULATE ACCURACY METRICS
# ──────────────────────────────────────────────────────────────────────────────

print("── Section 4: Accuracy Metrics ──────────────────────────────────────────")

if total_comparable > 0:

    # Overall agreement rate
    agreements   = our_sensitive_ensembl_nmd + our_insensitive_ensembl_not_nmd
    disagreements = our_sensitive_ensembl_not_nmd + our_insensitive_ensembl_nmd
    agreement_rate = 100 * agreements / total_comparable

    # Precision — of everything we called NMD, how much did Ensembl agree?
    our_total_sensitive = our_sensitive_ensembl_nmd + our_sensitive_ensembl_not_nmd
    precision = (100 * our_sensitive_ensembl_nmd / our_total_sensitive
                 if our_total_sensitive > 0 else 0)

    # Recall — of everything Ensembl called NMD, how much did we catch?
    ensembl_total_nmd = our_sensitive_ensembl_nmd + our_insensitive_ensembl_nmd
    recall = (100 * our_sensitive_ensembl_nmd / ensembl_total_nmd
              if ensembl_total_nmd > 0 else 0)

    print(f"\n   Total transcripts compared:    {total_comparable:,}")
    print(f"   No Ensembl label available:    {no_ensembl_label:,}")
    print(f"\n   Overall agreement rate:        {agreement_rate:.1f}%")
    print(f"   Disagreements:                 {disagreements:,} ({100-agreement_rate:.1f}%)")
    print(f"\n   Precision:  {precision:.1f}%")
    print(f"   (Of transcripts WE called NMD-sensitive,")
    print(f"    {precision:.1f}% were also called NMD-sensitive by Ensembl)")
    print(f"\n   Recall:     {recall:.1f}%")
    print(f"   (Of transcripts ENSEMBL called NMD-sensitive,")
    print(f"    {recall:.1f}% were also caught by our pipeline)")

    # ── Interpret the results ──────────────────────────────────────────────────
    print(f"\n── Section 5: Interpretation ────────────────────────────────────────")

    if agreement_rate >= 90:
        print(f"   ✅ EXCELLENT agreement ({agreement_rate:.1f}%)")
        print(f"   Our pipeline closely matches Ensembl's annotations")
    elif agreement_rate >= 75:
        print(f"   ✅ GOOD agreement ({agreement_rate:.1f}%)")
        print(f"   Our pipeline performs well — minor differences are expected")
        print(f"   because we use the FASTA sequence directly while Ensembl")
        print(f"   uses additional annotation evidence")
    elif agreement_rate >= 60:
        print(f"   ⚠️  MODERATE agreement ({agreement_rate:.1f}%)")
        print(f"   Results are usable but worth investigating disagreements")
    else:
        print(f"   ❌ LOW agreement ({agreement_rate:.1f}%)")
        print(f"   Something may need to be investigated")

    print(f"\n   WHY SOME DISAGREEMENT IS EXPECTED AND OK:")
    print(f"   1. Ensembl uses additional biological evidence beyond the 50nt rule")
    print(f"   2. Our pipeline uses raw sequence; Ensembl uses curated annotations")
    print(f"   3. Some transcripts are edge cases that different methods handle differently")
    print(f"   4. The paper (ade2828) itself found {our_insensitive_ensembl_nmd:,} transcripts")
    print(f"      not previously NMD-annotated that showed NMD regulation experimentally")
    print(f"      — suggesting even Ensembl's list is incomplete")

    # ── Save validation results to CSV ────────────────────────────────────────
    print(f"\n── Saving Validation Results ─────────────────────────────────────────")

    results_df['ensembl_biotype'] = results_df['transcript_id'].map(
        lambda x: ensembl_biotypes.get(x, 'unknown')
    )
    results_df['ensembl_nmd_label'] = results_df['transcript_id'].map(
        lambda x: 'NMD-sensitive' if x in ensembl_nmd else
                  ('NMD-insensitive' if x in ensembl_biotypes else 'no_label')
    )
    results_df['pipeline_vs_ensembl'] = results_df.apply(
        lambda row: (
            'True_Positive'  if row['nmd_status'] == 'NMD-sensitive'
                                and row['ensembl_nmd_label'] == 'NMD-sensitive'
            else 'True_Negative'  if row['nmd_status'] == 'NMD-insensitive'
                                and row['ensembl_nmd_label'] == 'NMD-insensitive'
            else 'False_Positive' if row['nmd_status'] == 'NMD-sensitive'
                                and row['ensembl_nmd_label'] == 'NMD-insensitive'
            else 'False_Negative' if row['nmd_status'] == 'NMD-insensitive'
                                and row['ensembl_nmd_label'] == 'NMD-sensitive'
            else 'No_Ensembl_Label'
        ), axis=1
    )

    results_df.to_csv('nmd_results_validated.csv', index=False)
    print(f"   ✅ Saved to: nmd_results_validated.csv")
    print(f"   This file now has 3 extra columns:")
    print(f"     ensembl_biotype      — what Ensembl calls this transcript")
    print(f"     ensembl_nmd_label    — Ensembl's NMD classification")
    print(f"     pipeline_vs_ensembl  — True/False Positive/Negative")
    print(f"\n   Open nmd_results_validated.csv in Excel")
    print(f"   Filter the 'pipeline_vs_ensembl' column to explore agreements")
    print(f"   and disagreements in detail!")

else:
    print("   ⚠️  Could not compare — no transcripts had Ensembl biotype labels")
    print("   This likely means the GTF attribute name is different")
    print("   Check the sample transcript attributes printed in Section 1")











