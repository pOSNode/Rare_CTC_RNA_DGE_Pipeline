#!/usr/bin/env python3
"""
generate_test_fastq.py
Generates minimal synthetic paired-end FASTQ files for pipeline testing.
These are NOT biologically meaningful — for architecture demonstration only.
"""

import gzip, os, random, string

random.seed(42)

SAMPLES = [
    "CTC001_R1", "CTC001_R2", "CTC002_R1", "CTC003_R1",
    "CTRL001_R1", "CTRL002_R1", "CTRL003_R1"
]

READ_LEN  = 75
NUM_READS = 1000  # Minimal — real CTC samples typically 5-20M reads
OUT_DIR   = os.path.join(os.path.dirname(__file__), "fastq")
os.makedirs(OUT_DIR, exist_ok=True)

BASES = "ACGT"
QUAL  = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"

def random_seq(length):
    return "".join(random.choices(BASES, k=length))

def write_fastq(path, sample_name, read_num, n_reads):
    with gzip.open(path, "wt") as fh:
        for i in range(n_reads):
            seq = random_seq(READ_LEN)
            fh.write(f"@{sample_name}_read{i}/{read_num}\n")
            fh.write(f"{seq}\n")
            fh.write("+\n")
            fh.write(f"{QUAL[:READ_LEN]}\n")

for sample in SAMPLES:
    for read_num in [1, 2]:
        out_path = os.path.join(OUT_DIR, f"{sample}_{read_num}.fastq.gz")
        write_fastq(out_path, sample, read_num, NUM_READS)
        print(f"  Written: {out_path}")

print(f"\nGenerated {len(SAMPLES) * 2} synthetic FASTQ files in {OUT_DIR}")
print("NOTE: These files contain random sequences for architecture testing only.")
