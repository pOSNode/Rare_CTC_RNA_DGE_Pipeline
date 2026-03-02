// ============================================================
// Module: VALIDATE_SAMPLESHEET
// Validates CSV structure, checks for required columns,
// duplicate run IDs, and missing files before pipeline starts.
// ============================================================

process VALIDATE_SAMPLESHEET {
    tag "validate"
    label 'process_low'

    input:
    path samplesheet

    output:
    path samplesheet

    script:
    """
    python3 - <<'EOF'
import csv, sys, os

required_cols = {'sample_id', 'condition', 'fastq_1', 'fastq_2'}
errors = []
warnings = []

with open("${samplesheet}") as fh:
    reader = csv.DictReader(fh)
    cols = set(reader.fieldnames or [])

    missing = required_cols - cols
    if missing:
        errors.append(f"Missing required columns: {missing}")

    if 'run_id' not in cols:
        warnings.append("Column 'run_id' not found - will default to sample_id")
    if 'batch' not in cols:
        warnings.append("Column 'batch' not found - will default to 'batch1'")

    run_ids = []
    rows = []
    for i, row in enumerate(reader, start=2):
        rows.append(row)
        rid = row.get('run_id', row.get('sample_id', ''))
        if rid in run_ids:
            errors.append(f"Row {i}: Duplicate run_id '{rid}'")
        run_ids.append(rid)

        for fq_col in ['fastq_1', 'fastq_2']:
            fq = row.get(fq_col, '')
            if fq and not os.path.exists(fq):
                errors.append(f"Row {i}: File not found: {fq}")

        cond = row.get('condition', '').strip()
        if not cond:
            errors.append(f"Row {i}: Empty condition field")

    conditions = set(r.get('condition','') for r in rows)
    if len(conditions) < 2:
        errors.append(f"At least 2 conditions required for DGE. Found: {conditions}")

for w in warnings:
    print(f"[WARN]  {w}", file=sys.stderr)

if errors:
    for e in errors:
        print(f"[ERROR] {e}", file=sys.stderr)
    sys.exit(1)

print(f"[OK] Sample sheet validated: {len(rows)} runs, {len(conditions)} conditions")
EOF
    """
}
