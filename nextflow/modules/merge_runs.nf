// ============================================================
// Module: MERGE_RUNS
// Collects all Salmon quant directories for a given sample_id
// and writes a paths TSV consumed by tximport in R.
// Merging happens at quantification level (not FASTQ) which
// preserves per-run QC traceability.
// ============================================================

process MERGE_RUNS {
    tag "${meta.sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/quant_merged", mode: 'copy'

    input:
    tuple val(meta), path(quant_dirs)

    output:
    path "quant_paths.tsv", emit: quant_paths_tsv

    script:
    """
    python3 - <<'EOF'
import os, glob

quant_dirs = "${quant_dirs}".split()
rows = []

for d in quant_dirs:
    quant_sf = os.path.join(d, "quant.sf")
    if not os.path.exists(quant_sf):
        raise FileNotFoundError(f"quant.sf not found in {d}")
    # Extract run_id from directory name (format: {run_id}_quant)
    run_id = os.path.basename(d).replace("_quant", "")
    rows.append(f"${meta.sample_id}\\t{run_id}\\t{os.path.abspath(quant_sf)}")

# Append to shared TSV (all samples)
with open("quant_paths.tsv", "a") as fh:
    for row in rows:
        fh.write(row + "\\n")

print(f"[MERGE] ${meta.sample_id}: {len(rows)} run(s) recorded in quant_paths.tsv")
EOF

    echo "[MERGE] $(date -u +%Y-%m-%dT%H:%M:%SZ) sample=${meta.sample_id} runs=${quant_dirs} status=COMPLETE" \\
        >> ${params.outdir}/audit/audit_trail.log
    """
}
