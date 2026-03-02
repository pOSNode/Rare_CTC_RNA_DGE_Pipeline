// ============================================================
// Module: WRITE_PROVENANCE
// Writes a machine-readable provenance JSON at pipeline end.
// Records: pipeline version, run ID, parameters, tool versions,
// input file checksums, and timestamp.
// ============================================================

process WRITE_PROVENANCE {
    tag "provenance"
    label 'process_low'

    publishDir "${params.outdir}/audit", mode: 'copy'

    input:
    val  run_id
    path completed_flag

    output:
    path "provenance.json", emit: provenance
    path "run_params.txt",  emit: params_dump

    script:
    """
    python3 - <<'EOF'
import json, subprocess, datetime, os, hashlib

def tool_version(cmd):
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        return (result.stdout + result.stderr).strip().split("\\n")[0]
    except Exception as e:
        return f"unavailable ({e})"

def md5(path):
    try:
        h = hashlib.md5()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()
    except:
        return "unavailable"

provenance = {
    "pipeline": {
        "name"    : "ctc-rnaseq-pipeline",
        "version" : "${params.pipeline_version}",
        "run_id"  : "${run_id}",
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z"
    },
    "parameters": {
        "low_input_mode"    : "${params.low_input_mode}",
        "min_cpm"           : "${params.min_cpm}",
        "min_sam"           : "${params.min_sam}",
        "fdr_cutoff"        : "${params.fdr_cutoff}",
        "lfc_treat"         : "${params.lfc_treat}",
        "prior_count"       : "${params.prior_count}",
        "exclusive_min_cpm" : "${params.exclusive_min_cpm}",
        "exclusive_min_sam" : "${params.exclusive_min_sam}",
        "priority_genes"    : "${params.priority_genes}",
        "threads"           : "${params.threads}"
    },
    "tool_versions": {
        "fastp"   : tool_version(["fastp", "--version"]),
        "salmon"  : tool_version(["salmon", "--version"]),
        "multiqc" : tool_version(["multiqc", "--version"]),
        "R"       : tool_version(["Rscript", "--version"])
    },
    "input_checksums": {
        "samplesheet" : md5("${params.samplesheet}"),
        "gtf"         : md5("${params.gtf}")
    },
    "output_directory": "${params.outdir}"
}

with open("provenance.json", "w") as fh:
    json.dump(provenance, fh, indent=2)

print("[PROVENANCE] Written to provenance.json")
EOF

    # Human-readable param dump alongside JSON
    echo "=== CTC RNA-seq Pipeline Run Parameters ===" > run_params.txt
    echo "Run ID       : ${run_id}"                   >> run_params.txt
    echo "Version      : ${params.pipeline_version}"  >> run_params.txt
    echo "Timestamp    : \$(date -u)"                 >> run_params.txt
    echo "Low-input    : ${params.low_input_mode}"    >> run_params.txt
    echo "FDR cutoff   : ${params.fdr_cutoff}"        >> run_params.txt
    echo "Min CPM      : ${params.min_cpm}"           >> run_params.txt
    echo "LFC threshold: ${params.lfc_treat}"         >> run_params.txt
    echo "Priority genes: ${params.priority_genes}"   >> run_params.txt
    echo "Samplesheet  : ${params.samplesheet}"       >> run_params.txt
    echo "GTF          : ${params.gtf}"               >> run_params.txt
    echo "Output dir   : ${params.outdir}"            >> run_params.txt

    echo "[PROVENANCE] $(date -u +%Y-%m-%dT%H:%M:%SZ) run_id=${run_id} status=COMPLETE" \\
        >> ${params.outdir}/audit/audit_trail.log
    """
}
