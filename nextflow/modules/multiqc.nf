// ============================================================
// Module: MULTIQC
// Aggregates fastp and Salmon logs into a single QC report.
// ============================================================

process MULTIQC {
    tag "multiqc"
    label 'process_low'

    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data_dir

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --title "CTC RNA-seq QC Report" \\
        --comment "Low-input CTC pipeline — interpret mapping rates in context of rare-cell input" \\
        --force

    echo "[MULTIQC] $(date -u +%Y-%m-%dT%H:%M:%SZ) status=COMPLETE report=multiqc_report.html" \\
        >> ${params.outdir}/audit/audit_trail.log
    """
}
