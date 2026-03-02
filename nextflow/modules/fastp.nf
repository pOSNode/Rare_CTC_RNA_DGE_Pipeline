// ============================================================
// Module: FASTP
// Paired-end QC and adapter trimming optimised for low-input
// CTC samples. Conservative trimming to preserve rare reads.
// ============================================================

process FASTP {
    tag "${meta.run_id}"
    label 'process_medium'

    publishDir "${params.outdir}/qc/fastp/${meta.sample_id}", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.html') ? fn : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.run_id}_{1,2}.trimmed.fastq.gz"), emit: trimmed_reads
    path "${meta.run_id}_fastp.json",                                emit: json_report
    path "${meta.run_id}_fastp.html",                                emit: html_report

    script:
    def (fq1, fq2) = reads
    def extra_args  = params.low_input_mode ? "--length_required 20 --low_complexity_filter" : ""
    """
    fastp \\
        --in1  ${fq1} \\
        --in2  ${fq2} \\
        --out1 ${meta.run_id}_1.trimmed.fastq.gz \\
        --out2 ${meta.run_id}_2.trimmed.fastq.gz \\
        --json ${meta.run_id}_fastp.json \\
        --html ${meta.run_id}_fastp.html \\
        --thread ${params.threads} \\
        --detect_adapter_for_pe \\
        --correction \\
        --overrepresentation_analysis \\
        --qualified_quality_phred 15 \\
        --unqualified_percent_limit 40 \\
        ${extra_args} \\
        2>&1 | tee ${meta.run_id}_fastp.log

    # Audit entry
    echo "[FASTP] $(date -u +%Y-%m-%dT%H:%M:%SZ) run_id=${meta.run_id} sample=${meta.sample_id} status=COMPLETE" \\
        >> ${params.outdir}/audit/audit_trail.log
    """
}
