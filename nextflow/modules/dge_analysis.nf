// ============================================================
// Module: DGE_ANALYSIS
// Calls dge_analysis.R with tximport + edgeR pipeline.
// Parameters are sourced from nextflow.config (CTC defaults).
// ============================================================

process DGE_ANALYSIS {
    tag "dge"
    label 'process_high'

    publishDir "${params.outdir}/dge", mode: 'copy'

    input:
    path quant_paths_tsv
    path samplesheet
    path gtf

    output:
    path "results/",          emit: dge_results
    path "dge_complete.flag", emit: completed_flag

    script:
    def priority = params.priority_genes ? "--priority_genes '${params.priority_genes}'" : ""
    def low_input = params.low_input_mode ? "--low_input_mode" : ""
    def use_treat  = params.use_treat     ? "--use_treat"      : ""
    """
    mkdir -p results

    Rscript /scripts/dge_analysis.R \\
        --quant_paths  ${quant_paths_tsv} \\
        --samples      ${samplesheet} \\
        --gtf          ${gtf} \\
        --outdir       results \\
        --min_cpm      ${params.min_cpm} \\
        --min_sam      ${params.min_sam} \\
        --fdr_cutoff   ${params.fdr_cutoff} \\
        --prior_count  ${params.prior_count} \\
        --lfc_treat    ${params.lfc_treat} \\
        --exclusive_min_cpm ${params.exclusive_min_cpm} \\
        --exclusive_min_sam ${params.exclusive_min_sam} \\
        --threads      ${params.threads} \\
        ${low_input} \\
        ${use_treat} \\
        ${priority}

    # Signal completion for provenance module
    echo "DGE_COMPLETE" > dge_complete.flag

    echo "[DGE] $(date -u +%Y-%m-%dT%H:%M:%SZ) status=COMPLETE outdir=results" \\
        >> ${params.outdir}/audit/audit_trail.log
    """
}
