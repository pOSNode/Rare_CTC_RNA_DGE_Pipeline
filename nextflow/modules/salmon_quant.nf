// ============================================================
// Module: SALMON_QUANT
// Quasi-mapping quantification via Salmon.
// validateMappings flag is intentionally relaxed for CTC
// samples where mapping rates are expected to be low (<30%).
// ============================================================

process SALMON_QUANT {
    tag "${meta.run_id}"
    label 'process_medium'

    publishDir "${params.outdir}/quant/${meta.sample_id}/${meta.run_id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("${meta.run_id}_quant"), emit: quant_dir
    path "${meta.run_id}_salmon.log",              emit: log

    script:
    def (fq1, fq2) = reads
    // Low-input: increase bootstrap samples for uncertainty quantification
    def bootstraps = params.low_input_mode ? "--numBootstraps 50" : "--numBootstraps 20"
    """
    salmon quant \\
        --index        ${index} \\
        --libType      A \\
        --mates1       ${fq1} \\
        --mates2       ${fq2} \\
        --output       ${meta.run_id}_quant \\
        --threads      ${params.threads} \\
        --gcBias \\
        --seqBias \\
        --posBias \\
        ${bootstraps} \\
        --validateMappings \\
        --minScoreFraction 0.65 \\
        2>&1 | tee ${meta.run_id}_salmon.log

    # Flag low mapping rates — expected in CTC but log for audit
    MAPPING_RATE=\$(grep "Mapping rate" ${meta.run_id}_salmon.log | \\
        grep -oP '[0-9]+\\.[0-9]+' | tail -1)

    echo "[SALMON] $(date -u +%Y-%m-%dT%H:%M:%SZ) run_id=${meta.run_id} sample=${meta.sample_id} mapping_rate=\${MAPPING_RATE}% status=COMPLETE" \\
        >> ${params.outdir}/audit/audit_trail.log

    if (( \$(echo "\${MAPPING_RATE} < 20" | bc -l) )); then
        echo "[WARN] $(date -u +%Y-%m-%dT%H:%M:%SZ) run_id=${meta.run_id} mapping_rate=\${MAPPING_RATE}% is below 20% - review sample quality" \\
            >> ${params.outdir}/audit/audit_trail.log
    fi
    """
}
