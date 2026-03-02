#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================
// CTC Low-Input RNA-seq Pipeline
// Author: Owais Siddiqi
// Description: End-to-end DGE pipeline optimised for rare-cell
//              circulating tumour cell (CTC) transcriptomics.
//              Handles paired-end FASTQ, multi-run merging,
//              low-input aware QC, quantification and DGE.
// ============================================================

include { VALIDATE_SAMPLESHEET  } from './modules/validate_samplesheet'
include { FASTP                  } from './modules/fastp'
include { SALMON_QUANT           } from './modules/salmon_quant'
include { MERGE_RUNS             } from './modules/merge_runs'
include { MULTIQC                } from './modules/multiqc'
include { DGE_ANALYSIS           } from './modules/dge_analysis'
include { WRITE_PROVENANCE       } from './modules/provenance'

// ------------------------------------------------------------
// Parameter validation
// ------------------------------------------------------------
def required_params = ['samplesheet', 'salmon_index', 'gtf', 'outdir']
required_params.each { p ->
    if (!params[p]) {
        error "Missing required parameter: --${p}\nRun with --help for usage."
    }
}

if (params.help) {
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║         CTC Low-Input RNA-seq Pipeline v${params.pipeline_version}          ║
    ╚══════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf [options]

    Required:
        --samplesheet       Path to sample sheet CSV
        --salmon_index      Path to Salmon index directory
        --gtf               Path to GTF annotation file (GRCh38/Gencode)
        --outdir            Output directory

    Optional:
        --low_input_mode    Enable low-input RNA optimisations (default: ${params.low_input_mode})
        --fdr_cutoff        FDR threshold for DGE (default: ${params.fdr_cutoff})
        --min_cpm           Minimum CPM filter (default: ${params.min_cpm})
        --lfc_treat         Log fold-change threshold (default: ${params.lfc_treat})
        --priority_genes    Comma-separated gene list (e.g. ERBB2,ESR1)
        --threads           Threads per process (default: ${params.threads})

    Sample sheet format (CSV):
        sample_id,condition,batch,run_id,fastq_1,fastq_2

    Full documentation: docs/methods.md
    """.stripIndent()
    exit 0
}

// ------------------------------------------------------------
// Workflow
// ------------------------------------------------------------
workflow {

    // -- Audit: record pipeline start
    run_id   = "RUN_" + new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
    log.info "Pipeline run ID : ${run_id}"
    log.info "Pipeline version: ${params.pipeline_version}"
    log.info "Low-input mode  : ${params.low_input_mode}"

    // -- 1. Validate & parse sample sheet
    ch_samplesheet = VALIDATE_SAMPLESHEET(
        Channel.fromPath(params.samplesheet, checkIfExists: true)
    )

    // -- 2. Build per-run input channel
    //    Each row: [meta, [fastq_1, fastq_2]]
    //    meta map keys: sample_id, condition, batch, run_id
    ch_reads = ch_samplesheet
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                sample_id : row.sample_id,
                condition : row.condition,
                batch     : row.batch     ?: "batch1",
                run_id    : row.run_id    ?: row.sample_id
            ]
            def fq1 = file(row.fastq_1, checkIfExists: true)
            def fq2 = file(row.fastq_2, checkIfExists: true)
            [ meta, [ fq1, fq2 ] ]
        }

    // -- 3. QC + trimming (per run)
    FASTP(ch_reads)

    // -- 4. Salmon quantification (per run)
    SALMON_QUANT(
        FASTP.out.trimmed_reads,
        Channel.fromPath(params.salmon_index, checkIfExists: true).collect()
    )

    // -- 5. Merge multi-run quants per sample_id
    //    Groups by sample_id; tximport handles list of quant dirs
    ch_quants_grouped = SALMON_QUANT.out.quant_dir
        .map { meta, dir -> [ meta.sample_id, meta, dir ] }
        .groupTuple(by: 0)
        .map { sample_id, metas, dirs ->
            def merged_meta = metas[0]  // representative meta for sample
            [ merged_meta, dirs ]
        }

    MERGE_RUNS(ch_quants_grouped)

    // -- 6. MultiQC
    ch_multiqc_files = FASTP.out.json_report
        .mix(SALMON_QUANT.out.log)
        .collect()

    MULTIQC(ch_multiqc_files)

    // -- 7. DGE analysis
    DGE_ANALYSIS(
        MERGE_RUNS.out.quant_paths_tsv,
        Channel.fromPath(params.samplesheet, checkIfExists: true),
        Channel.fromPath(params.gtf,         checkIfExists: true)
    )

    // -- 8. Write provenance / audit trail
    WRITE_PROVENANCE(
        run_id,
        DGE_ANALYSIS.out.completed_flag
    )
}

workflow.onComplete {
    log.info """
    ══════════════════════════════════════════
    Pipeline completed : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Run ID             : ${workflow.runName}
    Duration           : ${workflow.duration}
    Output directory   : ${params.outdir}
    ══════════════════════════════════════════
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
}
