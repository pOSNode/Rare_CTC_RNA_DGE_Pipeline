# CTC Low-Input RNA-seq Pipeline

[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-brightgreen)](https://nextflow.io)
[![Docker](https://img.shields.io/badge/Docker-enabled-blue)](https://docker.com)
[![R](https://img.shields.io/badge/R-4.3.2-276DC3)](https://r-project.org)

End-to-end differential gene expression pipeline optimised for **rare-cell circulating tumour cell (CTC) transcriptomics**. Standard RNA-seq pipelines fail for CTC data due to extremely low RNA input (<1 ng), sparse count matrices, and low mapping rates. This pipeline addresses those challenges through a series of biologically motivated design decisions.

---

## Background and Motivation

Circulating tumour cells are shed from solid tumours into peripheral blood and represent a minimally invasive window into tumour biology. Their rarity (1–10 CTCs per 10 mL blood) means that RNA extracted from CTC-enriched fractions is severely limited in quantity and quality.

Standard pipelines apply filtering thresholds (CPM > 0.5, min-sample > 3) and normalisation methods (TMM) that are calibrated for bulk tumour tissue or cell line data. Applied to CTC data, these pipelines silently discard genuine low-abundance transcripts, including clinically relevant markers such as ERBB2, ESR1, and VIM.

This pipeline was designed and validated on CTC transcriptomic data from a rare-cell liquid biopsy platform, with the following biological context in mind:

- Library sizes of 200K–2M reads (vs 20–50M in bulk RNA-seq)
- Mapping rates of 15–40% (expected given CTC purity and contaminating WBC RNA)
- Condition-exclusive expression of epithelial/mesenchymal markers defining CTC phenotype
- Multi-run samples from sequential enrichment experiments
- Highly variable cell input per sample (0 cells, 5 cells, 10 cells, positive controls) requiring cell-input-aware QC and scaling
- Batch effects from enrichment date, operator, and capture run that must be modelled additively rather than removed

---

## Pipeline Overview

```
FASTQ (paired-end, multi-run)
        │
        ▼
   [VALIDATE]  ← Sample sheet validation + audit initialisation
        │
        ▼
   [FASTP]     ← Adapter trimming, QC (conservative low-input settings)
        │
        ▼
   [SALMON]    ← Quasi-mapping quantification (GRCh38/Gencode)
        │
        ▼
   [MERGE]     ← Multi-run merging at quantification level (tximport)
        │
        ▼
   [MULTIQC]   ← Aggregated QC report
        │
        ▼
   [DGE]       ← edgeR: TMMwsp normalisation, exclusive gene retention, glmQLFTest
        │
        ▼
   [PROVENANCE] ← Audit trail, tool versions, parameter checksums
```

---

## Key Design Decisions

| Decision | Standard Pipeline | This Pipeline | Rationale |
|---|---|---|---|
| Normalisation | TMM | **TMMwsp** | TMMwsp handles zero-inflated sparse matrices better; TMM is unstable when many genes have zero counts |
| CPM filter | 0.5 | **0.1** | Preserves low-abundance CTC transcripts |
| Exclusive gene retention | Not implemented | **Enabled** | Condition-exclusive markers (e.g. epithelial genes absent in WBC fraction) are biologically meaningful even at CPM < 0.1 |
| Prior count | 2 | **0.5** | Reduced shrinkage towards mean preserves signal for lowly expressed genes |
| FDR threshold | 0.05 | **0.1** | Low replicate numbers reduce statistical power; 0.1 is appropriate given the exploratory nature of CTC biomarker work |
| Multi-run merging | FASTQ-level cat | **Quantification-level** | Preserves per-run QC traceability; tximport handles variance correctly |
| Bootstrap quantification | 20 | **50** (low-input mode) | Increased bootstraps improve TPM uncertainty estimates for sparse data |
| Batch correction | Not modelled / `removeBatchEffect` | **Additive in design matrix** | Batch is included as a covariate (`~ 0 + condition + batch`) rather than removed pre-analysis. This preserves biological variance, keeps the model honest, and avoids the artificial inflation of DE signal that comes from subtracting batch before testing |
| Cell input scaling | Not applicable | **cell_count-aware QC tiers** | Samples are annotated by input cell count (blank, 5-cell, 10-cell, positive control). QC thresholds, expected library sizes, and mapping rate warnings are scaled per tier rather than applying uniform cutoffs that would flag sparse-but-valid low-cell samples as failures |

---

## Requirements

- [Nextflow](https://nextflow.io/) >= 23.04.0
- [Docker](https://docker.com/) >= 20.10
- 8GB RAM minimum (16GB recommended)

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/ctc-rnaseq-pipeline.git
cd ctc-rnaseq-pipeline
```

### 2. Build the Docker image

```bash
docker build -t owaissiddiqi/ctc-rnaseq:1.0.0 docker/
```

### 3. Prepare your sample sheet

```csv
sample_id,condition,batch,run_id,cell_count,control_type,fastq_1,fastq_2
CTC001,tumour,batch1,CTC001_R1,10,,/data/CTC001_R1_1.fastq.gz,/data/CTC001_R1_2.fastq.gz
CTC001,tumour,batch1,CTC001_R2,10,,/data/CTC001_R2_1.fastq.gz,/data/CTC001_R2_2.fastq.gz
CTC002,tumour,batch1,CTC002_R1,5,,/data/CTC002_R1_1.fastq.gz,/data/CTC002_R1_2.fastq.gz
BLANK01,tumour,batch1,BLANK01_R1,0,blank,/data/BLANK01_R1_1.fastq.gz,/data/BLANK01_R1_2.fastq.gz
POSCTRL,tumour,batch1,POSCTRL_R1,10,positive,/data/POSCTRL_R1_1.fastq.gz,/data/POSCTRL_R1_2.fastq.gz
CTRL001,control,batch1,CTRL001_R1,10,,/data/CTRL001_R1_1.fastq.gz,/data/CTRL001_R1_2.fastq.gz
```

| Column | Required | Description |
|---|---|---|
| `sample_id` | Yes | Biological sample identifier. Multiple `run_id` rows per `sample_id` are merged at quantification level. |
| `condition` | Yes | Experimental group. Minimum 2 conditions required. |
| `batch` | No | Batch covariate. If provided, included **additively** in the design matrix (`~ 0 + condition + batch`). Never removed pre-analysis. |
| `run_id` | No | Technical run identifier. Defaults to `sample_id` if absent. |
| `cell_count` | No | Number of input cells (0, 5, 10, or any integer). Used to set per-sample QC thresholds. Blank samples should be 0. |
| `control_type` | No | `blank` or `positive`. Blank samples are QC-flagged and excluded from DGE. Positive controls are tracked separately to verify run quality. |
| `fastq_1` | Yes | Path to R1 FASTQ (gzipped). |
| `fastq_2` | Yes | Path to R2 FASTQ (gzipped). |

### 4. Run the pipeline

```bash
nextflow run nextflow/main.nf \
    --samplesheet sample_sheet.csv \
    --salmon_index /refs/salmon_index \
    --gtf /refs/gencode.v46.annotation.gtf \
    --outdir results \
    --low_input_mode
```

### 5. Run with test data

```bash
# Generate synthetic FASTQ files
python3 test_data/generate_test_fastq.py

# Build Salmon index from test reference (requires Salmon installed locally)
salmon index -t test_data/refs/transcriptome.fa -i test_data/refs/salmon_index

# Run pipeline in test profile
nextflow run nextflow/main.nf -profile test
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--low_input_mode` | `true` | Activates all low-input optimisations |
| `--min_cpm` | `0.1` | Minimum CPM for filterByExpr |
| `--fdr_cutoff` | `0.1` | FDR threshold for DE significance |
| `--lfc_treat` | `0.5` | Log fold-change threshold for glmTreat |
| `--exclusive_min_cpm` | `0.05` | CPM floor for condition-exclusive gene retention |
| `--priority_genes` | `""` | Comma-separated genes to retain regardless of expression level (e.g. `ERBB2,ESR1,VIM`) |
| `--prior_count` | `0.5` | Prior count for log-CPM transformation |
| `--threads` | `4` | Threads per process |

---

## QC Framework

CTC experiments have fundamentally heterogeneous sample inputs — a blank capture well, a 5-cell well, a 10-cell well, and a positive control cell line all look completely different by standard RNA-seq QC metrics. Applying uniform thresholds would either flag valid low-cell samples as failures or let blank contamination pass undetected. The pipeline uses a tiered QC system calibrated to expected cell input.

### Cell Input Tiers

| Tier | `cell_count` | Expected Library Size | Expected Mapping Rate | Genes Detected |
|---|---|---|---|---|
| Blank | 0 | < 50K reads | < 5% | < 200 (background only) |
| 5-cell | 5 | 200K – 800K reads | 15 – 35% | 1,000 – 4,000 |
| 10-cell | 10 | 500K – 2M reads | 20 – 40% | 3,000 – 8,000 |
| Positive control | 10 (cell line) | 1M – 5M reads | 40 – 70% | 8,000 – 15,000 |

QC warnings are issued relative to the expected range for each tier. A 5-cell sample with 300K reads and 22% mapping is **not flagged** — it is within expectation. The same metrics on a positive control cell line sample **would** be flagged as a potential capture failure.

### Blank Sample Handling

Samples where `control_type = blank` are:
1. Processed through the full pipeline (trimming, quantification, QC)
2. **Excluded from DGE analysis** automatically
3. Used to estimate background signal — any gene with CPM > 0.5 in a blank is flagged as a potential contamination signal and annotated in the results
4. Reported in `qc/blank_contamination_report.csv`

This is critical for CTC data. Genes that appear to be differentially expressed but are also detected in blank captures are candidates for WBC contamination or ambient RNA rather than genuine CTC signal.

### Positive Control Handling

Samples where `control_type = positive` are:
1. Processed and QC'd independently
2. Excluded from patient-facing DGE contrasts
3. Used to assess run-level quality: if the positive control fails expected QC thresholds, the entire run is flagged in the audit log
4. Reported in `qc/positive_control_report.csv` with pass/fail per metric

### Per-Sample QC Metrics

The pipeline generates the following QC metrics for every sample, written to `qc/sample_qc_metrics.csv` and surfaced in the MultiQC report:

| Metric | Source | Blank threshold | 5-cell threshold | 10-cell threshold |
|---|---|---|---|---|
| Total reads | fastp | Any | > 100K | > 300K |
| % reads passing QC | fastp | — | > 60% | > 70% |
| Mapping rate | Salmon | < 5% expected | 15 – 40% | 20 – 45% |
| Genes detected (CPM > 0.1) | edgeR | < 200 | > 500 | > 1,500 |
| Library size | edgeR | < 50K | > 150K | > 400K |
| Norm factor | edgeR | excluded | 0.3 – 3.0 | 0.3 – 3.0 |

Samples that fail thresholds for their tier are **not automatically removed** — they are flagged in the audit log and annotated in results, leaving the decision to the analyst.

---

## Batch Correction

### Why additive modelling, not removal

A common approach is to use `limma::removeBatchEffect` or `ComBat` to subtract batch effects from the count matrix before differential testing. This pipeline deliberately does not do this for two reasons:

1. **It inflates DE signal.** Removing batch effects before testing violates the assumption that the testing data is independent of the correction. This produces anti-conservative p-values and more false positives — a serious problem in a clinical biomarker context.

2. **It destroys variance information.** For CTC data with very few replicates per condition, batch-corrected counts can produce near-zero within-group variance, making the model overconfident.

Instead, batch is included as an **additive covariate** in the edgeR design matrix:

```r
# With batch (recommended when batch is known and balanced)
design <- model.matrix(~ 0 + condition + batch, data = samples)

# Without batch (used when batch is confounded with condition)
design <- model.matrix(~ 0 + condition, data = samples)
```

The pipeline automatically detects whether batch is confounded with condition (i.e. all samples from one condition are in the same batch). If confounding is detected, batch is dropped from the model and a warning is written to the audit log, since including a confounded batch term would make the model uninterpretable.

### Batch balance check

At pipeline start, the sample sheet validator checks whether batch is approximately balanced across conditions:

```
[WARN] Batch imbalance detected: condition 'tumour' has 80% of samples in batch1.
       Batch correction may be unreliable. Review experimental design.
```

This does not halt the pipeline but is recorded in `audit/audit_trail.log`.

## Output Structure

```
results/
├── audit/
│   ├── audit_trail.log              ← Timestamped per-process audit entries
│   ├── provenance.json              ← Machine-readable run provenance
│   ├── run_params.txt               ← Human-readable parameter dump
│   ├── nextflow_report.html         ← Nextflow execution report
│   ├── nextflow_timeline.html       ← Process timeline
│   ├── nextflow_trace.txt           ← Per-task resource usage
│   └── pipeline_dag.html            ← Pipeline DAG visualisation
├── qc/
│   ├── fastp/                       ← Per-sample fastp HTML reports
│   ├── multiqc/
│   │   └── multiqc_report.html      ← Aggregated QC report
│   ├── sample_qc_metrics.csv        ← Per-sample QC with tier-aware pass/fail flags
│   ├── blank_contamination_report.csv  ← Genes detected in blank captures
│   └── positive_control_report.csv  ← Run-level QC pass/fail from positive controls
├── quant/                           ← Per-run Salmon quant directories
├── quant_merged/
│   └── quant_paths.tsv              ← Sample → quant.sf path manifest
└── dge/
    └── YYYYMMDD_HHMMSS/
        ├── pipeline.log
        ├── *_all_genes.csv              ← Full DGE results per contrast
        ├── significant_DE_genes.csv     ← Filtered significant results
        ├── priority_genes_summary.csv
        ├── normalised_logCPM.csv
        └── sample_qc_metrics.csv
```

---

## Audit and Provenance

Every process writes a timestamped entry to `audit/audit_trail.log`:

```
[VALIDATE]   2024-11-15T09:10:01Z batch_balance=OK conditions=2 batches=2
[VALIDATE]   2024-11-15T09:10:01Z WARN batch imbalance: 'tumour' 80% in batch1 — review design
[FASTP]      2024-11-15T09:12:34Z run_id=CTC001_R1 sample=CTC001 tier=10-cell status=COMPLETE
[FASTP]      2024-11-15T09:12:40Z run_id=BLANK01_R1 sample=BLANK01 tier=blank status=COMPLETE
[SALMON]     2024-11-15T09:14:02Z run_id=CTC001_R1 sample=CTC001 tier=10-cell mapping_rate=28.4% status=PASS
[SALMON]     2024-11-15T09:14:10Z run_id=CTC002_R1 sample=CTC002 tier=5-cell mapping_rate=19.1% status=PASS
[SALMON]     2024-11-15T09:14:15Z run_id=BLANK01_R1 sample=BLANK01 tier=blank mapping_rate=2.1% status=PASS (expected for blank)
[SALMON]     2024-11-15T09:14:20Z run_id=POSCTRL_R1 sample=POSCTRL tier=positive_control mapping_rate=31.2% status=WARN (expected > 40%)
[QC]         2024-11-15T09:15:00Z blank_contamination: 47 genes detected in blank at CPM > 0.5 — annotated in DGE results
[DGE]        2024-11-15T09:22:11Z batch_in_model=TRUE design=~0+condition+batch status=COMPLETE
[PROVENANCE] 2024-11-15T09:22:15Z run_id=RUN_20241115_092215 status=COMPLETE
```

`audit/provenance.json` records tool versions, input file checksums, and all parameters used — sufficient for methods section reproduction and regulatory audit.

---

## Interpreting Low Mapping Rates

Mapping rates of 15–40% are **expected and not a failure** for CTC samples. CTCs are enriched from peripheral blood and the fraction invariably contains contaminating white blood cells (WBCs). WBC-derived reads do not map to a tumour-specific reference. If mapping rates fall below 15%, the pipeline will log a warning but continue — this may indicate sample quality issues rather than pipeline failure.

---

## Citation

If you use this pipeline, please cite:

> Siddiqi O. Rare_CTC_RNA_DGE_Pipeline: End-to-end Nextflow RNA-seq differential expression pipeline optimised for rare-cell CTC transcriptomics. Version 1.0.0. GitHub; 2026. Available from: https://github.com/pOSNode/Rare_CTC_RNA_DGE_Pipeline

---

## Author

**Owais Siddiqi** — Bioinformatics Data Scientist  
[owaissiddiqi.co.uk](https://owaissiddiqi.co.uk) | [omsiddiqi01@gmail.com](mailto:omsiddiqi01@gmail.com)
