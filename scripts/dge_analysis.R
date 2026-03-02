#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(tidyverse)
  library(tximport)
  library(edgeR)
  library(limma)
})

# ============================================================
# CTC Low-Input RNA-seq DGE Analysis
# Author: Owais Siddiqi
# Description: edgeR-based differential gene expression for
#              low-input rare-cell CTC transcriptomics.
#              Uses tximport for multi-run sample merging,
#              TMMwsp normalisation for sparse count data,
#              condition-exclusive gene retention to preserve
#              rare transcript signals, length-based offsets,
#              and a comprehensive QC report with spike-level
#              aware heatmaps for linearity assessment.
# ============================================================

# ---------- CLI ----------
option_list <- list(
  make_option("--quant_dir",  type="character", default="quant",
              help="Directory containing Salmon quantification folders"),
  make_option("--gtf",        type="character", default="ref/gencode.v46.annotation.gtf",
              help="GTF annotation file (GRCh38/Gencode)"),
  make_option("--samples",    type="character", default="sample_table.csv",
              help="Sample metadata CSV (sample_id, condition, batch, cell_count, control_type)"),
  make_option("--outdir",     type="character", default="results",
              help="Output directory"),
  # Filtering — relaxed for low-input CTC
  make_option("--min_cpm",    type="double",  default=0.1,
              help="Minimum CPM for filterByExpr [default: 0.1]"),
  make_option("--min_sam",    type="integer", default=2,
              help="Minimum samples for filterByExpr [default: 2]"),
  make_option("--low_input_mode", action="store_true", default=FALSE,
              help="Enable low-input RNA specific optimisations"),
  # Condition-exclusive retention for rare markers
  make_option("--keep_exclusive",      type="logical", default=TRUE,
              help="Keep genes exclusive to one condition"),
  make_option("--exclusive_min_cpm",   type="double",  default=0.05,
              help="Min CPM for exclusive genes [default: 0.05]"),
  make_option("--exclusive_min_sam",   type="integer", default=1,
              help="Min samples for exclusive genes [default: 1]"),
  # Statistical testing
  make_option("--lfc_treat",   type="double",  default=0.5,
              help="Log-fold change threshold for glmTreat [default: 0.5]"),
  make_option("--use_treat",   action="store_true", default=FALSE,
              help="Use glmTreat instead of glmQLFTest"),
  make_option("--fdr_cutoff",  type="double",  default=0.1,
              help="FDR cutoff for significance [default: 0.1]"),
  make_option("--prior_count", type="double",  default=0.5,
              help="Prior count for log transformation [default: 0.5]"),
  # Gene-specific
  make_option("--priority_genes", type="character", default="",
              help="Comma-separated genes to retain regardless of expression level"),
  make_option("--species",     type="character", default="human",
              help="Species: human or mouse [default: human]"),
  make_option("--threads",     type="integer",   default=1,
              help="Number of threads for parallel processing")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "Enhanced DGE pipeline optimised for low-input CTC RNA-seq"
))

# ---------- Output directory + logging ----------
timestamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")
outdir_run <- file.path(opt$outdir, timestamp)
dir.create(outdir_run, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(outdir_run, "pipeline.log")
log_msg  <- function(msg, ...) {
  full_msg <- sprintf(paste0("[%s] ", msg), format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"), ...)
  message(full_msg)
  cat(full_msg, "\n", file = log_file, append = TRUE)
}

log_msg("=== CTC DGE Pipeline Start ===")
log_msg("Output directory : %s", outdir_run)
log_msg("Samples file     : %s", opt$samples)
log_msg("GTF              : %s", opt$gtf)
log_msg("FDR cutoff       : %.2f", opt$fdr_cutoff)
log_msg("LFC threshold    : %.2f", opt$lfc_treat)
log_msg("Low-input mode   : %s", ifelse(opt$low_input_mode, "ENABLED", "disabled"))

# ---------- Low-input mode adjustments ----------
if (opt$low_input_mode) {
  opt$min_cpm           <- min(opt$min_cpm,           0.1)
  opt$exclusive_min_cpm <- min(opt$exclusive_min_cpm, 0.05)
  opt$exclusive_min_sam <- 1
  opt$prior_count       <- min(opt$prior_count,       0.5)
  opt$fdr_cutoff        <- max(opt$fdr_cutoff,        0.1)
  log_msg("Low-input adjustments: min_cpm=%.2f fdr=%.2f prior=%.2f",
          opt$min_cpm, opt$fdr_cutoff, opt$prior_count)
}

# ---------- Priority genes ----------
priority_genes <- if (nzchar(opt$priority_genes)) {
  trimws(unlist(strsplit(opt$priority_genes, ",")))
} else { character(0) }
if (length(priority_genes) > 0)
  log_msg("Priority genes: %s", paste(priority_genes, collapse=", "))

# ---------- Load annotation ----------
load_or_install <- function(pkgs) {
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
  if (length(to_install)) {
    log_msg("Installing missing packages: %s", paste(to_install, collapse=", "))
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager", repos="https://cloud.r-project.org")
    BiocManager::install(to_install, ask=FALSE, update=FALSE)
  }
  invisible(lapply(pkgs, require, character.only=TRUE))
}

if (tolower(opt$species) == "mouse") {
  load_or_install(c("AnnotationDbi","org.Mm.eg.db","pheatmap","ggrepel","RColorBrewer"))
  org_db <- org.Mm.eg.db::org.Mm.eg.db
  log_msg("Using mouse annotation database")
} else {
  load_or_install(c("AnnotationDbi","org.Hs.eg.db","pheatmap","ggrepel","RColorBrewer"))
  org_db <- org.Hs.eg.db::org.Hs.eg.db
  log_msg("Using human annotation database")
}

suppressPackageStartupMessages({
  library(AnnotationDbi); library(pheatmap)
  library(ggplot2); library(ggrepel); library(RColorBrewer)
})
select <- dplyr::select

# ---------- Breast / CTC marker set ----------
breast_marker_primary <- c("ESR1","PGR","ERBB2","MKI67","EPCAM","KRT19","FN1","PTPRC")
breast_marker_aliases <- c("HER2","KI67","EpCAM")
breast_marker_all     <- unique(toupper(c(breast_marker_primary, breast_marker_aliases)))
alias_to_canonical    <- c("HER2"="ERBB2","KI67"="MKI67","EPCAM"="EPCAM")

resolve_to_symbols <- function(keys_vec) {
  keys_vec <- unique(toupper(keys_vec))
  sym1 <- suppressMessages(AnnotationDbi::mapIds(org_db, keys=keys_vec, keytype="SYMBOL",
                                                  column="SYMBOL", multiVals="first"))
  sym2 <- suppressMessages(AnnotationDbi::mapIds(org_db, keys=keys_vec, keytype="ALIAS",
                                                  column="SYMBOL", multiVals="first"))
  unique(na.omit(c(unname(sym1), unname(sym2))))
}
breast_marker_symbols <- resolve_to_symbols(breast_marker_all)

# ---------- Helpers ----------
parse_cell_count <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("blank|control|ctrl|vehicle|untreated|none|zero|dmso", "0", x)
  as.numeric(stringr::str_extract(x, "\\d+\\.?\\d*"))
}

format_spike_labels <- function(spike_levels) {
  vapply(spike_levels, function(x) {
    if (is.na(x) || !is.finite(x)) "Unknown"
    else if (x == 0)                "Blank\n(0 cells)"
    else                            sprintf("%.0f cells", x)
  }, character(1))
}

create_sample_labels <- function(samples_df, include_batch=FALSE) {
  if (include_batch && nlevels(samples_df$batch) > 1) {
    paste0(samples_df$sample_id, "\n(", samples_df$condition, ", ", samples_df$batch, ")")
  } else {
    paste0(samples_df$sample_id, "\n(", samples_df$condition, ")")
  }
}

safe_discrete_cols <- function(n, base="Set1") {
  max_n <- RColorBrewer::brewer.pal.info[base, "maxcolors"]
  if (is.na(max_n)) max_n <- 9
  if (n <= max_n) brewer.pal(max(3, n), base)[seq_len(n)]
  else colorRampPalette(brewer.pal(max_n, base))(n)
}

annotate_genes <- function(df, gene_col="gene_id") {
  if (nrow(df) == 0) return(df)
  ens_ids <- sub("\\.\\d+$", "", df[[gene_col]])
  suppressMessages({
    symbols    <- AnnotationDbi::mapIds(org_db, keys=ens_ids, keytype="ENSEMBL",
                                        column="SYMBOL",   multiVals="first")
    names_full <- AnnotationDbi::mapIds(org_db, keys=ens_ids, keytype="ENSEMBL",
                                        column="GENENAME", multiVals="first")
    entrez     <- AnnotationDbi::mapIds(org_db, keys=ens_ids, keytype="ENSEMBL",
                                        column="ENTREZID", multiVals="first")
    alias_list <- AnnotationDbi::mapIds(org_db, keys=unname(entrez), keytype="ENTREZID",
                                        column="ALIAS",
                                        multiVals=function(x) paste(unique(x), collapse=","))
  })
  df$hgnc_symbol   <- unname(symbols)
  df$gene_name     <- unname(names_full)
  df$entrez_id     <- unname(entrez)
  df$alias_symbols <- unname(alias_list)
  df$symbol_search <- tolower(ifelse(
    is.na(df$hgnc_symbol) | df$hgnc_symbol == "",
    df$alias_symbols,
    ifelse(is.na(df$alias_symbols) | df$alias_symbols == "",
           df$hgnc_symbol,
           paste(df$hgnc_symbol, df$alias_symbols, sep=","))
  ))
  df
}

parse_contrast_name <- function(nm) {
  if      (grepl("_vs_", nm)) { parts <- strsplit(nm,"_vs_")[[1]]; list(A=parts[1],B=parts[2]) }
  else if (grepl(" vs ",  nm)) { parts <- strsplit(nm," vs ")[[1]];  list(A=parts[1],B=parts[2]) }
  else                          { list(A=nm, B="reference") }
}

.ens_to_symbol <- function(ens_ids, org_db) {
  ens_clean <- sub("\\.\\d+$","", ens_ids)
  syms <- suppressMessages(AnnotationDbi::mapIds(org_db, keys=ens_clean, keytype="ENSEMBL",
                                                  column="SYMBOL", multiVals="first"))
  ifelse(is.na(syms) | syms=="", ens_ids, unname(syms))
}

# ---------- Build / cache tx2gene from GTF ----------
tx2gene_path <- file.path(opt$outdir, "tx2gene.rds")
gtf_mtime    <- file.info(opt$gtf)$mtime

build_tx2gene <- function() {
  log_msg("Building tx2gene from GTF...")
  gtf_lines <- system(sprintf("grep -c -v '^#' %s", shQuote(opt$gtf)), intern=TRUE)
  log_msg("  Processing %s GTF lines...", gtf_lines)
  gtf <- fread(
    cmd=sprintf("grep -v '^#' %s | grep '\ttranscript\t'", shQuote(opt$gtf)),
    sep="\t", header=FALSE,
    col.names=c("chr","src","type","start","end","score","strand","phase","attr"),
    showProgress=FALSE
  )
  extract_attr <- function(attr_str, key) {
    stringr::str_match(attr_str, sprintf('%s "([^"]+)"', key))[,2]
  }
  tx2gene <- data.frame(
    transcript = extract_attr(gtf$attr, "transcript_id"),
    gene       = extract_attr(gtf$attr, "gene_id"),
    stringsAsFactors=FALSE
  ) %>%
    mutate(
      transcript = sub("\\|.*$","",  transcript),
      transcript = sub("\\.\\d+$","",transcript),
      gene       = sub("\\.\\d+$","",gene)
    ) %>%
    distinct() %>%
    filter(!is.na(transcript), !is.na(gene))
  saveRDS(list(tx2gene=tx2gene, gtf_path=opt$gtf, gtf_mtime=gtf_mtime, created=Sys.time()),
          tx2gene_path)
  log_msg("  Created tx2gene with %d mappings", nrow(tx2gene))
  tx2gene
}

if (file.exists(tx2gene_path)) {
  cached <- readRDS(tx2gene_path)
  if (!identical(cached$gtf_mtime, gtf_mtime)) {
    log_msg("GTF changed — rebuilding tx2gene...")
    tx2gene <- build_tx2gene()
  } else {
    tx2gene <- cached$tx2gene
    log_msg("Loaded cached tx2gene (%d mappings)", nrow(tx2gene))
  }
} else {
  tx2gene <- build_tx2gene()
}

# ---------- Sample metadata + validation ----------
log_msg("Loading sample metadata...")
samples <- readr::read_csv(opt$samples, show_col_types=FALSE)

required_cols <- c("sample_id","condition")
if (!all(required_cols %in% names(samples)))
  stop("Sample table must contain: ", paste(required_cols, collapse=", "))

if (!"batch"        %in% names(samples)) { samples$batch        <- "batch1";      log_msg("No batch column — assuming single batch") }
if (!"cell_count"   %in% names(samples))   samples$cell_count   <- NA_real_
if (!"control_type" %in% names(samples))   samples$control_type <- NA_character_

samples$sample_id <- trimws(as.character(samples$sample_id))
samples$condition <- factor(samples$condition)
samples$batch     <- factor(samples$batch)

if (any(duplicated(samples$sample_id)))
  stop("Duplicate sample IDs: ", paste(samples$sample_id[duplicated(samples$sample_id)], collapse=", "))

conds_tab <- table(samples$condition)
log_msg("Design: %s", paste(names(conds_tab), conds_tab, sep="=", collapse=", "))
if (length(conds_tab) < 2L) stop("Need ≥2 conditions; found: ", paste(names(conds_tab), collapse=", "))

low_rep      <- names(conds_tab)[conds_tab < 2L]
low_rep_warn <- names(conds_tab)[conds_tab < 3L]
if (length(low_rep))      stop("Conditions with <2 replicates: ", paste(low_rep, collapse=", "))
if (length(low_rep_warn)) log_msg("WARNING: Low replication (<3): %s", paste(low_rep_warn, collapse=", "))

# Batch balance check
if (nlevels(samples$batch) > 1) {
  bct <- table(samples$condition, samples$batch)
  if (sum(bct == 0) > 0)
    log_msg("WARNING: Batch imbalance — %d empty condition×batch cells. Review design.", sum(bct==0))
}

# Identify control samples
blank_idx   <- !is.na(samples$control_type) & tolower(samples$control_type) == "blank"
posctrl_idx <- !is.na(samples$control_type) & tolower(samples$control_type) == "positive"
dge_idx     <- !blank_idx & !posctrl_idx
log_msg("Samples: %d total | %d blanks | %d positive controls | %d DGE",
        nrow(samples), sum(blank_idx), sum(posctrl_idx), sum(dge_idx))

# ---------- Locate quant files ----------
sfns <- file.path(path.expand(opt$quant_dir), samples$sample_id, "quant.sf")
names(sfns) <- samples$sample_id
missing_files <- samples$sample_id[!file.exists(sfns)]
if (length(missing_files)) stop("Missing quant.sf: ", paste(missing_files, collapse=", "))
write.csv(samples, file.path(outdir_run,"sample_table_processed.csv"), row.names=FALSE)

# ---------- tximport ----------
log_msg("Importing Salmon quantification files...")
txi <- tryCatch(
  tximport(sfns, type="salmon", tx2gene=tx2gene,
           countsFromAbundance="no", ignoreAfterBar=TRUE, ignoreTxVersion=TRUE),
  error=function(e) { log_msg("ERROR in tximport: %s", e$message); stop("Import failed.") }
)
log_msg("tximport: %d genes x %d samples", nrow(txi$counts), ncol(txi$counts))

# ---------- Length-based offsets ----------
log_msg("Computing length-based offsets...")
len <- txi$length
if (is.null(len)) {
  log_msg("WARNING: No length info — skipping length correction")
  offset_mat <- matrix(0, nrow=nrow(txi$counts), ncol=ncol(txi$counts),
                       dimnames=dimnames(txi$counts))
} else {
  len[len <= 0 | !is.finite(len)] <- NA
  if (any(rowSums(is.na(len)) == ncol(len)))
    log_msg("WARNING: %d genes have no valid length info",
            sum(rowSums(is.na(len)) == ncol(len)))
  geomean <- exp(rowMeans(log(len), na.rm=TRUE))
  med_len <- median(geomean[is.finite(geomean)], na.rm=TRUE)
  geomean[!is.finite(geomean)] <- med_len
  offset_mat <- log(sweep(len, 1, geomean, "/"))
  offset_mat[!is.finite(offset_mat)] <- 0
  rownames(offset_mat) <- rownames(txi$counts)
  colnames(offset_mat) <- colnames(txi$counts)
}

# ---------- Align samples to count matrix ----------
samples <- samples[match(colnames(txi$counts), samples$sample_id), ]
stopifnot(all(colnames(txi$counts) == samples$sample_id))

samples_dge <- samples[dge_idx, ]
counts_dge  <- round(txi$counts)[, dge_idx, drop=FALSE]

# ---------- DGEList + filtering ----------
log_msg("Creating DGEList...")
y <- DGEList(counts=counts_dge)

# Batch confounding check
batch_confounded <- FALSE
if (nlevels(samples_dge$batch) > 1) {
  cc <- table(samples_dge$condition, samples_dge$batch)
  if (any(rowSums(cc > 0) == 1)) {
    log_msg("WARNING: Batch confounded with condition — dropping batch from model")
    batch_confounded <- TRUE
  }
}

design_base <- if (!batch_confounded && nlevels(samples_dge$batch) > 1) {
  log_msg("Design: ~ 0 + condition + batch")
  model.matrix(~ 0 + condition + batch, data=samples_dge)
} else {
  log_msg("Design: ~ 0 + condition")
  model.matrix(~ 0 + condition, data=samples_dge)
}

# Base filtering (relaxed for low-input)
keep_base <- filterByExpr(y, design=design_base,
                          min.count=opt$min_cpm,
                          min.total.count=opt$min_sam,
                          min.prop=0.5)

# Condition-exclusive gene retention for CTC rare markers
if (opt$keep_exclusive) {
  log_msg("Exclusive gene detection...")
  cpm_raw        <- cpm(y, log=FALSE)
  keep_exclusive <- rep(FALSE, nrow(y))
  for (cond in levels(samples_dge$condition)) {
    cond_samp  <- which(samples_dge$condition == cond)
    if (length(cond_samp) >= opt$exclusive_min_sam) {
      present_in  <- rowSums(cpm_raw[, cond_samp, drop=FALSE] >= opt$exclusive_min_cpm) >= opt$exclusive_min_sam
      other_samp  <- which(samples_dge$condition != cond)
      absent_else <- if (length(other_samp) > 0)
                       rowSums(cpm_raw[, other_samp, drop=FALSE] >= (opt$exclusive_min_cpm * 2)) == 0
                     else rep(TRUE, nrow(cpm_raw))
      excl <- present_in & absent_else
      keep_exclusive <- keep_exclusive | excl
      if (sum(excl) > 0) log_msg("  %d genes exclusive to '%s'", sum(excl), cond)
    }
  }
  if (length(priority_genes) > 0) {
    log_msg("Priority gene retention...")
    prio_sym       <- resolve_to_symbols(priority_genes)
    gene_ids_clean <- sub("\\.\\d+$","", rownames(y))
    ens_ids <- suppressMessages(AnnotationDbi::mapIds(org_db, keys=prio_sym, keytype="SYMBOL",
                                                      column="ENSEMBL", multiVals="first"))
    ens_ids   <- sub("\\.\\d+$","", unname(ens_ids))
    prio_rows <- which(gene_ids_clean %in% ens_ids)
    if (length(prio_rows) > 0) {
      keep_prio <- rowSums(cpm_raw[prio_rows, , drop=FALSE] >= 0.01) >= 1
      keep_exclusive[prio_rows[keep_prio]] <- TRUE
      log_msg("  Retained %d priority genes", sum(keep_prio))
    }
  }
  keep <- keep_base | keep_exclusive
  log_msg("Post-filter: %d genes (%d exclusive additions)",
          sum(keep), sum(keep_exclusive & !keep_base))
} else {
  keep <- keep_base
  log_msg("Post-filter: %d genes", sum(keep))
}
y <- y[keep,, keep.lib.sizes=FALSE]

# ---------- Blank contamination report ----------
if (any(blank_idx)) {
  log_msg("Blank contamination check...")
  blank_cpm  <- cpm(round(txi$counts)[, blank_idx, drop=FALSE], log=FALSE)
  blank_gens <- rownames(blank_cpm)[rowSums(blank_cpm >= 0.5) >= 1]
  blank_rep  <- data.frame(gene_id=blank_gens,
                           max_blank_cpm=apply(blank_cpm[blank_gens,,drop=FALSE],1,max),
                           in_dge_results=blank_gens %in% rownames(y))
  blank_rep  <- annotate_genes(blank_rep)
  write.csv(blank_rep, file.path(dirname(outdir_run),"blank_contamination_report.csv"), row.names=FALSE)
  log_msg("%d genes detected at CPM > 0.5 in blanks (annotated in results)", nrow(blank_rep))
}

# ---------- Positive control QC ----------
if (any(posctrl_idx)) {
  log_msg("Positive control QC...")
  pc_counts <- round(txi$counts)[, posctrl_idx, drop=FALSE]
  pc_ls     <- colSums(pc_counts)
  pc_rep    <- data.frame(
    sample_id      = samples$sample_id[posctrl_idx],
    lib_size       = pc_ls,
    genes_detected = colSums(cpm(pc_counts,log=FALSE) > 0.1),
    pass_lib_size  = pc_ls > 1e6,
    pass_genes     = colSums(cpm(pc_counts,log=FALSE) > 0.1) > 8000
  )
  write.csv(pc_rep, file.path(dirname(outdir_run),"positive_control_report.csv"), row.names=FALSE)
  if (any(!pc_rep$pass_lib_size | !pc_rep$pass_genes))
    log_msg("WARNING: Positive control(s) failed QC — review run quality")
}

# ---------- TMMwsp normalisation ----------
log_msg("TMMwsp normalisation...")
y <- calcNormFactors(y, method="TMMwsp")
nf_range <- range(y$samples$norm.factors)
if (nf_range[2] / nf_range[1] > 5) {
  log_msg("WARNING: Extreme norm factor range (%.2f - %.2f) — evaluating UQ", nf_range[1], nf_range[2])
  y_uq    <- calcNormFactors(y, method="upperquartile")
  uq_range <- range(y_uq$samples$norm.factors)
  if ((uq_range[2]/uq_range[1]) < (nf_range[2]/nf_range[1])) {
    log_msg("Switching to upper-quartile normalisation")
    y <- y_uq
  }
}
log_msg("Norm factors: %.3f - %.3f", nf_range[1], nf_range[2])

# ---------- Design matrix ----------
design <- if (!batch_confounded && nlevels(samples_dge$batch) > 1)
  model.matrix(~ 0 + condition + batch, data=samples_dge)
else
  model.matrix(~ 0 + condition, data=samples_dge)
colnames(design) <- make.names(colnames(design))
write.csv(design, file.path(outdir_run,"design_matrix.csv"))

# ---------- Apply length-based offsets ----------
offset_aligned <- offset_mat[rownames(y), colnames(y), drop=FALSE]
lib_factor      <- log((y$samples$norm.factors * y$samples$lib.size) / mean(y$samples$lib.size))
y$offset        <- sweep(offset_aligned, 2, lib_factor, "+")
target          <- log(y$samples$lib.size / mean(y$samples$lib.size))
delta_off       <- colMeans(y$offset) - target
y$offset        <- sweep(y$offset, 2, delta_off, "-")

# ---------- Dispersion estimation ----------
log_msg("Estimating dispersions (robust=TRUE)...")
y <- estimateDisp(y, design, robust=TRUE)
if (opt$low_input_mode && median(y$samples$lib.size) < 5e5) {
  log_msg("Very low library sizes (median %.0f) — increased prior.df", median(y$samples$lib.size))
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design, prior.df=20)
}
log_msg("Common dispersion  : %.4f", y$common.dispersion)
log_msg("Trended range      : %.4f - %.4f",
        range(y$trended.dispersion)[1], range(y$trended.dispersion)[2])

# ---------- Model fit ----------
log_msg("Fitting QL model (robust=TRUE)...")
fit <- glmQLFit(y, design, robust=TRUE)

# QL dispersion outlier shrinkage for low-input data
# Floors var.post at 10% of prior to prevent noisy genes dominating variance estimates
if (opt$low_input_mode) {
  ql_range <- range(fit$var.post)
  log_msg("QL dispersion range: %.4f - %.4f", ql_range[1], ql_range[2])
  if (ql_range[2] / ql_range[1] > 100) {
    log_msg("Large QL range — applying shrinkage floor (10%% of prior)")
    fit$var.post <- pmax(fit$var.post, fit$var.prior * 0.1)
    log_msg("Post-shrinkage QL: %.4f - %.4f", min(fit$var.post), max(fit$var.post))
  }
}

# ---------- Normalised expression matrices ----------
# Computed once; referenced throughout remainder of script
log_msg("Computing normalised expression matrices...")
y_cpm       <- y; y_cpm$offset <- NULL
cpm_norm    <- cpm(y_cpm, normalized.lib.sizes=TRUE, log=FALSE)
logcpm_norm <- cpm(y_cpm, normalized.lib.sizes=TRUE, log=TRUE, prior.count=opt$prior_count)

# Raw (unnormalised) CPM for detection rate calculations across all samples
cnt_all    <- round(txi$counts)[, dge_idx, drop=FALSE]
lib_all    <- colSums(cnt_all)
cpm_all    <- cpm(cnt_all, log=FALSE, normalized.lib.sizes=FALSE, lib.size=lib_all)
logcpm_all <- cpm(cnt_all, log=TRUE,  prior.count=2,
                  normalized.lib.sizes=FALSE, lib.size=lib_all)

# Heatmap base matrix (blank-centred if reference identified, else global z-score)
ref_patterns <- c("control","ctrl","blank","wt","wildtype","untreated","dmso","vehicle")
ref_cond     <- NULL
for (pattern in ref_patterns) {
  m <- grep(pattern, levels(samples_dge$condition), ignore.case=TRUE, value=TRUE)
  if (length(m) > 0) { ref_cond <- m[1]; break }
}
H_base <- if (!is.null(ref_cond)) {
  ref_idx <- which(samples_dge$condition == ref_cond)
  sweep(logcpm_norm, 1, rowMeans(logcpm_norm[, ref_idx, drop=FALSE]), "-")
} else {
  t(scale(t(logcpm_norm)))
}
zlim <- c(-2.5, 2.5)
H_base[H_base < zlim[1]] <- zlim[1]
H_base[H_base > zlim[2]] <- zlim[2]

# ---------- Contrasts ----------
log_msg("Setting up contrasts...")
conds         <- levels(samples_dge$condition)
contrast_list <- list()
ref_cond_des  <- NULL
for (pattern in ref_patterns) {
  m <- grep(pattern, conds, ignore.case=TRUE, value=TRUE)
  if (length(m) > 0) { ref_cond_des <- m[1]; break }
}
if (!is.null(ref_cond_des)) {
  log_msg("Reference condition: %s", ref_cond_des)
  for (c in setdiff(conds, ref_cond_des)) {
    c1 <- make.names(paste0("condition",c)); c0 <- make.names(paste0("condition",ref_cond_des))
    nm <- paste(c,"vs",ref_cond_des,sep="_")
    contrast_list[[nm]] <- limma::makeContrasts(contrasts=paste0(c1," - ",c0), levels=design)
  }
} else {
  log_msg("No reference — all pairwise comparisons")
  for (i in 1:(length(conds)-1)) for (j in (i+1):length(conds)) {
    c1 <- make.names(paste0("condition",conds[i])); c2 <- make.names(paste0("condition",conds[j]))
    nm <- paste(conds[i],"vs",conds[j],sep="_")
    contrast_list[[nm]] <- limma::makeContrasts(contrasts=paste0(c1," - ",c2), levels=design)
  }
}
log_msg("%d contrast(s)", length(contrast_list))

# ---------- DGE testing ----------
log_msg("Differential expression analysis...")
res_all <- list()

for (contrast_name in names(contrast_list)) {
  log_msg("  Testing: %s", contrast_name)
  contr <- contrast_list[[contrast_name]]
  qres  <- if (opt$use_treat) {
    log_msg("    glmTreat (LFC=%.2f)", opt$lfc_treat)
    glmTreat(fit, contrast=contr, lfc=opt$lfc_treat)
  } else {
    glmQLFTest(fit, contrast=contr)
  }
  tt <- topTags(qres, n=Inf, adjust.method="BH", sort.by="PValue")$table %>%
    rownames_to_column("gene_id") %>%
    mutate(gene_id_clean=sub("\\.\\d+$","",gene_id), contrast=contrast_name)

  AB     <- parse_contrast_name(contrast_name)
  A_cols <- which(samples_dge$condition == AB$A)
  B_cols <- which(samples_dge$condition == AB$B)

  tt$meanCPM_A    <- if (length(A_cols)>0) rowMeans(cpm_norm[match(tt$gene_id,rownames(cpm_norm)),    A_cols,drop=FALSE]) else NA
  tt$meanLogCPM_A <- if (length(A_cols)>0) rowMeans(logcpm_norm[match(tt$gene_id,rownames(logcpm_norm)),A_cols,drop=FALSE]) else NA
  tt$meanCPM_B    <- if (length(B_cols)>0) rowMeans(cpm_norm[match(tt$gene_id,rownames(cpm_norm)),    B_cols,drop=FALSE]) else NA
  tt$meanLogCPM_B <- if (length(B_cols)>0) rowMeans(logcpm_norm[match(tt$gene_id,rownames(logcpm_norm)),B_cols,drop=FALSE]) else NA
  tt$group_A      <- AB$A; tt$group_B <- AB$B

  gene_idx          <- match(tt$gene_id, rownames(logcpm_all))
  tt$pct_detected_A <- if (length(A_cols)>0) rowMeans(logcpm_all[gene_idx,A_cols,drop=FALSE] > 0)*100 else 0
  tt$pct_detected_B <- if (length(B_cols)>0) rowMeans(logcpm_all[gene_idx,B_cols,drop=FALSE] > 0)*100 else 0
  tt$exclusive_to_A <- (tt$pct_detected_A > 50) & (tt$pct_detected_B == 0)
  tt$exclusive_to_B <- (tt$pct_detected_B > 50) & (tt$pct_detected_A == 0)

  tt$significant <- tt$FDR < opt$fdr_cutoff
  tt$sig_up      <- tt$significant & (tt$logFC > 0)
  tt$sig_down    <- tt$significant & (tt$logFC < 0)

  tt$priority_gene <- FALSE
  if (length(priority_genes) > 0) {
    prio_sym    <- resolve_to_symbols(priority_genes)
    sym_up      <- toupper(ifelse(is.na(tt$hgnc_symbol),"",tt$hgnc_symbol))
    sym_matches <- which(sym_up %in% prio_sym)
    if (length(sym_matches) > 0) {
      tt$priority_gene[sym_matches] <- TRUE
      prio_sig <- (tt$FDR[sym_matches] < opt$fdr_cutoff*2) & (abs(tt$logFC[sym_matches]) > 0.25)
      tt$significant[sym_matches[prio_sig]] <- TRUE
      log_msg("    Priority: %d found, %d sig (relaxed)", length(sym_matches), sum(prio_sig))
    }
  }

  tt <- annotate_genes(tt)

  n_sig <- sum(tt$significant, na.rm=TRUE)
  log_msg("    %d sig (FDR<%.2f): %d up, %d down | exclusive: %d to %s, %d to %s",
          n_sig, opt$fdr_cutoff,
          sum(tt$sig_up,na.rm=TRUE), sum(tt$sig_down,na.rm=TRUE),
          sum(tt$exclusive_to_A & tt$significant,na.rm=TRUE), AB$A,
          sum(tt$exclusive_to_B & tt$significant,na.rm=TRUE), AB$B)

  up_genes <- tt %>% filter(sig_up) %>% arrange(FDR, desc(abs(logFC))) %>%
    select(gene_id, hgnc_symbol, alias_symbols, symbol_search, gene_name,
           logFC, logCPM, PValue, FDR, meanCPM_A, meanCPM_B,
           pct_detected_A, pct_detected_B, exclusive_to_A, group_A, group_B)
  down_genes <- tt %>% filter(sig_down) %>% arrange(FDR, desc(abs(logFC))) %>%
    select(gene_id, hgnc_symbol, alias_symbols, symbol_search, gene_name,
           logFC, logCPM, PValue, FDR, meanCPM_A, meanCPM_B,
           pct_detected_A, pct_detected_B, exclusive_to_B, group_A, group_B)

  if (nrow(up_genes)   > 0) write.csv(up_genes,   file.path(outdir_run, paste0("UP_",   contrast_name,".csv")), row.names=FALSE)
  if (nrow(down_genes) > 0) write.csv(down_genes, file.path(outdir_run, paste0("DOWN_", contrast_name,".csv")), row.names=FALSE)

  res_all[[contrast_name]] <- tt
}

res_combined <- bind_rows(res_all)

# ---------- Integrity check + directed union ----------
suppressPackageStartupMessages({ library(readr); library(purrr) })
out_dir  <- outdir_run
all_path <- file.path(out_dir,"DE_all_contrasts.csv")

extract_contrast <- function(path)
  sub("\\.csv$","", sub("^.*?(UP_|DOWN_)","", basename(path)))

write.csv(res_combined, all_path, row.names=FALSE)

if (file.exists(all_path)) {
  all_df  <- read_csv(all_path, show_col_types=FALSE)
  idx_all <- all_df %>% distinct(gene_id, contrast) %>% mutate(in_all=TRUE)
  bad <- Filter(function(x) x$n > 0, lapply(
    c(Sys.glob(file.path(out_dir,"UP_*.csv")), Sys.glob(file.path(out_dir,"DOWN_*.csv"))),
    function(fp) {
      df <- suppressWarnings(read_csv(fp, show_col_types=FALSE))
      if (!nrow(df)) return(list(n=0))
      ct   <- extract_contrast(fp)
      miss <- df %>% transmute(gene_id, contrast=ct) %>% distinct() %>%
              left_join(idx_all, by=c("gene_id","contrast")) %>% filter(is.na(in_all))
      list(file=fp, n=nrow(miss))
    }
  ))
  if (length(bad)) log_msg("WARNING: Integrity check — entries missing from DE_all_contrasts")
  else             log_msg("Integrity check passed")
}

read_sig <- function(fp, dirn) {
  if (!file.exists(fp)) return(NULL)
  df <- suppressWarnings(read_csv(fp, show_col_types=FALSE))
  if (!nrow(df)) return(NULL)
  for (nm in setdiff(c("hgnc_symbol","alias_symbols","symbol_search","gene_name"), names(df)))
    df[[nm]] <- NA_character_
  df %>% mutate(direction=dirn, contrast=extract_contrast(fp)) %>%
    select(gene_id, hgnc_symbol, alias_symbols, symbol_search, gene_name,
           logFC, logCPM, PValue, FDR, meanCPM_A, meanCPM_B,
           pct_detected_A, pct_detected_B, group_A, group_B, contrast, direction)
}
sig_union <- bind_rows(
  lapply(Sys.glob(file.path(out_dir,"UP_*.csv")),   read_sig, dirn="up"),
  lapply(Sys.glob(file.path(out_dir,"DOWN_*.csv")), read_sig, dirn="down")
) %>% arrange(contrast, direction, FDR, desc(abs(logFC)))
if (nrow(sig_union)) write_csv(sig_union, file.path(out_dir,"DE_all_contrasts_SIGNIFICANT.csv"))

# ---------- Summary statistics ----------
summary_stats <- res_combined %>%
  group_by(contrast) %>%
  summarise(
    n_tested           = n(),
    n_sig_FDR05        = sum(FDR < 0.05,                   na.rm=TRUE),
    n_sig_FDR01        = sum(FDR < 0.01,                   na.rm=TRUE),
    n_up_FDR05         = sum(FDR < 0.05 & logFC > 0,       na.rm=TRUE),
    n_down_FDR05       = sum(FDR < 0.05 & logFC < 0,       na.rm=TRUE),
    n_sig_LFC1         = sum(FDR < 0.05 & abs(logFC) >= 1, na.rm=TRUE),
    n_sig_LFC2         = sum(FDR < 0.05 & abs(logFC) >= 2, na.rm=TRUE),
    mean_abs_logFC_sig = mean(abs(logFC[FDR<0.05]),         na.rm=TRUE),
    .groups="drop"
  )
write.csv(summary_stats, file.path(out_dir,"DE_summary.csv"), row.names=FALSE)

# ---------- Output matrices ----------
log_msg("Writing output files...")
lib_info <- y$samples %>% rownames_to_column("sample_id") %>%
  left_join(samples_dge, by="sample_id")
write.csv(lib_info,    file.path(out_dir,"library_info.csv"),     row.names=FALSE)
write.csv(cpm_norm,    file.path(out_dir,"CPM_normalized.csv"),    row.names=TRUE)
write.csv(logcpm_norm, file.path(out_dir,"logCPM_normalized.csv"), row.names=TRUE)

# ---------- QC plots ----------
log_msg("Generating QC report...")
cond_levels  <- levels(samples_dge$condition)
batch_levels <- levels(samples_dge$batch)
cond_cols    <- setNames(safe_discrete_cols(length(cond_levels)),         cond_levels)
batch_cols   <- setNames(safe_discrete_cols(length(batch_levels),"Set2"), batch_levels)

spike_num           <- parse_cell_count(samples_dge$condition)
spike_num[is.na(spike_num)] <- 0
spike_levels_sorted <- sort(unique(spike_num))
spike_labels_sorted <- format_spike_labels(spike_levels_sorted)
spike_label         <- factor(spike_labels_sorted[match(spike_num, spike_levels_sorted)],
                              levels=spike_labels_sorted)
spike_cols          <- setNames(safe_discrete_cols(length(spike_levels_sorted)), spike_labels_sorted)

pdf(file.path(out_dir,"QC_report.pdf"), width=16, height=12)

# Title page
plot.new()
text(0.5,0.7,"CTC Low-Input RNA-seq DGE Analysis", cex=2, font=2)
text(0.5,0.5, format(Sys.time(),"%Y-%m-%d %H:%M:%S"), cex=1)
text(0.5,0.3, sprintf("%d samples | %d conditions | %d genes",
                      nrow(samples_dge), length(cond_levels), nrow(y)), cex=1)

# Sample summary table
plot.new(); par(mar=c(0,0,2,0)); plot.window(xlim=c(0,1),ylim=c(0,1))
title("Sample Summary", font.main=2)
ss <- table(samples_dge$condition, samples_dge$batch)
text(0.5,0.8,"Samples per condition × batch:", font=2)
ty <- 0.65
for (i in 1:nrow(ss)) { text(0.5,ty,paste(rownames(ss)[i],":",paste(ss[i,],collapse=", ")),cex=0.9); ty <- ty-0.1 }

# Library sizes, norm factors, density, boxplot
par(mfrow=c(2,2), mar=c(5,4,4,2))
barplot(y$samples$lib.size/1e6, names.arg=rownames(y$samples),
        col=cond_cols[as.character(samples_dge$condition)],
        las=2, main="Library Sizes", ylab="Million reads", cex.names=0.7)
abline(h=median(y$samples$lib.size)/1e6, lty=2, col="red")
legend("topright", legend=cond_levels, fill=cond_cols, cex=0.8)

barplot(y$samples$norm.factors, names.arg=rownames(y$samples),
        col=cond_cols[as.character(samples_dge$condition)],
        las=2, main="TMM Normalisation Factors", ylab="Factor", cex.names=0.7)
abline(h=1, lty=2, col="red")

plot(density(logcpm_norm[,1]), lwd=2,
     col=cond_cols[as.character(samples_dge$condition[1])],
     main="Expression Density", xlab="log2 CPM", ylab="Density",
     ylim=c(0, max(apply(logcpm_norm,2,function(x) max(density(x)$y)))))
for (i in 2:ncol(logcpm_norm))
  lines(density(logcpm_norm[,i]), col=cond_cols[as.character(samples_dge$condition[i])], lwd=2)
legend("topright", legend=cond_levels, col=cond_cols, lty=1, lwd=2, cex=0.8)

boxplot(logcpm_norm, col=cond_cols[as.character(samples_dge$condition)],
        las=2, main="Expression Distribution", ylab="log2 CPM", cex.axis=0.7)

# PCA — uniform shape, coloured by spike level, labelled by sample_id
gene_var <- apply(logcpm_norm, 1, var, na.rm=TRUE)
keep_pca <- which(is.finite(gene_var) & gene_var > 0)
sel_idx  <- keep_pca[order(gene_var[keep_pca], decreasing=TRUE)][seq_len(min(2000,length(keep_pca)))]
X_pca    <- t(scale(t(logcpm_norm[sel_idx,,drop=FALSE]), center=TRUE, scale=TRUE))
X_pca[!is.finite(X_pca)] <- 0
pca      <- prcomp(t(X_pca), center=FALSE, scale.=FALSE)
ve       <- (pca$sdev^2) / sum(pca$sdev^2)
pc_lab   <- function(k) sprintf("PC%d (%.1f%%)", k, 100*ve[k])

pca_scores <- as.data.frame(pca$x[,1:min(5,ncol(pca$x)),drop=FALSE]) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(samples_dge %>% select(sample_id,condition,batch), by="sample_id") %>%
  mutate(spike=spike_num, spike_label=spike_label)
write.csv(pca_scores, file.path(out_dir,"PCA_scores.csv"), row.names=FALSE)

print(
  ggplot(pca_scores, aes(PC1,PC2,colour=spike_label)) +
    geom_hline(yintercept=0,linetype="dashed",linewidth=0.2) +
    geom_vline(xintercept=0,linetype="dashed",linewidth=0.2) +
    geom_point(shape=16, size=3.8, alpha=0.9) +
    ggrepel::geom_text_repel(aes(label=sample_id), min.segment.length=0.05,
                              box.padding=0.4, max.overlaps=Inf, seed=42, size=3) +
    scale_colour_manual(values=spike_cols, name="Spike level") +
    labs(title="PCA — PC1 vs PC2", x=pc_lab(1), y=pc_lab(2)) +
    theme_minimal(base_size=12) + theme(legend.position="right")
)
if (ncol(pca$x) >= 3)
  print(
    ggplot(pca_scores, aes(PC2,PC3,colour=spike_label)) +
      geom_hline(yintercept=0,linetype="dashed",linewidth=0.2) +
      geom_vline(xintercept=0,linetype="dashed",linewidth=0.2) +
      geom_point(shape=16, size=3.8, alpha=0.9) +
      ggrepel::geom_text_repel(aes(label=sample_id), min.segment.length=0.05,
                                box.padding=0.4, max.overlaps=Inf, seed=42, size=3) +
      scale_colour_manual(values=spike_cols, name="Spike level") +
      labs(title="PCA — PC2 vs PC3", x=pc_lab(2), y=pc_lab(3)) +
      theme_minimal(base_size=12) + theme(legend.position="right")
  )

# BCV & mean-variance
par(mfrow=c(1,2))
plotBCV(y, main="BCV", cex=0.5)
plotMeanVar(y, show.raw.vars=TRUE, show.tagwise.vars=TRUE,
            show.binned.common.disp.vars=FALSE, main="Mean-Variance")

# Per-contrast: MA, volcano, p-value hist, QQ, effect sizes, exclusive panel
par(mfrow=c(2,3))
for (contrast_name in names(res_all)) {
  tt  <- res_all[[contrast_name]]
  AB  <- parse_contrast_name(contrast_name)
  sig <- tt$FDR < opt$fdr_cutoff

  # MA
  col_ma <- ifelse(sig & tt$logFC>0,"red", ifelse(sig & tt$logFC<0,"blue","grey80"))
  plot(tt$logCPM, tt$logFC, pch=16, cex=0.3, col=col_ma,
       main=paste("MA:",contrast_name), xlab="Avg log CPM", ylab="log2 FC")
  abline(h=0,lty=2); abline(h=c(-opt$lfc_treat,opt$lfc_treat),lty=3,col="grey50")
  legend("topright", legend=c(paste0("Up (n=",sum(sig&tt$logFC>0),")"),
                               paste0("Down (n=",sum(sig&tt$logFC<0),")"), "NS"),
         col=c("red","blue","grey80"), pch=16, cex=0.8, bty="n")

  # Volcano
  yv <- -log10(tt$PValue)
  bg <- rep("grey85",nrow(tt)); bd <- rep(adjustcolor("grey40",0.6),nrow(tt)); lw <- rep(0.3,nrow(tt))
  bg[sig & tt$logFC>0] <- adjustcolor("tomato",0.6);    bd[sig & tt$logFC>0] <- "tomato4"
  bg[sig & tt$logFC<0] <- adjustcolor("royalblue",0.6); bd[sig & tt$logFC<0] <- "royalblue4"
  plot(tt$logFC, yv, pch=21, bg=bg, col=bd, lwd=lw, cex=0.6,
       main=paste("Volcano:",contrast_name), xlab="log2 FC", ylab="-log10(p)")
  abline(v=c(-opt$lfc_treat,opt$lfc_treat),lty=3,col="grey50")
  abline(h=-log10(0.05),lty=3,col="grey50")
  # Breast marker overlay — one row per canonical symbol
  hgnc_v   <- toupper(ifelse(is.na(tt$hgnc_symbol),"",tt$hgnc_symbol))
  cand_sym <- ifelse(hgnc_v %in% breast_marker_symbols, hgnc_v, NA_character_)
  eligible <- !is.na(cand_sym)
  if (any(eligible)) {
    idx      <- which(eligible)
    ord      <- order(cand_sym[idx], tt$FDR[idx], -abs(tt$logFC[idx]))
    idx_best <- idx[ord][!duplicated(cand_sym[idx[ord]])]
    points(tt$logFC[idx_best], yv[idx_best], pch=21, bg=NA, col="black", lwd=1.4, cex=1.0)
    text(tt$logFC[idx_best], yv[idx_best], labels=cand_sym[idx_best], pos=3, cex=0.7, offset=0.35)
  }
  legend("topright", legend=c(paste0("Up (n=",sum(sig&tt$logFC>0),")"),
                               paste0("Down (n=",sum(sig&tt$logFC<0),")"), "Breast markers"),
         pt.bg=c(adjustcolor("tomato",0.6),adjustcolor("royalblue",0.6),NA),
         pch=c(21,21,21), col=c("tomato4","royalblue4","black"),
         pt.lwd=c(0.3,0.3,1.4), bty="n", cex=0.8)

  # P-value histogram
  hist(tt$PValue, breaks=50, main=paste("P-value dist:",contrast_name),
       xlab="P-value", col="lightblue", border="darkblue")
  abline(v=0.05, lty=2, col="red", lwd=2)

  # QQ plot
  obs <- -log10(sort(tt$PValue)); exp_q <- -log10(ppoints(length(tt$PValue)))
  plot(exp_q, obs, pch=16, cex=0.3, main=paste("QQ:",contrast_name),
       xlab="Expected -log10(p)", ylab="Observed -log10(p)")
  abline(0,1,col="red",lty=2)

  # Effect size distribution
  hist(tt$logFC[sig], breaks=30,
       main=paste("Effect sizes (FDR<",opt$fdr_cutoff,"):",contrast_name),
       xlab="log2 FC", col="lightgreen", border="darkgreen")
  abline(v=c(-opt$lfc_treat,opt$lfc_treat),lty=2,col="red",lwd=2)

  # Exclusive genes panel
  if (sum(tt$exclusive_to_A | tt$exclusive_to_B, na.rm=TRUE) > 0) {
    bg_e <- rep("grey90",nrow(tt)); bd_e <- rep("grey50",nrow(tt)); lw_e <- rep(0.3,nrow(tt))
    a_ex <- tt$exclusive_to_A & sig; b_ex <- tt$exclusive_to_B & sig
    bg_e[a_ex] <- "deepskyblue"; bd_e[a_ex] <- "black"; lw_e[a_ex] <- 1.0
    bg_e[b_ex] <- "yellow";      bd_e[b_ex] <- "black"; lw_e[b_ex] <- 1.0
    plot(tt$logFC, yv, pch=21, bg=bg_e, col=bd_e, lwd=lw_e, cex=0.6,
         main=paste("Exclusive genes:",contrast_name), xlab="log2 FC", ylab="-log10(p)")
    cand_ex <- ifelse(hgnc_v %in% breast_marker_symbols, hgnc_v, NA_character_)
    elig_ex <- !is.na(cand_ex) & (tt$exclusive_to_A | tt$exclusive_to_B)
    if (any(elig_ex)) {
      idx  <- which(elig_ex); ord <- order(cand_ex[idx],-abs(tt$logFC[idx]),tt$FDR[idx])
      idx  <- idx[ord]; idx_best <- idx[!duplicated(cand_ex[idx])]
      points(tt$logFC[idx_best],yv[idx_best],pch=21,bg=NA,col="black",lwd=1.4,cex=1.0)
      text(tt$logFC[idx_best],yv[idx_best],labels=cand_ex[idx_best],pos=3,cex=0.7,offset=0.35)
    }
    legend("topright",
           legend=c(paste("Exclusive to",AB$A),paste("Exclusive to",AB$B),"Breast markers"),
           pt.bg=c("deepskyblue","yellow",NA), col=c("black","black","black"),
           pch=c(21,21,21), pt.lwd=c(1,1,1.4), cex=0.8, bty="n")
  }
}

# Sample correlation heatmap
if (ncol(logcpm_norm) <= 50) {
  so       <- order(samples_dge$condition, samples_dge$sample_id)
  lc_ord   <- logcpm_norm[,so]; sa_ord <- samples_dge[so,]
  cor_mat  <- cor(lc_ord, method="spearman")
  slabs    <- create_sample_labels(sa_ord, include_batch=TRUE)
  colnames(cor_mat) <- rownames(cor_mat) <- slabs
  heatmap(cor_mat, col=colorRampPalette(c("blue","white","red"))(100),
          margins=c(12,12), main="Sample Correlation (Spearman)", cexRow=0.6, cexCol=0.6)
}

# Summary page
plot.new(); par(mar=c(0,0,2,0)); plot.window(xlim=c(0,1),ylim=c(0,1))
title("Analysis Summary", font.main=2)
total_de <- sum(summary_stats$n_sig_FDR05)
text(0.5,0.85,sprintf("Total DE genes (FDR < %.2f): %d",opt$fdr_cutoff,total_de),cex=1.2,font=2)
ty2 <- 0.70
for (i in 1:min(10,nrow(summary_stats))) {
  text(0.5,ty2,sprintf("%s: %d DE (%d up / %d down)", summary_stats$contrast[i],
                        summary_stats$n_sig_FDR05[i], summary_stats$n_up_FDR05[i],
                        summary_stats$n_down_FDR05[i]),cex=0.85)
  ty2 <- ty2-0.05
}

# ---- Breast marker linearity heatmaps (spike-level scaling assessment) ----
genes_mat_clean <- sub("\\.\\d+$","",rownames(logcpm_norm))
ens_from_sym    <- suppressMessages(AnnotationDbi::mapIds(org_db, keys=breast_marker_symbols,
                    keytype="SYMBOL", column="ENSEMBL", multiVals="first"))
ens_from_sym    <- sub("\\.\\d+$","",unname(ens_from_sym))
marker_rows     <- which(genes_mat_clean %in% ens_from_sym)

if (length(marker_rows) > 0) {
  X_hm    <- logcpm_norm[marker_rows,,drop=FALSE]
  ord_hm  <- order(spike_num, samples_dge$sample_id)
  X_hm    <- X_hm[,ord_hm,drop=FALSE]
  samp_hm <- samples_dge[ord_hm,]; spike_hm <- spike_num[ord_hm]
  colnames(X_hm) <- paste0(samp_hm$sample_id,"\n",format_spike_labels(spike_hm))
  symb_hm <- suppressMessages(AnnotationDbi::mapIds(org_db,keys=genes_mat_clean[marker_rows],
               keytype="ENSEMBL",column="SYMBOL",multiVals="first"))
  rownames(X_hm) <- ifelse(is.na(symb_hm)|symb_hm=="",rownames(X_hm),symb_hm)

  # Delta-based row ordering (Up > Down > Flat)
  df_long_hm <- as.data.frame(X_hm) %>% tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene,names_to="sample",values_to="log2cpm") %>%
    mutate(spike=spike_hm[match(sample,colnames(X_hm))]) %>% filter(is.finite(spike))
  gm_hm  <- df_long_hm %>% group_by(gene,spike) %>% summarise(mean=mean(log2cpm),.groups="drop")
  levs_hm <- sort(unique(gm_hm$spike))
  delta_hm <- gm_hm %>%
    filter(spike %in% c(levs_hm[1],levs_hm[length(levs_hm)])) %>%
    mutate(wh=ifelse(spike==levs_hm[length(levs_hm)],"high","low")) %>%
    select(gene,wh,mean) %>% pivot_wider(names_from=wh,values_from=mean) %>%
    mutate(delta=high-low,
           Regulation=case_when(delta>=0.2~"Up",delta<=-0.2~"Down",TRUE~"Flat"))
  row_ord <- delta_hm %>%
    mutate(ok=case_when(Regulation=="Up"~1L,Regulation=="Down"~2L,TRUE~3L)) %>%
    arrange(ok,desc(delta),gene) %>% pull(gene)
  row_ord <- row_ord[row_ord %in% rownames(X_hm)]
  X_hm    <- X_hm[row_ord,,drop=FALSE]

  ann_col_hm    <- data.frame(Condition=factor(samp_hm$condition))
  rownames(ann_col_hm) <- colnames(X_hm)
  ann_colors_hm <- list(Condition=setNames(safe_discrete_cols(nlevels(ann_col_hm$Condition),"Set2"),
                                           levels(ann_col_hm$Condition)))
  rle_c   <- rle(as.character(ann_col_hm$Condition))
  gaps_hm <- cumsum(rle_c$lengths); gaps_hm <- gaps_hm[-length(gaps_hm)]

  # Heatmap A: absolute log2CPM ordered by spike level
  lims_abs   <- quantile(as.numeric(X_hm), probs=c(0.02,0.98), na.rm=TRUE)
  breaks_abs <- seq(lims_abs[1], lims_abs[2], length.out=101)
  pheatmap(X_hm, color=colorRampPalette(brewer.pal(9,"YlOrRd"))(100),
           breaks=breaks_abs, scale="none",
           cluster_rows=FALSE, cluster_cols=FALSE,
           annotation_col=ann_col_hm, annotation_colors=ann_colors_hm, gaps_col=gaps_hm,
           show_colnames=TRUE, show_rownames=TRUE, border_color=NA,
           main=sprintf("Heatmap A — Breast markers (absolute log2CPM)\nOrdered by spike level; rows: Up → Down → Flat"),
           fontsize_row=10, fontsize_col=8, angle_col=90, cellheight=22)

  # Heatmap B: logFC vs lowest spike — reveals linearity of response
  gm_wide <- gm_hm %>% pivot_wider(id_cols=gene,names_from=spike,values_from=mean) %>% as.data.frame()
  rownames(gm_wide) <- gm_wide$gene; gm_wide$gene <- NULL
  cs      <- as.numeric(colnames(gm_wide)); oc <- order(cs)
  gm_wide <- as.matrix(gm_wide[,oc,drop=FALSE]); cs <- cs[oc]
  gm_wide <- sweep(gm_wide,1,gm_wide[,1],"-")
  gm_wide <- gm_wide[rownames(X_hm),,drop=FALSE]
  colnames(gm_wide) <- format_spike_labels(cs)
  mx_lfc  <- max(abs(gm_wide), na.rm=TRUE); if (mx_lfc < 0.01) mx_lfc <- 0.1
  pheatmap(gm_wide, color=colorRampPalette(c("blue","white","red"))(100),
           breaks=seq(-mx_lfc,mx_lfc,length.out=101), scale="none",
           cluster_rows=TRUE, cluster_cols=FALSE,
           show_colnames=TRUE, show_rownames=TRUE, border_color=NA,
           main=sprintf("Heatmap B — Breast markers LogFC vs lowest spike (%.0f cells)\nBlue = down, Red = up",cs[1]),
           fontsize_row=10, fontsize_col=10, angle_col=90, cellheight=22)
}

# ---- Top-at-max-spike heatmaps (global gene-level linearity) ----
plot_top_at_max_spike <- function(logcpm_norm, samples, org_db,
                                   direction="up", top_n=30, value="logfc",
                                   main_prefix="Top-at-max spike") {
  spike_n <- parse_cell_count(samples$condition); spike_n[is.na(spike_n)] <- 0
  df_long <- as.data.frame(logcpm_norm) %>% tibble::rownames_to_column("ens") %>%
    pivot_longer(-ens,names_to="sample",values_to="log2cpm") %>%
    mutate(spike=spike_n[match(sample,colnames(logcpm_norm))]) %>% filter(is.finite(spike))
  gm <- df_long %>% group_by(ens,spike) %>% summarise(mean=mean(log2cpm),.groups="drop")
  levs <- sort(unique(gm$spike))
  if (length(levs) < 2) { plot.new(); title(paste(main_prefix,": need ≥2 spike levels")); return(invisible(NULL)) }
  deltas <- full_join(gm %>% filter(spike==levs[1])             %>% select(ens,mean_low=mean),
                      gm %>% filter(spike==levs[length(levs)]) %>% select(ens,mean_high=mean),
                      by="ens") %>%
    mutate(delta=mean_high-mean_low) %>% filter(is.finite(delta))
  top_genes <- if (direction=="up") deltas %>% arrange(desc(delta)) %>% slice_head(n=top_n) %>% pull(ens)
               else                 deltas %>% arrange(delta)       %>% slice_head(n=top_n) %>% pull(ens)
  if (!length(top_genes)) { plot.new(); title(paste(main_prefix,"-",direction,": no genes")); return(invisible(NULL)) }
  mat <- gm %>% filter(ens %in% top_genes) %>%
    pivot_wider(id_cols=ens,names_from=spike,values_from=mean) %>% as.data.frame()
  rownames(mat) <- mat$ens; mat$ens <- NULL
  cs  <- as.numeric(colnames(mat)); oc <- order(cs)
  mat <- as.matrix(mat[,oc,drop=FALSE]); cs <- cs[oc]
  rownames(mat) <- .ens_to_symbol(rownames(mat), org_db)
  if (value == "logfc") {
    mat <- sweep(mat,1,mat[,1],"-")
    mabs <- max(abs(mat),na.rm=TRUE); if (mabs<0.01) mabs<-0.1
    brks <- seq(-mabs,mabs,length.out=101); cols <- colorRampPalette(c("blue","white","red"))(100)
    scale_title <- sprintf("LogFC vs lowest spike (%.0f)",cs[1])
  } else {
    rng  <- quantile(as.numeric(mat),probs=c(0.02,0.98),na.rm=TRUE)
    brks <- seq(rng[1],rng[2],length.out=101); cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
    scale_title <- "Group-mean log2CPM"
  }
  colnames(mat) <- format_spike_labels(cs)
  pheatmap(mat, color=cols, breaks=brks, scale="none",
           cluster_rows=TRUE, cluster_cols=FALSE,
           show_colnames=TRUE, show_rownames=TRUE, border_color=NA,
           main=sprintf("%s (%s)\n%s | Max=%.0f cells | Min=%.0f cells",
                        main_prefix, ifelse(direction=="up","Upregulated","Downregulated"),
                        scale_title, max(cs), min(cs)),
           fontsize_row=9, fontsize_col=10, angle_col=90, cellwidth=100, cellheight=20)
}

plot_top_at_max_spike(logcpm_norm, samples_dge, org_db, direction="up",   top_n=30, value="logfc")
plot_top_at_max_spike(logcpm_norm, samples_dge, org_db, direction="down", top_n=30, value="logfc")

dev.off()
log_msg("QC report saved: QC_report.pdf")

# ---------- Gene sets for enrichment ----------
log_msg("Writing gene sets...")
gene_sets_dir <- file.path(out_dir,"gene_sets")
dir.create(gene_sets_dir, showWarnings=FALSE)
for (contrast_name in names(res_all)) {
  tt <- res_all[[contrast_name]]
  gene_sets <- list(
    up_fdr05   = tt %>% filter(FDR<0.05,logFC>0)  %>% pull(hgnc_symbol) %>% na.omit(),
    down_fdr05 = tt %>% filter(FDR<0.05,logFC<0)  %>% pull(hgnc_symbol) %>% na.omit(),
    up_fdr01   = tt %>% filter(FDR<0.01,logFC>0)  %>% pull(hgnc_symbol) %>% na.omit(),
    down_fdr01 = tt %>% filter(FDR<0.01,logFC<0)  %>% pull(hgnc_symbol) %>% na.omit(),
    up_lfc2    = tt %>% filter(FDR<0.05,logFC>2)  %>% pull(hgnc_symbol) %>% na.omit(),
    down_lfc2  = tt %>% filter(FDR<0.05,logFC< -2)%>% pull(hgnc_symbol) %>% na.omit()
  )
  for (sn in names(gene_sets))
    if (length(gene_sets[[sn]]) > 0)
      writeLines(as.character(gene_sets[[sn]]),
                 file.path(gene_sets_dir, paste0(contrast_name,"_",sn,".txt")))
  ranked <- tt %>% filter(!is.na(hgnc_symbol)) %>%
    mutate(rank_score=-log10(PValue)*sign(logFC)) %>%
    arrange(desc(rank_score)) %>% select(hgnc_symbol,rank_score)
  write.table(ranked, file.path(gene_sets_dir,paste0(contrast_name,"_ranked.rnk")),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# ---------- Session info + README ----------
writeLines(capture.output(sessionInfo()), file.path(out_dir,"session_info.txt"))

writeLines(sprintf("
CTC Low-Input RNA-seq DGE Analysis Results
Generated : %s

INPUT:
  Samples   : %s
  GTF       : %s
  Quant dir : %s

PARAMETERS:
  FDR cutoff  : %.2f
  LFC thresh  : %.1f
  Test type   : %s
  Min CPM     : %.2f
  Low-input   : %s
  Species     : %s

SAMPLES:
  Total          : %d
  Blanks         : %d
  Positive ctrls : %d
  DGE samples    : %d
  Conditions     : %s
  Batches        : %s

ANALYSIS:
  Genes tested : %d
  Contrasts    : %d
  Total DE     : %d (FDR < %.2f)

KEY OUTPUTS:
  DE_all_contrasts.csv              All results (every gene, every contrast)
  DE_all_contrasts_SIGNIFICANT.csv  Significant results only
  DE_summary.csv                    Summary stats per contrast
  UP_*.csv / DOWN_*.csv             Per-contrast gene lists
  CPM_normalized.csv                Normalised CPM matrix
  logCPM_normalized.csv             log2 normalised CPM matrix
  QC_report.pdf                     Full QC + heatmap report
  PCA_scores.csv                    PCA coordinates
  gene_sets/                        Gene lists + ranked files for GSEA
  library_info.csv                  Sample + library metadata
  design_matrix.csv                 edgeR design matrix
  blank_contamination_report.csv    Genes detected in blank captures
  positive_control_report.csv       Positive control pass/fail
  session_info.txt                  R package versions
  pipeline.log                      Full run log
",
  format(Sys.time(),"%Y-%m-%d %H:%M:%S"),
  opt$samples, opt$gtf, opt$quant_dir,
  opt$fdr_cutoff, opt$lfc_treat,
  ifelse(opt$use_treat,"glmTreat","glmQLFTest"),
  opt$min_cpm, ifelse(opt$low_input_mode,"YES","NO"), opt$species,
  nrow(samples), sum(blank_idx), sum(posctrl_idx), sum(dge_idx),
  paste(levels(samples$condition),collapse=", "),
  paste(levels(samples$batch),collapse=", "),
  nrow(y), length(contrast_list), sum(summary_stats$n_sig_FDR05), opt$fdr_cutoff
), file.path(out_dir,"README.txt"))

# ---------- Completion ----------
total_sig <- sum(summary_stats$n_sig_FDR05)
log_msg("=== Pipeline Complete ===")
log_msg("Total significant DE genes: %d", total_sig)
log_msg("Results: %s", out_dir)

if (total_sig == 0) {
  log_msg("WARNING: No significant DE genes.")
  log_msg("For low-input CTC data, consider:")
  log_msg("  --low_input_mode")
  log_msg("  --fdr_cutoff 0.2")
  log_msg("  --min_cpm 0.05")
  log_msg("  --priority_genes 'ERBB2,ESR1,PGR,MKI67,VIM,CDH1'")
}
