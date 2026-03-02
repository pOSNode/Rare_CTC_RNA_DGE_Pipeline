# Methods

## Read QC and Trimming

Paired-end reads were quality-trimmed using fastp (v0.23.4) with automatic adapter detection (`--detect_adapter_for_pe`), mismatch correction (`--correction`), and overrepresentation analysis enabled. In low-input mode, a minimum read length of 20 bp is enforced and low-complexity reads are removed to reduce noise from degraded RNA fragments.

## Quantification

Transcript-level quantification was performed using Salmon (v1.10.3) in quasi-mapping mode against the GRCh38 reference genome (Gencode v46). GC bias, sequence bias, and positional bias corrections were applied (`--gcBias --seqBias --posBias`). For low-input samples, 50 bootstrap replicates were used to characterise quantification uncertainty. Minimum score fraction was set to 0.65 to permit alignment of partially degraded reads.

## Multi-run Merging

For samples sequenced across multiple runs, per-run Salmon quantification directories were merged at the quantification level using tximport (v1.30.0). Length-scaled TPM values (`countsFromAbundance = "lengthScaledTPM"`) were used to correct for transcript length differences between runs.

## Filtering

Gene-level filtering was performed using edgeR's `filterByExpr` with a minimum CPM threshold of 0.1 (reduced from the standard 0.5) and minimum sample threshold of 2. In addition, genes showing condition-exclusive expression (CPM ≥ 0.05 in at least 1 sample of one condition and absent from all other conditions) were retained irrespective of the base filter. This step preserves epithelial and mesenchymal marker genes that are expressed specifically in CTCs and absent in contaminating WBCs.

## Normalisation

Library size normalisation was performed using TMMwsp (Trimmed Mean of M-values with singleton pairing), which is the recommended method for sparse count matrices with many zero counts. Standard TMM normalisation is unstable in this context because it requires shared non-zero counts across samples for reference selection. Normalisation factors outside a 5-fold range trigger an automatic switch to upper-quartile normalisation with a logged warning.

## Differential Expression

Quasi-likelihood F-tests were performed using edgeR's `glmQLFit` and `glmQLFTest` with robust dispersion estimation (`robust = TRUE`). Where batch information was available, batch was included as a covariate in the design matrix. All pairwise condition contrasts were tested. A log-CPM transformation with prior count of 0.5 was used for visualisation, reducing shrinkage of low-count genes relative to the standard prior count of 2. Significance was defined as FDR < 0.1 (Benjamini-Hochberg correction). Priority genes specified by the user were tested with relaxed criteria (FDR < 0.2, |logFC| > 0.25).

## Audit and Reproducibility

Every pipeline step writes a timestamped entry to a central audit log. At pipeline completion, a provenance JSON file records pipeline version, all runtime parameters, input file MD5 checksums, and tool versions. Nextflow's built-in tracing records per-task CPU, memory, and runtime. Together these outputs satisfy the reproducibility requirements of a GxP-adjacent research environment.
