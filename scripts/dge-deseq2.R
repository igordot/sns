#!/usr/bin/env Rscript


##
## Differential gene expression with DESeq2.
##
## usage: Rscript --vanilla dge-deseq2.R genome_build genes.gtf counts.txt groups.csv
##


# increase output width
options(width = 120)
# print warnings as they occur
options(warn = 1)
# java heap size
options(java.parameters = "-Xmx8G")

# get scripts directory (directory of this file) and load relevant functions
args_all = commandArgs(trailingOnly = FALSE)
scripts_dir = normalizePath(dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(paste0(scripts_dir, "/load-install-packages.R"))
source(paste0(scripts_dir, "/deseq2-pca.R"))
source(paste0(scripts_dir, "/deseq2-compare.R"))
source(paste0(scripts_dir, "/plot-volcano.R"))
source(paste0(scripts_dir, "/plot-heatmap.R"))
source(paste0(scripts_dir, "/gse-fgsea.R"))

# relevent arguments
args = commandArgs(trailingOnly = TRUE)
genome_build = args[1]
genes_gtf = args[2]
counts_table_file = args[3]
groups_table_file = args[4]

# check for arguments
if (length(args) < 4) stop("not enough arguments provided")

# check that input files exist
if (!file.exists(counts_table_file)) stop("file does not exist: ", counts_table_file)
if (!file.exists(groups_table_file)) stop("file does not exist: ", groups_table_file)

# create separate directories for certain output files
r_dir = "r-data"
if (!dir.exists(r_dir)) dir.create(r_dir)

# for general data manipulation
load_install_packages("magrittr")
load_install_packages("tibble")
load_install_packages("dplyr")
load_install_packages("tidyr")
load_install_packages("readr")
load_install_packages("glue")
# for differenial expression
load_install_packages("DESeq2")
load_install_packages("ashr")
load_install_packages("genefilter")
# for processing GTF (for gene lengths for FPKMs)
load_install_packages("rtracklayer")
# for exporting Excel xlsx files
load_install_packages("writexl")
# for color scheme
load_install_packages("RColorBrewer")
# for standard plotting
load_install_packages("ggplot2")
load_install_packages("ggrepel")
load_install_packages("cowplot")
# for heatmaps
load_install_packages("pheatmap")
# for gene set enrichment (pathways)
load_install_packages("msigdbr")
load_install_packages("fgsea")

message(" ========== import inputs ========== ")

# import counts table
counts_table = read.delim(file = counts_table_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
message("input counts table gene num:      ", nrow(counts_table))
message("input counts table sample num:    ", ncol(counts_table))
message("input counts table sample names:  ", toString(colnames(counts_table)))
message("")

# import groups table
groups_table = read.csv(file = groups_table_file, header = TRUE, row.names = 1, colClasses = "factor")
message("sample groups table sample num:   ", nrow(groups_table))
message("sample groups table sample names: ", toString(rownames(groups_table)))
message("sample groups table group names:  ", toString(colnames(groups_table)))
message("")

# check that all samples from the groups table are found in the counts table
diff_samples = setdiff(rownames(groups_table), colnames(counts_table))
if (length(diff_samples)) stop("some samples not in counts table: ", toString(diff_samples))

# subset to samples in groups table (also sets samples to be in the same order)
counts_table = counts_table[, rownames(groups_table)]
message("subset counts table gene num:     ", nrow(counts_table))
message("subset counts table sample num:   ", ncol(counts_table))
message("")

# group info (use the first column for grouped comparisons)
group_name = colnames(groups_table)[1]
message("group name: ", group_name)
# not using levels() to preserve order
group_levels = groups_table[, group_name] %>% as.character() %>% unique()
message("group levels: ", toString(group_levels))
message("")

# design formula
design_formula = formula(glue("~ {group_name}"))
if (length(group_levels) == 1) { design_formula = formula("~ 1") }
message("design formula: ", design_formula)

message(" ========== import GTF genes annotations ========== ")

# import genes GTF file
genes_granges = import(genes_gtf)
message("GTF total entries:      ", length(genes_granges))

# extract gene lengths (sum of exons)
exons_granges = genes_granges[genes_granges$type == "exon"]
exons_by_gene = split(exons_granges, exons_granges$gene_name)
message("GTF num genes:          ", length(exons_by_gene))
gene_lengths = exons_by_gene %>% reduce() %>% width() %>% sum()
message("GTF mean gene length:   ", round(mean(gene_lengths), 2))
message("GTF median gene length: ", median(gene_lengths))
message("")

# save gene annotations
genes_tbl = as_tibble(genes_granges)
if ("gene" %in% genes_tbl$type) {
  genes_tbl = dplyr::filter(genes_tbl, type == "gene")
  genes_tbl = dplyr::rename(genes_tbl, chr = seqnames)
  genes_tbl = dplyr::select(genes_tbl, one_of(c("gene_name", "gene_id", "chr", "start", "end", "strand", "gene_type")))
  genes_tbl = dplyr::select(genes_tbl, gene_name, gene_id, everything())
  genes_tbl = dplyr::distinct(genes_tbl)
  genes_tbl = dplyr::arrange(genes_tbl, gene_name, gene_id)
  write_csv(genes_tbl, "genes.csv")
}

message(" ========== normalize ========== ")

# import raw counts and create DESeq object
# since v1.16 (11/2016), betaPrior is set to FALSE and shrunken LFCs are obtained afterwards using lfcShrink
dds = DESeqDataSetFromMatrix(countData = counts_table, colData = groups_table, design = design_formula)
dds = DESeq(dds, parallel = FALSE)

# add gene lengths (used to generate FPKM values)
if (identical(sort(names(gene_lengths)), sort(rownames(dds)))) {
  mcols(dds)$basepairs = gene_lengths[rownames(dds)]
} else {
  message("GTF num genes:          ", length(gene_lengths))
  message("counts table num genes: ", nrow(counts_table))
  message("dds genes:              ", nrow(dds))
  stop("genes in the GTF and the DESeq dataset object do not match")
}

# VST
vsd = varianceStabilizingTransformation(dds, blind = TRUE)

# save DESeqDataSet and VST DESeqTransform objects
saveRDS(dds, file = glue("{r_dir}/deseq2.dds.rds"))
saveRDS(vsd, file = glue("{r_dir}/deseq2.vsd.rds"))
Sys.sleep(1)

# export counts
raw_counts_table = counts(dds, normalized = FALSE) %>% as_tibble(rownames = "gene")
write_csv(raw_counts_table, "counts.raw.csv.gz")
norm_counts_table = counts(dds, normalized = TRUE) %>% round(3) %>% as_tibble(rownames = "gene")
write_csv(norm_counts_table, "counts.normalized.csv.gz")
write_xlsx(list(normalized_counts = norm_counts_table), "counts.normalized.xlsx")
Sys.sleep(1)

# export average counts per group
if (length(group_levels) > 1) {
  norm_counts_mat = counts(dds, normalized = TRUE)
  norm_counts_table = sapply(group_levels, function(x) rowMeans(norm_counts_mat[, colData(dds)[, group_name] == x]))
  norm_counts_table = norm_counts_table %>% round(3) %>% as_tibble(rownames = "gene")
  write_csv(norm_counts_table, glue("counts.normalized.{group_name}.csv"))
  Sys.sleep(1)
}

# export FPMs/CPMs (fragments/counts per million mapped fragments)
# robust version uses size factors to normalize rather than taking the column sums of the raw counts
# not using the robust median ratio method to generate the classic values (comparable across experiments)
cpm_matrix = fpm(dds, robust = FALSE)
cpm_table = cpm_matrix %>% round(3) %>% as_tibble(rownames = "gene")
write_csv(cpm_table, "counts.cpm.csv.gz")
write_xlsx(list(CPMs = cpm_table), "counts.cpm.xlsx")
Sys.sleep(1)

# export FPKMs (fragment counts normalized per kilobase of feature length per million mapped fragments)
fpkm_matrix = fpkm(dds, robust = FALSE)
fpkm_table = fpkm_matrix %>% round(3) %>% as_tibble(rownames = "gene")
write_csv(fpkm_table, "counts.fpkm.csv.gz")
write_xlsx(list(FPKMs = fpkm_table), "counts.fpkm.xlsx")
Sys.sleep(1)

# export TPMs (transcripts per million)
tpm_matrix = apply(fpkm_matrix, 2, function(x) { exp(log(x) - log(sum(x)) + log(1e6)) })
tpm_table = tpm_matrix %>% round(3) %>% as_tibble(rownames = "gene")
write_csv(tpm_table, "counts.tpm.csv.gz")
write_xlsx(list(TPMs = tpm_table), "counts.tpm.xlsx")
Sys.sleep(1)

# export variance stabilized counts
vsd_table = assay(vsd) %>% round(3) %>% as_tibble(rownames = "gene")
write_csv(vsd_table, "counts.vst.csv.gz")
Sys.sleep(1)

message(" ========== QC ========== ")

# sparsity plot
png("plot.sparsity.png", width = 6, height = 6, units = "in", res = 300)
print(plotSparsity(dds, normalized = TRUE))
dev.off()
Sys.sleep(1)

# PCA plot
pca_plot = deseq2_pca(vsd, intgroup = group_name, ntop = 1000, point_labels = TRUE)
save_plot("plot.pca.png", pca_plot, base_height = 6, base_width = 8, units = "in")
Sys.sleep(1)
save_plot("plot.pca.pdf", pca_plot, base_height = 6, base_width = 8, units = "in")
Sys.sleep(1)

# PCA plot without labels for larger projects
if (ncol(dds) > 10) {
  pca_plot = deseq2_pca(vsd, intgroup = group_name, ntop = 1000, point_labels = FALSE)
  save_plot("plot.pca.nolabels.png", pca_plot, base_height = 6, base_width = 8, units = "in")
  Sys.sleep(1)
  save_plot("plot.pca.nolabels.pdf", pca_plot, base_height = 6, base_width = 8, units = "in")
  Sys.sleep(1)
}

message(" ========== differential expression ========== ")

# perform comparisons for all combinations of group levels
if (length(group_levels) > 1) {
  group_levels_combinations = combn(group_levels, m = 2, simplify = TRUE)
  for (combination_num in 1:ncol(group_levels_combinations)) {
    # numerator is second in order (order should match the input table group order)
    level_numerator = group_levels_combinations[2, combination_num]
    level_denominator = group_levels_combinations[1, combination_num]
    message(glue("comparison : {group_name} : {level_numerator} vs {level_denominator}"))
    deseq2_compare(deseq_dataset = dds, contrast = c(group_name, level_numerator, level_denominator), genome = genome_build)
  }
}

# delete Rplots.pdf (left by some plotting functions)
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
