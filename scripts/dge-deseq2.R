#!/usr/bin/env Rscript


##
## Differential gene expression with DESeq2.
##
## usage: Rscript --vanilla dge-deseq2.R genes.gtf counts_table.txt groups_table.csv
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
source(paste0(scripts_dir, "/deseq2-results.R"))
source(paste0(scripts_dir, "/deseq2-heatmap.R"))

# relevent arguments
args = commandArgs(trailingOnly = TRUE)
genes_gtf = args[1]
counts_table_file = args[2]
groups_table_file = args[3]

# check for arguments
if (length(args) < 3) stop("not enough arguments provided")

# check that input files exist
if (!file.exists(counts_table_file)) stop("file does not exist: ", counts_table_file)
if (!file.exists(groups_table_file)) stop("file does not exist: ", groups_table_file)

# create separate directories for certain output files
r_dir = "r-data"
if (!dir.exists(r_dir)) dir.create(r_dir)

# load relevant packages
load_install_packages("magrittr")
load_install_packages("DESeq2")
load_install_packages("genefilter")
load_install_packages("rtracklayer")
load_install_packages("RColorBrewer")
load_install_packages("ggplot2")
load_install_packages("ggrepel")
load_install_packages("cowplot")
load_install_packages("xlsx")
load_install_packages("pheatmap")

message(" ========== import inputs ========== ")

# import genes GTF file
genes_granges = import(genes_gtf)
message("GTF total entries: ", length(genes_granges))

# extract gene lengths (sum of exons)
exons_granges = genes_granges[genes_granges$type == "exon"]
exons_by_gene = split(exons_granges, exons_granges$gene_name)
message("GTF genes: ", length(exons_by_gene))
gene_lengths = exons_by_gene %>% reduce %>% width %>% sum
message("GTF mean gene length: ", round(mean(gene_lengths), 1))
message("GTF median gene length: ", median(gene_lengths))
message("")

# import counts table
counts_table = read.delim(file = counts_table_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
message("input counts table gene num: ", nrow(counts_table))
message("input counts table sample num: ", ncol(counts_table))
message("input counts table samples: ", toString(colnames(counts_table)))
message("")

# import groups table
groups_table = read.csv(file = groups_table_file, header = TRUE, row.names = 1, colClasses = "factor")
message("sample groups table sample num: ", nrow(groups_table))
message("sample groups table samples: ", toString(rownames(groups_table)))
message("sample groups table groups: ", toString(colnames(groups_table)))
message("")

# check that all samples from the groups table are found in the counts table
diff_samples = setdiff(rownames(groups_table), colnames(counts_table))
if (length(diff_samples)) stop("some samples not in counts table: ", toString(diff_samples))

# subset to samples in groups table (also sets samples to be in the same order)
counts_table = counts_table[, rownames(groups_table)]
message("subset counts table gene num: ", nrow(counts_table))
message("subset counts table sample num: ", ncol(counts_table))
message("")

# group info (use the first column for grouped comparisons)
group_name = colnames(groups_table)[1]
message("group name: ", group_name)
group_levels = levels(groups_table[, group_name])
message("group levels: ", toString(group_levels))
message("")

# design formula
design_formula = formula(paste("~", group_name))
message("design formula: ", design_formula)

message(" ========== normalize ========== ")

# import raw counts and create DESeq object
# since v1.16 (11/2016), betaPrior is set to FALSE and shrunken LFCs are obtained afterwards using lfcShrink
dds = DESeqDataSetFromMatrix(countData = counts_table, colData = groups_table, design = design_formula)
dds = DESeq(dds, betaPrior = TRUE, parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 4))

# add gene lengths (used to generate FPKM values)
if (identical(sort(names(gene_lengths)), sort(rownames(dds)))) {
  mcols(dds)$basepairs = gene_lengths[rownames(dds)]
} else {
  message("GTF num genes: ", length(gene_lengths))
  message("counts table num genes: ", nrow(counts_table))
  message("dds genes: ", nrow(dds))
  stop("genes in the GTF and the DESeq dataset object do not match")
}

# VST
vsd = varianceStabilizingTransformation(dds, blind = TRUE)

# save DESeqDataSet and VST DESeqTransform objects
saveRDS(dds, file = paste0(r_dir, "/deseq2.dds.rds"))
saveRDS(vsd, file = paste0(r_dir, "/deseq2.vsd.rds"))

# export counts
write.csv(counts(dds, normalized = FALSE), file = "counts.raw.csv")
norm_counts_table = counts(dds, normalized = TRUE) %>% round(2)
write.csv(norm_counts_table, file = "counts.normalized.csv")
write.xlsx2(x = norm_counts_table, file = "counts.normalized.xlsx", sheetName = "normalized counts")

# export FPKMs (fragment counts normalized per kilobase of feature length per million mapped fragments)
fpkm_table = fpkm(dds, robust = TRUE) %>% round(2)
write.csv(fpkm_table, file = "counts.fpkm.csv")
write.xlsx2(x = fpkm_table, file = "counts.fpkm.xlsx", sheetName = "FPKMs")

# export variance stabilized counts
write.csv(round(assay(vsd), digits = 2), file = "counts.vst.csv")

message(" ========== QC ========== ")

# sparsity plot
png("plot.sparsity.png", width = 6, height = 6, units = "in", res = 300)
print(plotSparsity(dds, normalized = TRUE))
dev.off()

# PCA plot
pca_plot = deseq2_pca(vsd, intgroup = group_name, ntop = 1000)
save_plot("plot.pca.pdf", pca_plot, base_width = 8, base_height = 5, units = "in", dpi = 300)
save_plot("plot.pca.png", pca_plot, base_width = 8, base_height = 5, units = "in", dpi = 300)

message(" ========== differential expression ========== ")

# perform comparisons for all combinations of group levels
group_levels_combinations = combn(group_levels, m = 2, simplify = TRUE)
for (combination_num in 1:ncol(group_levels_combinations)) {
  # numerator is second in order (usually second alphabetically, at least for timepoints)
  level_numerator = group_levels_combinations[2, combination_num]
  level_denominator = group_levels_combinations[1, combination_num]
  message("comparison: ", paste(group_name, ":", level_numerator, "vs", level_denominator))
  deseq2_results(deseq_dataset = dds, contrast = c(group_name, level_numerator, level_denominator))
}

# delete Rplots.pdf (generated by pheatmap)
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
