#!/usr/bin/env Rscript


##
## Differential gene expression with DESeq2.
##
## usage: Rscript --vanilla dge-deseq2.R counts_table.txt groups_table.csv
##


# increase output width
options(width = 150)

# java heap size
options(java.parameters="-Xmx8G")

# get scripts directory (directory of this file) and load relevant functions
args_all = commandArgs(trailingOnly = FALSE)
scripts_dir = normalizePath(dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(paste0(scripts_dir, "/load-install-packages.R"))
source(paste0(scripts_dir, "/deseq2-pca.R"))
source(paste0(scripts_dir, "/deseq2-results.R"))

# relevent arguments
args = commandArgs(trailingOnly = TRUE)
counts_table_file = args[1]
groups_table_file = args[2]

# check that input files exist
if (!file.exists(counts_table_file)) stop("file does not exist: ", counts_table_file)
if (!file.exists(groups_table_file)) stop("file does not exist: ", groups_table_file)

# load relevant packages
load_install_packages("DESeq2")
load_install_packages("RColorBrewer")
load_install_packages("gplots")
load_install_packages("genefilter")
load_install_packages("lattice")
load_install_packages("xlsx")

message(" ========== prepare inputs ========== ")

# import groups table
groups_table = read.csv(file = groups_table_file, header = TRUE, row.names = 1, colClasses = "factor")
message("groups table sample num: ", nrow(groups_table))
message("groups table groups: ", colnames(groups_table))

# import counts table
counts_table = read.delim(file = counts_table_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
message("full counts table gene num: ", nrow(counts_table))
message("full counts table sample num: ", ncol(counts_table))

# subset to samples in groups table (also sets samples to be in the same order)
counts_table = counts_table[,rownames(groups_table)]
message("group-subset counts table gene num: ", nrow(counts_table))
message("group-subset counts table sample num: ", ncol(counts_table))

# group info
group_name = colnames(groups_table)[1]
message("group name: ", group_name)
group_levels = levels(groups_table[,group_name])
message("group levels: ", toString(group_levels))

# design formula
design_formula = formula(paste("~", group_name))
message("design formula: ", design_formula)

message(" ========== normalization ========== ")

# import and normalize
dds = DESeqDataSetFromMatrix(countData = counts_table, colData = groups_table, design = design_formula)
dds = DESeq(dds, parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 4))
save(dds, file = "deseq2.dds.RData")

vsd = varianceStabilizingTransformation(dds)
save(vsd, file = "deseq2.vsd.RData")

# export counts
write.csv(counts(dds, normalized = FALSE), file = "counts.raw.csv")
write.csv(round(counts(dds, normalized = TRUE), digits = 3), file = "counts.norm.csv")
write.xlsx2(x = round(counts(dds, normalized = TRUE), digits = 3), file = "counts.norm.xlsx", sheetName = "normalized counts")

message(" ========== QC ========== ")

png("plot.sparsity.png", width = 6, height = 6, units = "in", res = 300)
print(plotSparsity(dds, normalized = TRUE))
dev.off()

pdf("plot.pca.pdf", width = 7, height = 5, family = "Palatino", pointsize = 10)
deseq2_pca(vsd, intgroup = c("group"), ntop = 1000)
dev.off()

png("plot.pca.png", width = 7, height = 5, units = "in", res = 300)
deseq2_pca(vsd, intgroup = c("group"), ntop = 1000)
dev.off()

message(" ========== differential expression ========== ")

# perform comparisons for all combinations of group levels
group_levels_combinations = combn(group_levels, m = 2, simplify = TRUE)
for (combination_num in 1:ncol(group_levels_combinations)) {
  level_numerator = group_levels_combinations[1, combination_num]
  level_denominator = group_levels_combinations[2, combination_num]
  message("comparison: ", paste(group_name, ":", level_numerator, "vs", level_denominator))
  deseq2_results(DESeqDataSet = dds, contrast = c(group_name, level_numerator, level_denominator))
}



# end