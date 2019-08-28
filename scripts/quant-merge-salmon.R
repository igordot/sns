#!/usr/bin/env Rscript


##
## Summarize Salmon transcript abundances on a gene level for individual samples and combine into a single table.
##
## usage: Rscript --vanilla quant-merge-salmon.R genes.gtf quant_sf_dir out_base
##


# increase output width
options(width = 120)
# print warnings as they occur
options(warn = 1)

# get scripts directory (directory of this file) and load relevant functions
args_all = commandArgs(trailingOnly = FALSE)
scripts_dir = normalizePath(dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(paste0(scripts_dir, "/load-install-packages.R"))

# relevant arguments
args = commandArgs(trailingOnly = TRUE)
genes_gtf = args[1]
quant_sf_dir = args[2]
out_base = args[3]

# check that inputs are valid
if (length(args) != 3) stop("usage: Rscript --vanilla quant-merge-salmon.R genes.gtf quant_sf_dir out_base")
if (!file.exists(genes_gtf)) stop("file does not exist: ", genes_gtf)
if (!dir.exists(quant_sf_dir)) stop("directory does not exist: ", quant_sf_dir)

# load relevant packages
load_install_packages("tximport")
load_install_packages("GenomicFeatures")
load_install_packages("glue")
load_install_packages("stringr")
load_install_packages("readr")

# import GTF
txdb = makeTxDbFromGFF(genes_gtf, format = "gtf")
txkey = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, txkey, "GENEID", "TXNAME")

# find quant.sf files in a specified directory
quant_files = list.files(path = quant_sf_dir, pattern = "quant.sf.gz", full.names = TRUE, recursive = TRUE)

message("num found quant files: ", length(quant_files))

# determine sample names from quant.sf file names
quant_files_names = quant_files
# renamed quant.sf file name with sample name in the file name (if files were moved)
quant_files_names = str_remove(quant_files_names, ".quant.sf.gz")
# original quant.sf file name with sample name as the directory (if files are in the default output directory)
quant_files_names = str_remove(quant_files_names, "/quant.sf.gz")
quant_files_names = str_remove(quant_files_names, ".*/")
names(quant_files) = quant_files_names

# import transcript-level estimates and summarizes to the gene-level
# "salmon" software type uses "TPM" column as abundances and "NumReads" as estimated counts
# scale using the average transcript length over samples and the library size (lengthScaledTPM)
txi = NULL
try(txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM"))

# if tximport failed and GTF transcripts did not have version decimals, try ignoring version in Salmon files
if (is.null(txi)) {
  message("tximport failed")
  if (length(str_which(tx2gene$TXNAME, "\\.")) == 0) {
    message("repeating tximport ignoring transcript version")
    txi = tximport(quant_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM",
                   ignoreTxVersion = TRUE)
  }

}

message("num imported samples: ", ncol(txi$counts))
message("num imported genes:   ", nrow(txi$counts))

# check that the counts table has a reasonable number of genes
if (nrow(txi$counts) < 1000) stop("tximport counts table too small")

# importing estimates for use with differential gene expression methods
# use the tximport argument countsFromAbundance="lengthScaledTPM" or "scaledTPM"
# use the gene-level matrix txi$counts as you would a regular count matrix ("bias corrected counts without an offset")

# extract bias corrected gene counts
txi_counts = txi$counts
txi_counts = round(txi_counts, 3)

# extract TPMs
txi_tpms = txi$abundance
txi_tpms = round(txi_tpms, 3)

# generate random string so the output files do not conflict if running multiple samples in parallel
rand_str = paste0(sample(LETTERS, 10, replace = TRUE), collapse = "")

# save tximport list (for downstream Bioconductor DGE packages such as edgeR or DESeq2)
saveRDS(txi, file = glue("temp.{rand_str}.rds"))
system(glue("mv -vf temp.{rand_str}.rds {out_base}.tximport.rds"))
Sys.sleep(1)

# save counts table
write.table(txi_counts, file = glue("temp.{rand_str}.counts.txt"), quote = FALSE, sep = "\t", col.names = NA)
system(glue("mv -vf temp.{rand_str}.counts.txt {out_base}.tximport.counts.txt"))
Sys.sleep(1)

# save counts table
write.table(txi_tpms, file = glue("temp.{rand_str}.tpms.txt"), quote = FALSE, sep = "\t", col.names = NA)
system(glue("mv -vf temp.{rand_str}.tpms.txt {out_base}.tximport.tpms.txt"))
Sys.sleep(1)


# end
