#!/usr/bin/env Rscript


##
## Try loading (installing if not available) a few R packages.
##
## usage: Rscript --vanilla test-packages.R
##


# get scripts directory (directory of this file) and load relevant functions
args_all = commandArgs(trailingOnly = FALSE)
scripts_dir = normalizePath(dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(paste0(scripts_dir, "/load-install-packages.R"))

# load some packages (some will probably need to be installed)
# used by GATK: ggplot2, gplots, reshape, gsalib
test_packages = c("bitops", "mnormt", "reshape", "gplots", "ggplot2", "gsalib", "dplyr", "limma")
load_install_packages(test_packages)



# end
