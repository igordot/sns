#!/usr/bin/env Rscript


##
## Try loading a few packages to test R settings.
##
## usage: Rscript --vanilla test-packages.R
##


# get scripts directory (directory of this file) and load relevant functions
args_all = commandArgs(trailingOnly = FALSE)
scripts_dir = normalizePath(dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(paste0(scripts_dir, "/load-install-packages.R"))

# load some packages
test_packages = c("bitops", "mnormt", "gplots", "xtable", "ggplot2", "dplyr", "gsalib", "limma")
load_install_packages(test_packages)



# end
