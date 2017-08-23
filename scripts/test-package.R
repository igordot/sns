#!/usr/bin/env Rscript


##
## Try loading an R package. Install to local user library if not present.
##
## usage: Rscript --vanilla test-package.R package
##


# get scripts directory (directory of this file) and load relevant functions
args_all = commandArgs(trailingOnly = FALSE)
scripts_dir = normalizePath(dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(paste0(scripts_dir, "/load-install-packages.R"))

# relevent arguments
args = commandArgs(trailingOnly = TRUE)
package_name = args[1]

# check for arguments
if (length(args) < 1) stop("not enough arguments provided")

# load relevant packages
load_install_packages(package_name)



# end
