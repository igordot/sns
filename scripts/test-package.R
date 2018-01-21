#!/usr/bin/env Rscript


##
## Try loading an R package. Install to local user library if not present.
## This standalone script allows for basic package management without modifying any code that requires those packages.
##
## usage: Rscript --vanilla test-package.R package
##


# get scripts directory (directory of this file) to add load_install_packages() function
args_all = commandArgs(trailingOnly = FALSE)
scripts_dir = normalizePath(dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(paste0(scripts_dir, "/load-install-packages.R"))

# relevant arguments
args = commandArgs(trailingOnly = TRUE)
package_name = args[1]

# check for arguments
if (length(args) < 1) stop("not enough arguments provided")

# load relevant packages
message("checking R package: ", package_name)
load_install_packages(package_name)



# end
