##
## Load packages. If missing, try installing from CRAN and Bioconductor.
##


load_install_packages = function(package_list) {

  for (p in package_list) {

    message("loading package: ", p)

    # create local user library path (not present by default)
    dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

    # try to install using biocLite() (installs both Bioconductor and CRAN packages)
    # use require() instead of installed.packages() to check that package is both installed and usable
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      message("installing package: ", p)
      source("https://bioconductor.org/biocLite.R")
      biocLite(p, suppressUpdates = TRUE, lib = Sys.getenv("R_LIBS_USER"))
    }

    # package still not available (after CRAN and Bioconductor)
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      message("package not available on CRAN or Bioconductor: ", p)
    }

    # load package
    suppressPackageStartupMessages(library(p, character.only = TRUE))

  }

}



# end
