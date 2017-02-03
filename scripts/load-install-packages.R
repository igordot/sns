##
## Load packages. If missing, try installing from CRAN and Bioconductor.
##


load_install_packages = function(package_list) {

  for (p in package_list) {

    message("loading package: ", p)

    # create local user library path (not present by default)
    dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

    # try to install from CRAN if package is not already installed
    # use require() instead of installed.packages() to check that package is both installed and usable
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      message("try installing package from CRAN: ", p)
      # install to local user library path (destination needs to be specified the first time)
      install.packages(p, lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/", quiet = TRUE)
    }

    # try to install from Bioconductor if package is not already installed
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      message("package not available on CRAN: ", p)
      message("try installing package from Bioconductor: ", p)
      source("https://bioconductor.org/biocLite.R")
      biocLite(p, suppressUpdates = TRUE)
    }

    # package still not available (after CRAN and Bioconductor)
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      message("package not available on Bioconductor: ", p)
    }

    # force load package
    library(p, character.only = TRUE, quietly = TRUE)

  }

}



# end
