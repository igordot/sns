##
## Load packages. If missing, try installing from CRAN and Bioconductor.
##


load_install_packages = function(package_list) {

  for (p in package_list) {

    message("loading package: ", p)

    # create personal library
    dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

    # try to install from CRAN if package is not already installed
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      message("try installing package from CRAN: ", p)
      # install.packages(p, quiet = TRUE, repos = "https://cran.rstudio.com/")
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
