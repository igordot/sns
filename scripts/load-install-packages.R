##
## Load packages. If missing, try installing from CRAN and Bioconductor.
##


load_install_packages = function(package_list) {

  for (p in package_list) {

    message("test loading package: ", p)

    # try to install from CRAN if package is not already installed
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      message("installing package from CRAN: ", p)
      install.packages(p, quiet = TRUE, repos = "http://cran.rstudio.com")
    }

    # try to install from Bioconductor if package is not already installed
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      message("installing package from Bioconductor: ", p)
      source("https://bioconductor.org/biocLite.R")
      biocLite(p, suppressUpdates = TRUE)
    }

    # force load package
    library(p, character.only = TRUE, quietly = TRUE)

  }

}



# end
