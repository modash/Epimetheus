if (!require("limma",character.only = TRUE))
{
source("http://bioconductor.org/biocLite.R")
biocLite("limma")
if(!require("limma",character.only = TRUE)) stop("Package not found")
}

