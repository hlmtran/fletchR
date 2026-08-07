# archRiSEE

Turns ArchR projects into iSEE-able SingleCellExperiments.  Also allows some control over dimension reduction (i.e. improvements on iterative LSI). 

[![Build Status](https://travis-ci.org/ttriche/archRiSEE.png?branch=master)](https://travis-ci.org/ttriche/archRiSEE)  [![codecov](https://codecov.io/gh/ttriche/archRiSEE/branch/master/graph/badge.svg)](https://codecov.io/gh/ttriche/archRiSEE)


## Installing

The devel version of the package can be pulled from GitHub with BiocManager:

    # install.packages(c("BiocManager", "remotes"))
    BiocManager::install("trichelab/fletcheR", build_vignettes=TRUE)

## Using it

The core of the package: 

    # convert an ArchR project
    library(archR)
    library(fletcheR)
    proj <- getTestProject()
    proj <- addIterativeLSI(proj, dimsToUse=1:5, varFeatures=1000, force=TRUE)
    SCE <- archRtoSCE(proj) 

    # explore it in iSEE
    iSEEarchR(SCE) 
   
Convenience functions for genomic-coordinates plotting:
 
    # plot a track 
    library(igvR)
    mcols(SCE)$aTrack <- rpois(nrow(SCE), lambda=50)
    igv <- igvR() 
    addIgvTrack(SCE, "aTrack", igv) 

The latter works best for things like pseudobulk coverage, ChromHMM, NMF, etc.


## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](https://github.com/hadley/testthat) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `test_package` will do a regular-expression pattern match within the file names. See its documentation in the `testthat` package.

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
