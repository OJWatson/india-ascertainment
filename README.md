
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.1.2-brightgreen.svg)](https://cran.r-project.org/)
[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)

## Research compendium for COVID-19 analysis in India

This is a working R compendium (think R package but for reproducible
analysis). The analysis directory contains R scripts used to generate
the results.

### Installation

    git clone https://github.com/mrc-ide/india-ascertainment.git
    cd india-rascertainment
    open india-rascertainment.Rproj
    devtools::install_deps()

### Overview

The structure within analysis is as follows:

    analysis/
        |
        ├── 01_xxxxx /           # analysis scripts used for generating figures
        |
        ├── figures/              # location of figures produced by the analysis scripts
        |
        ├── data/
        │   ├── DO-NOT-EDIT-ANY-FILES-IN-HERE-BY-HAND
        │   ├── raw_data/       # data obtained from elsewhere
        │   └── derived_data/   # data generated during the analysis
        |
        ├── src/                # orderly files

### Compendium DOI:

<http://dx.doi.org/xxxxxxx>

The files at the URL above will generate the results as found in the
publication.

### The R package

This repository is organized as an R package. There are no/negligable R
functions exported in this package - the majority of the R code is in
the analysis and src directory. The R package structure is here to help
manage dependencies, to take advantage of continuous integration, and so
we can keep file and data management simple.

To download the package source as you see it on GitHub, for offline
browsing, use this line at the shell prompt (assuming you have Git
installed on your computer):

``` r
git clone https://github.com/OJWatson/india-ascertainment.git
```

Once the download is complete, open the `india-ascertainment.Rproj` in
RStudio to begin working with the package and compendium files. We will
endeavour to keep all package dependencies required listed in the
DESCRIPTION. This has the advantage of allowing
`devtools::install_dev_deps()` to install the required R packages needed
to run the code in this repository

### Licenses

Code: [MIT](http://opensource.org/licenses/MIT) year: 2021, copyright
holder: OJ Watson

Data: [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse
