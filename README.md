# cancerSeqStudy

Identifying genes with more mutations then expected has been central methodology for identifying putative cancer driver genes in exome sequencing studies of cancer samples. Identifying significantly mutated genes (SMG) fundamentally relies on estimating a background mutation rate. Mutation rate varies over more than 2 orders of magnitude providing a substantial statistical estimation challenge. However, recent methods have taken an alternative approach known as "ratio-metric" method. Ratio-metric methods examine specific compositions of mutations normalized by the total number of mutations occurring in the gene. Regardless of methodology, analysis not accounting for the uncertainty in mutation parameters yields overly optimistic assessments. In this package, we examine statistical power (either with known or uncertain mutation rate) and false positives induced by unaccounted variation in mutation rate.

## Documentation

A vignette describing the usage of this package is available on github, [here](https://github.com/KarchinLab/cancerSeqStudy/blob/master/vignettes/cancerSeqStudy.Rmd).

## Installation

[![Build Status](https://travis-ci.org/KarchinLab/cancerSeqStudy.svg?branch=master)](https://travis-ci.org/KarchinLab/cancerSeqStudy)

Install this package via github using the `devtools` package.

```R
> devtools::install_github('KarchinLab/cancerSeqStudy')
```
