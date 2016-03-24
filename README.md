# cancerSeqStudy

Identifying genes with more mutations then expected has been central methodology for identifying putative cancer driver genes in exome sequencing studies of cancer samples. Identifying significantly mutated genes (SMG) fundamentally relies on estimating a background mutation rate. Mutation rate varies over more than 2 orders of magnitude providing a substantial statistical estimation challenge. Analysis not accounting for the uncertainty in mutation rate yields overly optimistic assessments. In this package, we examine statistical power (either with known or uncertain mutation rate) and false positives induced by unaccounted variation in mutation rate.

## Installation

[![Build Status](https://travis-ci.com/ctokheim/cancerSeqStudy.svg?token=KhnctpTdxNuuZ9Z1kcsg&branch=master)](https://travis-ci.com/ctokheim/cancerSeqStudy)

Install this package via github using the `devtools` package.

```R
> devtools::install_github('ctokheim/cancerSeqStudy')
```
