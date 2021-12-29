---
title: "Multi-Omics Sleep Disorder Breathing"
author: "Tamar Sofer & Nuzulul Kurniansyah"
date: "2/18/2021"
output: 
  pdf_document:
    toc: true
  html_document:
    toc: true
---


## Introduction

This repository provides information and instruction how to perform the
analysis in the manuscripts Multi-omics analysis of sleep-disorder
sleeping traits across multiple blood tissues.

Second, this repository also provides the codes the we used for the
analyses (see folder “Code”).

## Required packages

We used [PRSice 2.3.1.e](https://www.prsice.info "PRSice 2.3.1.e") to
generate PRS. We provide example code that also uses PRSice to construct
PRS based on the provided summary statistics in folder
“Summary\_Statistics\_for\_PRS\_construction”.

Other software and packages that we used, but may not be necessary for
others to construct the PRS, are as follows:  
1. We performed the analysis using R version 4.0.2.  
2. We used the following packages from CRAN: dplyr, tidyverse,
data.table, purrr, pROC.  
3. We used the following packages from BioConductor: GENESIS,
GWASTools.  

    install.packages("dplyr")
    library(dplyr)
    library(bigsnpr)
