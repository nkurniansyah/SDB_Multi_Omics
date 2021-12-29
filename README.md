## Introduction

This repository provides information and instruction how to perform the
analysis in the manuscripts Multi-omics analysis of sleep-disorder
sleeping traits across multiple blood tissues.

Second, this repository also provides the codes the we used for the
analyses (see folder “Code”).

## RNASeq Analyis

We performed association analysis of SDB phenotype using the
[Olivia](https://github.com/nkurniansyah/Olivia "Olivia") with RNA-seq
data in MESA (Multi-Ethnic Study of Atherosclerosis).

    install.packages("dplyr")
    library(dplyr)
    library(bigsnpr)
