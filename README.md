## Introduction

This repository provides information and instruction how to perform the
analysis in the manuscripts Multi-omics analysis of sleep-disorder
sleeping traits across multiple blood tissues.

Second, this repository also provides the codes the we used for the
analyses (see folder “Code”).

## RNASeq Analyis

We performed association analysis of SDB phenotype (AvgO2, MinO2, AHI
and multiple SDB traits)with RNA-seq in MESA (Multi-Ethnic Study of
Atherosclerosis) using
[Olivia](https://github.com/nkurniansyah/Olivia "Olivia") R packge.

we used 01\_Transcript\_association\_analysis.R (single trait/exposure)
and 02\_Transcript\_association\_analysis\_multi\_exposure.R (mutiple
exposure) to run the RNASeq analyis.

This scripts takes 5 arguments:  
1. Phenotype files  
2. Rnaseq count matrix (\*.RData)  
3. Trait /exposure  
4. Covariates to adjust  
5. Outpu files  

see example below to excute the code:


    Rscript .Code/01_Transcript_association_analysis.R \
            ./Data/SDB_phenotype.csv \
            ./Data/MESA_RNASeq.RData \
            AvgO2 \
            'age,sex,study_site,race,shipment,plate,BroadUW'\
            ./output_rna_seq/AvgO2_unadjBMI.csv
            
