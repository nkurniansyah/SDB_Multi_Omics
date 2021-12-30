## Introduction

This repository provides information and instruction how to perform the
analysis in the manuscripts Multi-omics analysis of sleep-disorder
sleeping traits across multiple blood tissues.

Second, this repository also provides the codes the we used for the
analyses (see folder “Code”).

## RNASeq Analyis

We performed Transcript wide association analysis of SDB phenotype
(AvgO2, MinO2, AHI and multiple SDB traits) in multiple blood tissue in
MESA (Multi-Ethnic Study of Atherosclerosis) using
[Olivia](https://github.com/nkurniansyah/Olivia "Olivia") R package.

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


    # Single trait


    Rscript .Code/01_Transcript_association_analysis.R \
            ./Data/SDB_phenotype.csv \
            ./Data/MESA_RNASeq.RData \
            AvgO2 \
            'age,sex,study_site,race,shipment,plate,BroadUW'\
            ./output_rna_seq/AvgO2_unadjBMI.csv
            

    # Multi trait

    Rscript .Code/02_Transcript_association_analysis_multi_exposure.R \
            ./Data/SDB_phenotype.csv \
            ./Data/MESA_RNASeq.RData \
            'AvgO2,MinO2,AHI' \
            'age,sex,study_site,race,shipment,plate,BroadUW'\
            ./output_rna_seq/Multi_sdb_unadjBMI.csv

There are 104 transcripts for un-adjusted BMI and 28 transcripts with
FDR p-value &lt; 0.1 accross multiple blood tissues and traits for un.
We used all these transcript for additional analyses.

## Genome-wide association study (GWAS)

Next step is perform Genome-wide Assocition Study for each transcript
using MESA refrence panel. We followed the guideline to perform GWAS
using
[TopmedPipeline](https://github.com/UW-GAC/analysis_pipeline "TopmedPipeline").

We first generate NullModel then performe assiaction test.
(03\_Transcript\_NullModel.R and 04\_Transcript\_GWAS.R)
