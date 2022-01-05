## Introduction

This repository provides information and instruction on how to perform
the analysis in the manuscripts Multi-omics analysis of sleep-disorder
sleeping traits across multiple blood tissues.

Second, this repository also provides the codes we used for the analyses
(see folder “Code”).

## RNASeq Analyis

We performed Transcript wide association analysis of SDB phenotype
(AvgO2, MinO2, AHI, and multiple SDB traits) in multiple blood tissue in
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
    for trait in AvgO2 MinO2 AHI; do \ 

    bsub -J $trait -q medium -n 2 -o ./output_rna_seq/log/$trait.log \

    Rscript .Code/01_Transcript_association_analysis.R \
            ./Data/SDB_phenotype.csv \
            ./Data/MESA_RNASeq.RData \
            $trait \
            'age,sex,study_site,race,shipment,plate,BroadUW'\
            ./output_rna_seq/$trait_unadjBMI.csv; \
    done

    # Multi trait

    bsub -J Multi_SDB -q medium -n 2 -o ./output_rna_seq/log/Multi_SDB.log \

    Rscript .Code/02_Transcript_association_analysis_multi_exposure.R \
            ./Data/SDB_phenotype.csv \
            ./Data/MESA_RNASeq.RData \
            'AvgO2,MinO2,AHI' \
            'age,sex,study_site,race,shipment,plate,BroadUW'\
            ./output_rna_seq/Multi_sdb_unadjBMI.csv
            

There are 104 transcripts for un-adjusted BMI and 28 transcripts with
FDR p-value &lt; 0.1 across multiple blood tissues and traits for un. We
used all these transcripts for additional analyses.

## Genome-wide association study (GWAS)

The next step is to perform a Genome-wide Association Study for each
transcript using the MESA reference panel. We followed the guideline to
perform GWAS using
[TopmedPipeline](https://github.com/UW-GAC/analysis_pipeline "TopmedPipeline").

We first generate NullModel then performe assiaction test.
(03\_Transcript\_NullModel.R and 04\_Transcript\_GWAS.R)

see the example below how to execute the code:


    #STEP 1. NullModel


    for transcript in FAM20A .... AJUBA ; do \ 

    bsub -J $trait -q medium -n 2 -o ./GWAS/$trait/log/$transcript.log \

    Rscript .Code/03_Transcript_NullModel.R \
            ./Data/SDB_phenotype_with_trancsripts_log_transform.csv \
            'age,sex,study_site,race,shipment,plate,BroadUW,PC1..PC10'\
            $trancript \
            ./Data/TOPMed_kinship_matrix.RData \
            ./GWAS/$Trait/Data/$transcript\_NUllModel.RData; \
    done
            


    #STEP 2. GWAS

    for transcript FAM20A .... AJUBA ; do \ 

      for chr in {1..22}{
      
          bsub -J $transcript\_$chr -q medium -n 2 -o ./GWAS/$transcript/log/$transcript\_chr$chr.log \
      
                Rscript .Code/04_Transcript_GWAS.R \
                ./Data/SDB_phenotype_with_trancsripts_log_transform.csv \
                $chr \
                ./GWAS/$transcript/Data/$transcript\_NullModel.RData; \
                ./GWAS/$transcript/$transcript\_chr$chr.RData ; \

            done \
    ;done

            

After completing the GWAS, we combine the GWAS for each transcript using
QC using a clumping approach were implemented in
[PLINK](https://www.cog-genomics.org/plink/2.0/ "PLINK"). We used
clumping parameters R2=0.95, Distance 500kb, and P=1. We used HCHS/SOL
reference panel to perform clumping because we want to prioritize the
existing SNP in HCHS/SOL to generate PRS. See the example below to
perform clumping.


    for transcript in transcripts ; do \

    plink \
        --bfile HCHS_SOL_geno \
        --clump-p1 1 \
        --clump-r2 0.95 \
        --clump-kb 500 \
        --clump ./GWAS/$transcript/$transcript\_combine.txt \
        --clump-snp-field SNP \
        --clump-field P \
        --out ./GWAS/$transcript/$transcript\_clean_R_0.95_500_kb.txt \
        
    ;done

## Constructing PRS

We use each of summary statistics from each of the Tranccript GWAS after
QC to develop PRS weights for the corresponding transcripts. We used the
LDPRed2-auto model implemented in
[bigsnpr](https://privefl.github.io/bigsnpr/articles/LDpred2.html#computing-ldpred2-scores-genome-wide-1 "bigsnpr")
R package. Below are the step to construct PRS:

### Step 1: Correlation Matrix

First, we computed correlation matrix between variants. We used MESA
reference panel and same set of people to generate correlation matrix.
We first merge all the summary statistics for each tissues(PBMC,
Monocyte and T-cell) then create correlation matrix. See example below
to excute the code:


    for tiss PBMC Mono Tcell  ; do \ 

      for chr in {1..22}{
      
          bsub -J $transcript\_$chr -q big -n 4 -o ./PRS/Data/$tissue/log/$tissue\_corr_chr$chr.log \
      
                Rscript .Code/05_Correlation_matrix.R \
                $chr \
                ./Genotype/$tissue\_MESA.bed \
                ./Data/sample_to_include_$tissue.txt \
                ./Data/summary_stat_$tissue.txt \
                ./PRS/Data/$tissue/
                ./PRS/Data/$tissue/corr_chr$chr.rds; \
            done \
    ;done
      

### Step 2: LDpred2-auto: automatic model

After completed the correlation matrix, we compute weight based on set
of people that we used to generate GWAS. See example below to excute the
code:

    for transcript in transcripts ; do \

     bsub -J $transcript\_beta -q big -n 4 -o ./PRS/Data/PBMC/$transcript/log/$transcript\_beta.log \
      
          Rscript .Code/06_Compute_beta_auto.R \
                  ./GWAS/$transcript/$transcript\_clean_R_0.95_500_kb.txt\
                  ./Genotype/$tissue\_MESA.rds \
                  ./Data/sample_to_include_PBMC.txt \
                  ./PRS/Data/$tissue \
                  ./PRS/$transcript \
                  4 \
                   ./PRS/$transcript\_MESA_pred_auto_beta.txt \
        
    ;done

### Step 3: Generate PRS

Finally, We used new weighted from MESA and constructed PRS using
HCHS/SOL references panel. see example below to run PRS:

    for transcript in transcripts ; do \

     bsub -J $transcript\_beta -q big -n 4 -o ./PRS/Data/PBMC/$transcript/log/$transcript\_PRS.log \
      
          Rscript .Code/07_Construct_PRS.R \
                  ./PRS/$transcript\_MESA_pred_auto_beta.txt \
                  ./Genotype/HCHS_SOL_F5_imputed.rds \
                  4 \
                   ./PRS/$transcript\_HCHS_SOL_PRS_auto.txt \
        
    ;done

## Association test transcript PRS with SDB phenothpe in HCHS/SOL

We performed association analysis using mixed models (two stage
procedure) implemented in the
[GENESIS](https://github.com/UW-GAC/GENESIS "GENESIS") R package. Below
is an example how to run the association test.



    for transcript in transcripts ; do \

     bsub -J $transcript\_Assoc_avgSatO2 -q big -n 2 -o ./PRS_Asoc/PBMC/log/$transcript\_PRS_assoc.log \
      
          Rscript .Code/08_PRS_assoc_test.R \
                  ./Data/HCHS_SOL_sleep_pheno.txt \
                  ./PRS/$transcript\_HCHS_SOL_PRS_auto.txt \
                  ./PRS/$transcript\_MESA_pred_auto_beta.txt \
                  'Age,Sex,Center,log_weight_norm,Background,PC1,PC2,PC3,PC4,PC5' \
                  AvgSatO2 \
                  ./Data/HCHS_SOL_CovMatlist.RData \
                  ./PRS_Asoc/PBMC/$transcript\_assoc_with_AvgSatO2.RData
                 
    ;done

## Association test transcript PRS with Metabolomics HCHS/SOL
