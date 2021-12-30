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
        --clump ./GWAS/$transcript/$transcrip\_combine.txt \
        --clump-snp-field SNP \
        --clump-field P \
        --out ./GWAS/$transcript/$transcrip\_clean_R_0.95_500_kb.txt \
        
    ;done
