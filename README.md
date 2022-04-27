## Introduction

This repository provides information regarding transcriptome
wide-association test and the construction of a transcript polygenic
risk score (tPRS) for sleep disorder breathing (SDB) that we developed
in manuscript “An integrated multi-omics analysis of sleep-disordered
breathing traits across multiple blood cell types” (link to be added).

First, it provides instructions for constructing the tPRS based on
transcript GWAS. We provide the relevant for each transcript summary
statistics (see folder “Summary\_Statistics\_for\_tPRS\_construction”),
as well as code for using them to construct the tPRS.

Second, this repository also provides code the we used for the analyses
in the manuscript (see folder “Code”).

## Transcriptome wide-association test

We performed Transcript wide association analysis of SDB phenotype
(AvgO2, MinO2, AHI, and multiple SDB traits) in multiple blood tissue in
MESA (Multi-Ethnic Study of Atherosclerosis) using
[Olivia](https://github.com/nkurniansyah/Olivia "Olivia") R package that
we developed our previous works in manuscript [Benchmarking association
analyses of continuous exposures with RNA-seq in observational
studies](https://academic.oup.com/bib/article-abstract/22/6/bbab194/6278609).

Briefly, we performed median normalization for gene-level expected
counts, and then filtered lowly expressed gene transcripts defined by
removing transcripts with proportion of zero higher than 0.5, median
value lower than 1, maximum expression range value lower than 5, and
maximum expression value lower than 10. Transcript counts were log
transformed after counts of zero were replaced with half the minimum of
the observed transcript count in the sample The analyses were adjusted
for age, sex, study center, race/ethnic group, and batch variables:
plates, shipment batch, and study site. Because BMI is a strong risk
factor for SDB and is assumed to be part of the causal chain, we
conducted additional analyses adjusting for BMI (BMIadj). We computed
empirical p-values to account for the highly skewed distribution of SDB
phenotypes which implemented in Olivia R package, which may lead to
false negative associations if ignored. Finally, we accounted for
multiple testing by applying False-Discovery Rate (FDR) correction to
each of the association analyses using the Benjamini-Hochberg (BH)
procedure. We carried forward transcript associations with FDR
p-value&lt;0.1

## Transcript genome-wide association study (GWAS)

The next step is to perform a Genome-wide Association Study for each
transcript with FDR p-value &lt; 0.1 (96 transcript unadjusted BMI and
24 transcriptadjusted BMI analysis) using the MESA reference panel. We
followed the guideline to perform GWAS using
[TopmedPipeline](https://github.com/UW-GAC/analysis_pipeline "TopmedPipeline").

## Constructing tPRS

We used [PRSice 2.3.1.e](https://www.prsice.info "PRSice 2.3.1.e") to
construct tPRS for each transcript GWAS. This command below is to
construct tPRS using the transcript GWAS that we provide
(./Summary\_Statistics\_for\_tPRS\_construction). No clumping is needed
and no selection of SNPs. The summary statistics are already based on
the specific set of SNPs selected after validation step in idependent
data set(see manuscript for the detail).Note that genetic data files
need to be specified in the –target argument.



    Rscript ./PRSice.R \
     --dir ./PRS_Output \
     --prsice ./PRSice_linux/PRSice_linux \
     --base ./Summary_Statistics_for_tPRS_construction/. \
     --target ./Genotype \
     --thread 2 \
     --chr Chromosome 
     --bp Position 
     --A1 Allele1 
     --A2 Allele2 
     --pvalue PValue \
     --bar-levels Threshold \
     --stat BETA 
     --all-score T \
     --out ./out_prs \
     --no-clump T
     --print-snp T \
     --ignore-fid T 
     --no-regress T 
     --fastscore T 
     --model add 
     --no-full T 
     --chr-id c:l:a:b

## Example code for association analsis

We performed association analysis for each tPRS using mixed models
implemented in the GENESIS R package. Below is an example code. It uses
function that we provide in the folder “Code”.

    library(GENESIS)
    library(GWASTools)
    library(data.table)

    source("./Code/*")


    #phenotype

    pheno<- fread(phenotype_file, data.table=F)


    # merge tPRS with phenotype

    pheno_df<-left_join(pheno,tPRS_transcript1, by="sample.id" )


    covarites_prs<- c("age","sex","site","race",paste0("PC_",1:5),"tPRS")

    #example outcome
    outcome<-"AvgO2"

    ## Kinship matrix

    covMatlist<-getobj(covMatlist)


    assoc_df<- run_assoc_mixmodel(pheno=pheno,
                                  outcome=outcome,
                                  covars_prs=covarites_prs, 
                                  covmat=covMatlist,
                                  group.var=NULL)

## Example code for metebolome wide-association test

We performed metebolome wide-association test for each tPRS with
metabolomics data. see example below:
