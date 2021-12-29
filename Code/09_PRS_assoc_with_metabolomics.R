library(survey)
library(pastecs)
library(dplyr)
library(tidyverse)
library(GWASTools)
library(plyr)
library(factoextra)
library(purrr)
library(qvalue)
library(scales)
library(data.table)
source("Code/Utils.R")



metabolomics_files<-args[1]

metabolomics_df<-fread(metabolomics_files, data.table = F)
metabolomics_df[,3:ncol(metabolomics_df)]<-sapply(metabolomics_df[,3:ncol(metabolomics_df)], as.numeric)

pheno<-read.csv(pheno_file_name,header = TRUE, sep = ",")
head(pheno)


prs_file<-args[3]

prs_df<- fread(prs_file, data.table = F, header = T)
head(prs_df)

colnames(prs_df)<-c("sample.id","prs_auto")


prs_df[,"prs_auto"]<-(prs_df[,"prs_auto"]-mean(prs_df[,"prs_auto"]))/sd(prs_df[,"prs_auto"])


covariates_string<-args[4]

cov<- unlist(str_split(as.character(covariates_string),","))


clean_pheno<-pheno[,c("sample.id","ID","WEIGHT_FINAL_NORM_OVERALL","STRAT","PSU_ID",cov)]

clean_pheno<-clean_pheno[!is.na(clean_pheno$PSU_ID),]


complete_pheno<- full_join(clean_pheno, prs_df, by="sample.id")


# Drop samples with more than ...% missing (25%)
drop_sample_ind<-which(rowSums(is.na(metabolomics_df)) < (ncol(metabolomics_df)-2)* 0.25)
metabolomics_df<-metabolomics_df[drop_sample_ind,]


# Drop metabolites with more than ...% missing (75 %)
drop_metab_ind <- which(colSums(is.na(metabolomics_df)) < nrow(metabolomics_df)*0.75)
metabolomics_df<-metabolomics_df[,drop_metab_ind]

#rannk normalize the metab value

metabolomics_df[,3:ncol(metabolomics_df)]<-sapply(metabolomics_df[,3:ncol(metabolomics_df)], rank_normalization)


combine_pheno<- left_join(complete_pheno,metabolomics_df, by= "ID")
dim(combine_pheno)


list_metab<-colnames(metabolomics_df)[2:ncol(metabolomics_df)]



metab_assoc<-run_metab_assoc(phenotype = combine_pheno,
                             sample_weight = "WEIGHT_FINAL_NORM_OVERALL",
                             psu_id = "PSU_ID",
                             strata = "STRAT",
                             covars = covars,
                             exposure = "prs_auto",
                             metab_include =list_metab )


# remove unknow metabolite
identify_metab<-metab_assoc%>% dplyr::filter(!str_detect(metabolite,"^X"))

## Add FDR

assoc_clean_df<- identify_metab %>% mutate(FDR_BH=p.adjust(P.value,method="BH"))


metab_output<-args[]
write.csv(assoc_clean_df, file=metab_output, row.names = F)

