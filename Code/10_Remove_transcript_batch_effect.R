
args<-commandArgs(trailingOnly = T)

library(GWASTools)
library(dplyr)
library(tidyverse)
library(plyr)
library(data.table)
library(lme4)
library(Olivia)




phenotype_file<-args[1]


pheno<- fread(phenotype_file, data.table = F)
head(pheno)


rnaseq_count_matrix<-args[2]


rnaseq_count<-getobj(rnaseq_count_matrix)


IDs_both <- intersect(pheno$TOR_ID, colnames(rnaseq_count))
rnaseq_matrix <- rnaseq_count[, IDs_both] 

phenotypes <- pheno[match(IDs_both,pheno$TOR_ID),]
head(phenotypes)

rownames(phenotypes)<-phenotypes$TOR_ID
dim(phenotypes)
head(phenotypes)
median_norm <- median_normalization(rnaseq_matrix)
dim(median_norm)



clean_count_matrix <- apply_filters(count_matrix = median_norm, 
                                    median_min = 1, 
                                    expression_sum_min = 10, 
                                    max_min = 10, 
                                    range_min = 5, 
                                    prop_zero_max = 0.5)




log_count_matrix<- log_replace_half_min(clean_count_matrix)


list_transcripts<-args[3]

transcripts<-unlist(str_split(list_transcripts,","))


#list_transcripts<-"ENSG00000000457,ENSG00000000419"
transcripts<-unlist(str_split(list_transcripts,","))

count_matrix<- log_count_matrix[rownames(log_count_matrix) %in% transcripts,,drop=FALSE]


count_matrix<-data.frame(t(count_matrix))

head(count_matrix)
count_matrix<-count_matrix %>% rownames_to_column(var="TOR_ID")


combine_pheno<-left_join(pheno,count_matrix, by="TOR_ID")

trait<-args[4]


covariates_string<-as.character(args[5])

non_batch_effect_covars<- unlist(str_split(covariates_string,","))
#covariates_string<-"as.factor(sex),as.factor(study_site),as.factor(Shipment),as.factor(BroadUW),age,as.factor(race),as.factor(Plate)"

fixed_batch_effect<-as.character(args[6])

fixed_batch_effect<- unlist(str_split(fixed_batch_effect,","))



clean_transcripts<-matrix(NA, nrow = nrow(combine_pheno), ncol = length(transcripts))
colnames(clean_transcripts)<-transcripts

for(transcript in transcripts){
  
  
  
  fit<-random_effect_mix_model(phenotype = combine_pheno,
                               covariates_string = non_batch_effect_covars,
                               fixed_batch_effect =fixed_batch_effect,
                               trait = trait, 
                               outcome = transcript )
  
  
  batch_effect<- remove_fixed_effect(regression_model = fit,
                                     fixed_batch_effect = fixed_batch_effect)
  
  
  clean_transcripts[,transcript]<-  as.matrix(combine_pheno[,transcript])- batch_effect
  
  
  
  
}



clean_transcripts<-data.frame(clean_transcripts)

output_transcript<-args[7]

write.csv(output_transcript, file = output_transcript, row.names = F)
