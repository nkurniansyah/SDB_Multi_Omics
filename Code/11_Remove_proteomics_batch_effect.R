args<-commandArgs(trailingOnly = T)
library(GWASTools)
library(dplyr)
library(tidyverse)
library(plyr)
library(data.table)
library(lme4)


proteomics_file<-args[1]



proteomics_dat_df<-getobj(proteomics_file)
#change proteomics into numeric
proteomics_dat_df[,1:ncol(proteomics_dat_df)]<-sapply(proteomics_dat_df[,1:ncol(proteomics_dat_df)],as.character)

proteomics_dat_df[,1:ncol(proteomics_dat_df)]<-sapply(proteomics_dat_df[,1:ncol(proteomics_dat_df)],as.numeric)


#normalize to log transform
proteomics_dat_df[,1:ncol(proteomics_dat_df)]<-sapply(proteomics_dat_df[,1:ncol(proteomics_dat_df)],log)
  
phenotype<-args[2]
pheno<- fread(phenotype, data.table = F)



protein_list<-args[3]

proteins<-as.character(str_split(protein_list,","))
log_proteomics<-proteomics_dat_df[rownames(proteomics_dat_df) %in% proteins ,, drop=FALSe]



log_proteomics<-data.frame(t(log_proteomics))

head(log_proteomics)
log_proteomics<-log_proteomics %>% rownames_to_column(var="TOP_ID")


combine_pheno<-left_join(pheno,log_proteomics, by="TOP_ID")

trait<-args[4]


#covariates_string<-"gender1,site1c,race1c,age5c,Subarray"
covariates_string<-as.character(args[5])

non_batch_effect_covars<- unlist(str_split(covariates_string,","))
#covariates_string<-"as.factor(sex),as.factor(study_site),as.factor(Shipment),as.factor(BroadUW),age,as.factor(race),as.factor(Plate)"

fixed_batch_effect<-as.character(args[6])

fixed_batch_effect<- unlist(str_split(fixed_batch_effect,","))


random_batch_effect<-args[7]

random_batch_effect<- unlist(str_split(random_batch_effect,","))


clean_proteins<-matrix(NA, nrow = nrow(combine_pheno), ncol = length(proteins))
colnames(clean_proteins)<-proteins


for(protein in proteins){
  
  
  model<-random_effect_mix_model(phenotype = combine_pheno,
                               random_batch_effect=random_batch_effect,
                               covariates_string = non_batch_effect_covars,
                               fixed_batch_effect =fixed_batch_effect,
                               trait = trait, 
                               outcome = protein )
  
  fit<-model[[1]]
  fixed_batch_effect<- remove_fixed_effect(regression_model = fit,
                                          fixed_batch_effect = fixed_batch_effect)
  
  
  clean_proteins[,protein]<-  as.matrix(combine_pheno[,protein])- fixed_batch_effect
  
  
  if (model[[2]] == "MixedOK") {
    message("remove random effect")
    random_batch_effect<-remove_random_effect(regression_model = fit,
                                              random_batch_effect = random_batch_effect,
                                              phenotype = combine_pheno)
    clean_proteins[,protein]<-  clean_proteins[,protein]- random_batch_effect
    
  }
  
  
  
  
  
  
}
clean_proteins_df<-data.frame(clean_proteins)
clean_proteins_df$TOP_ID<-combine_pheno$TOP_ID

out_file<-config["out_file"]
write.csv(clean_proteins_df, file = out_file, row.names = F)

