args<-commandArgs(trailingOnly = T)
library(data.table)
library(dplyr)
library(GENESIS)
library(SeqArray)



#We performed median normalization and log half min transcripts in RNA-seq counts, then we extract significant transcripts fdr_bh < 0.1.
#finaly we merge with our phenotypes. 




phenotype_file<-args[1]

pheno<- fread(phenotype_file, data.table = F)
head(pheno)


covariates_string<-as.character(args[2])

covars<- unlist(str_split(as.character(covariates_string),","))

#covars<-c("shipment","BroadUW","plate","age","sex","race","Study site","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

outcome<as.character(args[3])

cov.mat<-args[4]
cov.mat<- getobj(cov.mat)



pheno_annot<-AnnotatedDataFrame(pheno)

nullmodel_stage1<-fitNullModel(pheno, 
                              outcome= outcome, 
                              covars=covars,
                              group.var = "race",
                              family = "gaussian", 
                              cov.mat=cov.mat)



message("Refitting model with inverse normal transformation.")
message(length(nullmodel_stage1$outcome), " samples")

nullmodel_stage2 <- nullModelInvNorm(nullmodel_stage1, cov.mat=cov.mat,
                                    norm.option="all",
                                    rescale="residSD")


output_file<-args[5]
save(nullmodel_stage2, file = output_file)










