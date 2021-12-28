
args<-commandArgs(trailingOnly = T)
library(data.table)
library(dplyr)
library(GWASTools)
library(GENESIS)
library(SeqArray)

source("Code/Utils.R")

phenotype_file<-args[1]

pheno<-fread(phenotype_file,data.table = F)

head(pheno)



prs_file<- args[2]


prs_df<- fread(prs_file, data.table = F)

head(prs_df)

colnames(prs_df)<-c("sample.id","prs_auto")

prs_df[,"prs_auto"]<- (prs_df[,"prs_auto"]-(mean(prs_df[,"prs_auto"])))/sd(prs_df[,"prs_auto"])

beta_file<-config[3]

beta_df<- fread(beta_file, data.table = F)

all_pheno<- left_join(pheno,prs_df, by="sample.id" )
head(all_pheno)


covars<-args[4]

cov<- unlist(str_split(as.character(covars),","))

cov_prs<-c(cov,"prs_auto")



outcome<-args[5]


pheno_annot <- AnnotatedDataFrame(all_pheno)


print(outcome)

if(all(all_pheno[,outcome] %in% 0:1) ){
  family<-"binomial"
}else{
  family<-"gaussian"
}

print(family)

covMatList<- args[6]

assoc<-run_assoc_mixmodel(pheno=pheno_annot,
                          outcome=outcome, 
                          covars_prs=cov_prs,
                          covmat=covMatList,
                          group.var="race")

output_file<-args[7]

save(assoc, file = output_file)


