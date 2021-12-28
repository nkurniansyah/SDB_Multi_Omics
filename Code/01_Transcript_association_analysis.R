args<-commandArgs(trailingOnly = T)
library(Olivia)
library(dplyr)
library(data.table)
library(EnsDb.Hsapiens.v86)
library(GWASTools)



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



#outcomes<-c("AvgO2","MinO2","AHI")

outcome<-as.character(args[3])

covariates_string<-as.character(args[4])

covars<- gsub(",","+",as.character(covariates_string))
#covariates_string<-"as.factor(sex),as.factor(study_site),as.factor(Shipment),as.factor(BroadUW),age,as.factor(race),as.factor(Plate)"

head(phenotypes)
dim(phenotypes)

set.seed(123)
head(phenotypes)


message(paste0("outcome: ",outcome," and covar to adjust: ", covars))

quantile_emp_trascript<-lm_count_mat_emp_pval(count_matrix = clean_count_matrix, 
                                              pheno=phenotypes, 
                                              trait=outcome,
                                              covariates_string=covars,
                                              n_permute=100,
                                              log_transform = "log_replace_half_min",
                                              outcome_type ="continuous",
                                              gene_IDs=NULL)

head(quantile_emp_trascript)

add_annotation<-function(deg_res){
  
  gene_symbol<- select(EnsDb.Hsapiens.v86, 
                       keys =as.character(deg_res$geneID) ,
                       keytype = "GENEID",
                       columns = c("GENEID", "GENENAME"))
  
  colnames(gene_symbol)<- c("geneID","geneName")
  
  annot_deg<-left_join(deg_res,gene_symbol, by="geneID")
  annot_deg<- annot_deg %>% dplyr::rename(IDs=geneID,
                                          geneID= geneName)
  return(annot_deg)
  
}


res<- add_annotation(quantile_emp_trascript)
head(res)
output_file<-args[5]
write.csv(res, file =output_file, row.names = F )

