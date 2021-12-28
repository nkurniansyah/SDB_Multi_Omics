#' Run association test
#'
#' @param pheno : data frame of the phenotype 
#' @param outcome : as.numeric , outcome to test 
#' @param covars_prs : covariates to adjust
#' @param group.var : 
#' @param covmat: specifies the covariance structures for the random effects included in the model 
#'
#' @return association test results in data frame
#' @export
#'
#' @examples
#' 

run_assoc_mixmodel<-function(pheno, outcome, covars_prs, covmat=NULL,group.var=NULL){
  
  
  if(all(pheno[,outcome] %in% c(0,1,NA))){
    family<-"binomial"
  }else{
    family="gaussian"
  }
  
  nullmod_stage1 <- fitNullModel(pheno,
                          outcome = outcome,
                          covars = covars_prs,
                          family = family,
                          cov.mat=covmat,
                          group.var = group.var,
                          verbose=TRUE)
  
  nullmod_stage2<-nullModelInvNorm(nullmod_stage1, cov.mat=cov.mat,
                                   norm.option="all",
                                   rescale="residSD")
  
  
  

  nullmod_other_fixef<- nullmod_stage2$fixef
  nullmod_other_fixef<-nullmod_other_fixef[grepl("prs",rownames(nullmod_other_fixef)),]
  sample_size<-length(nullmod$outcome)
  assoc_df<- data.frame(nullmod_other_fixef, sample_size)
  colnames(assoc_df)<- c("prs_effect","prs_se","prs_stat","pval_assoc_with_outcome","sample_size")
  assoc_df<- assoc_df %>% rownames_to_column(var="exposure")
  
  return(assoc_df)
  
}



#' Title: rank normalization
#'
#' @param x : vector and as numeric  outcome
#'
#' @return normalize value
#' @export
#'
#' @examples
rank_normalization<-function(x){
  qnorm((rank(x,na.last="keep")-0.5)/length(x))
}





#' Title: Metaboomics association test
#'
#' @param phenotype :data frame of the phenotype include all the metabolomics
#' @param sample_weight : weight of sample
#' @param psu_id  : PSU ID
#' @param strata : stratum 
#' @param covars : covars to adjust
#' @param exposure : exposure 
#' @param metab_include : vector of the name of metabolomics
#'
#' @return data frame of metabolomics association test
#' @export
#'
#' @examples
run_metab_assoc<-function(phenotype, sample_weight, psu_id, strata, covars, exposure,metab_include){
  
  survey_design=svydesign(id=as.formula(paste0("~", psu_id,collapse = "")),
                          strata=as.formula(paste0("~", strata, collapse = "")),
                          weights=as.formula(paste0("~", sample_weight, collapse = "")), 
                          data=phenotype)
  
  
  out<-list()
  
  for(metab in metab_include) {
   
    model<- svyglm(as.formula(paste0(metab,"~",exposure,"+",paste(cov,collapse= "+"))),design=survey_design)
    
   
    summary<-summary(model)$coef
    
    
    sample.size<- length(model$residuals)
    val<- summary[rownames(summary)==exposure,]
    matab_names<- colnames(metab_imp_cont)[i]
    val<- c(matab_names,sample.size,val,exposure )
    
    val_df<- data.frame(t(val))
    
    colnames(val_df) <- c("Metabolite","SampleSize" ,"Beta", "SE", "Stat", "P.value")
    

    # replace clase into numeric--> first into character first then numeric
    val_df[, 2:ncol(val_df)] <- sapply(val_df[, 2:ncol(val_df)], as.character)
    val_df[, 2:ncol(val_df)] <- sapply(val_df[, 2:ncol(val_df)], as.numeric)
    
    out[[metab]]<- val_df
    
  }
  
  
  res<- do.call(rbind, out)
  res<-data.frame(res)
  
  
  
  return(res)
  
}

