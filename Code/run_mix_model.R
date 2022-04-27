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
run_assoc_mixmodel<-function(pheno, outcome, covars_prs, covmat=NULL,group.var=NULL,inverse_normal=TRUE ){
  
  
  if(all(pheno[,outcome] %in% c(0,1,NA))){
    family<-"binomial"
  }else{
    family="gaussian"
  }
  
  nullmod <- fitNullModel(pheno,
                          outcome = outcome,
                          covars = covars_prs,
                          family = family,
                          cov.mat=covmat,
                          group.var = group.var,
                          verbose=TRUE)
  
  
  
  if (as.logical(config["inverse_normal"]) & family == "gaussian") {
    

    
    message("Refitting model with inverse normal transformation.")
    message(length(nullmod$outcome), " samples")
    
    nullmod<-nullModelInvNorm(nullmod_stage1, cov.mat=cov.mat,
                              norm.option="all",
                              rescale="residSD")
  }
  
 
  

  
  nullmod_other_fixef<- nullmod$fixef
  nullmod_other_fixef<-nullmod_other_fixef[grepl("tPRS",rownames(nullmod_other_fixef)),]
  sample_size<-length(nullmod$outcome)
  assoc_df<- data.frame(nullmod_other_fixef, sample_size)
  colnames(assoc_df)<- c("prs_effect","prs_se","prs_stat","pval_assoc_with_outcome","sample_size")
  assoc_df<- assoc_df %>% rownames_to_column(var="exposure")
  
  return(assoc_df)
  
}