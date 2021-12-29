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
    metab_res<- summary[rownames(summary)==exposure,]
    matab_names<- as.character(metab)
    metab_res<- c(matab_names,sample.size,metab_res)
    
    metab_res<- data.frame(t(metab_res))
    
    colnames(metab_res) <- c("Metabolite","SampleSize" ,"Beta", "SE", "Stat", "P.value")
    

    # replace clase into numeric--> first into character first then numeric
    metab_res[, 3:ncol(metab_res)] <- sapply(metab_res[, 3:ncol(metab_res)], as.character)
    metab_res[, 3:ncol(metab_res)] <- sapply(metab_res[, 3:ncol(metab_res)], as.numeric)
    
    out[[metab]]<- metab_res
    
  }
  
  
  res<- do.call(rbind, out)
  res<-data.frame(res)
  
  
  
  return(res)
  
}




#' Title : Run mix model for random effect or lmer
#'
#' @param phenotype: data frame of the phenotypes, include all proteomics or transcriptomics (outcome)
#' @param covariates_string : covariates to adjust (without batch effect)
#' @param trait : trait (SDB) exposure
#' @param fixed_batch_effect : vector Fix batch effect 
#' @param random_batch_effect " vector of random effecr
#'
#' @return regression model
#' 
#' @export
#'
#' @examples
#' 
random_effect_mix_model<-  function(phenotype, covariates_string, trait, fixed_batch_effect, random_batch_effect=NULL, outcome){
  
  fixed_covariates<- c(covariates_string, trait, fixed_batch_effect)
  fixed_covariates<-paste(outcome,"~",paste(fixed_covariates,collapse = "+" ))

  if(!is.null(random_batch_effect)){
    random_effect<- paste('(1|',random_batch_effect , ')', collapse = "+")
    
  }
  
  fit.warn <- tryCatch(
    {list(lmer(formula(paste(fixed_covariates, random_effect, sep = " + ")), data =phenotype ),"MixedOK")},
    warning = function(Warn){
      print(paste("MY_WARNING:  ",Warn))
      fit <- lm(fixed_covariates , data = phenotype)
      return (list(fit,"Warn"))},
    error = function(err) 
    {print(paste("MY_ERROR:  ",err))
      fit <- lm(fixed_covariates, data = phenotype)
      return (list(fit,"err"))})
  
  fit <- fit.warn[[1]]
  return(fit)
  
}


#' Title Removing Fixed effect 
#'
#' @param regression_model : regression model (routput of random_effect_mix_model)
#' @param fixed_batch_effect : vector of fixed batch effect
#'
#' @return value of fix batch effect to remove
#' @export
#'
#' @examples
#' 
remove_fixed_effect <- function(regression_model,fixed_batch_effect ){
  
  #whicgrepl(batch, colnames(model_mat))
  model_mat<- model.matrix(regression_model)
  out<-list()
  for(batch in fixed_batch_effect){
    message(paste0("remove ", batch))
    fixed_effect_value <- as.matrix(model_mat[,grepl(batch, colnames(model_mat),fixed = TRUE)])
    colnames(fixed_effect_value)<- colnames(model_mat)[grepl(batch, colnames(model_mat), fixed = T)]
    fixed_effects_estimate <- summary(regression_model)$coefficients[grepl(batch,rownames(summary(regression_model)$coefficients), fixed = TRUE),"Estimate"]
    names(fixed_effects_estimate)<- rownames(summary(regression_model)$coefficients)[grepl(batch,rownames(summary(regression_model)$coefficients), fixed = TRUE)]
    fixed_effect_value<-fixed_effect_value[,names(fixed_effects_estimate)]
    fixed_batch_effect_to_remove<-as.matrix(fixed_effect_value)%*%fixed_effects_estimate
    out[[batch]]<-fixed_batch_effect_to_remove
    
  }
  
  fixed_effect_to_remove<--Reduce(`+`, out)
  
  return(fixed_effect_to_remove)
  
}




#' Title Removing random effect
#'
#' @param regression_model : regression model (routput of random_effect_mix_model)
#' @param random_batch_effect :vector of random batch effect
#' @param phenotype : : data frame of the phenotypes, include all proteomics or transcriptomics (outcome)
#'
#' @return value of random batch effect to remove
#' 
#' @export
#'
#' @examples
remove_random_effect <- function(regression_model,random_batch_effect, phenotype){
  
  out<-list()
  for(random_effect in random_batch_effect){
    random_eff_runef<- ranef(regression_model)[[random_effect]][match(phenotype[[random_effect]], rownames(ranef(regression_model)[[random_effect]])),1]
    
    out[[random_effect]]<- random_eff_runef
  }
  

  random_batch_effect_to_remove<-Reduce(`+`, out)

  return(random_batch_effect_to_remove)
}


