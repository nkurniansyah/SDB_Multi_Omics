

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

