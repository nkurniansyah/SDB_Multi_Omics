library(Olivia)
library(dplyr)
library(data.table)
library(EnsDb.Hsapiens.v86)
library(GWASTools)



phenotype_file<-"./SDB_phenotype.csv"
phenotype_df<-read.csv(phenotype_file)
#phenotype_annot<-AnnotatedDataFrame(phenotype_df)


## Use cluster to run this GWAS
chr<-1:22


gdsfile<-paste0("/TOPMed_chr",chr,".gds")


gds <- seqOpen(gdsfile)


filterByPass(gds)



null_model_file<- "./FAM106A_nullmodel.RData"
nullModel <- getobj(null_model_file)


annot <- matchAnnotGds(gds, phenotype_df)
seqData <- SeqVarData(gds, sampleData=annot)

# get null model

# get samples included in null model
nullModel <- GENESIS:::.updateNullModelFormat(nullModel)
sample.id <- nullModel$fit$sample.id


filterByMAF(seqData, sample.id, maf.min=0.05, build="hg38")


# create iterator
block.size <- as.integer(1024)
iterator <- SeqVarBlockIterator(seqData, variantBlock=block.size)


assoc <- assocTestSingle(iterator, nullModel,
                         test="score")


output<-paste0(".GWAS/FAM106A/FAM106A_assoc_chr",chr,".RData")
save(assoc,file = output)

seqClose(seqData)

# mem stats

