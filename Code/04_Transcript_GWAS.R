library(Olivia)
library(dplyr)
library(data.table)
library(EnsDb.Hsapiens.v86)
library(GWASTools)



phenotype_file<-args[1]

phenotype_df<-read.csv(phenotype_file)
#phenotype_annot<-AnnotatedDataFrame(phenotype_df)


## Use cluster to run this GWAS
chr<-args[2]


gdsfile<-paste0("/TOPMed_chr",chr,".gds")


gds <- seqOpen(gdsfile)


filterByPass(gds)



null_model_file<- args[3]
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


output<-args[4]
save(assoc,file = output)

seqClose(seqData)

# mem stats

