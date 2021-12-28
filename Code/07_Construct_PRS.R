args<-commandArgs(trailingOnly = TRUE)
library(bigsnpr)
library(GWASTools)
library(dplyr)
library(data.table)
library(bigsnpr)
library(argparser)



beta_file<-args[1]

beta_df<- bigreadr::fread2(beta_file, data.table = F)
#chr	pos	a0	a1	BETA	pvalue	ID

colnames(beta_df)<- c("chr","pos","a0","a1","beta", "pvalue","rsid")

# assume all chr are combine into 1 single bed files

bed_file <- args[2]

bed_rds <- paste(bed_file,sep="")

obj.bigSNP <- snp_attach(bed_rds)

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
sample.id  <- obj.bigSNP$fam


map <- obj.bigSNP$map[-3]
head(map)
names(map) <- c("chr", "rsid", "pos", "a0", "a1")


snp_df_matched<- snp_match(beta_df, map)

ncores<-as.integer(args[3])

G2 <- snp_fastImputeSimple(G, ncores=ncores)

message("genotyped imputed done")

pred_auto <- big_prodMat(G2, snp_df_matched[["beta"]],ind.col = snp_df_matched[["_NUM_ID_"]])

message("done pred auto")

sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_pred_auto <- rowMeans(pred_auto[, keep])


prs <- cbind(obj.bigSNP$fam$sample.ID,final_pred_auto)
colnames(prs) <- c("sample.id","prs_auto")


output_file<-args[4]


write.table(prs, file=output_file, row.names = F, quote = F, sep="\t")


