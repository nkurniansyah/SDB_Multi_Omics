args<-commandArgs(trailingOnly = TRUE)
library(bigsnpr)
library(GWASTools)
library(dplyr)
library(data.table)
library(bigsnpr)
library(argparser)



summary_stat_file<-args[1]

summary_stat_df<- bigreadr::fread2(summary_stat_file, data.table = F)



summary_stat_df<-summary_stat_df %>% dplyr::select(Chromosome,SNP,Position,Allele1,Allele2,N,Est.SE,BETA,PValue)

colnames(summary_stat_df)<- c("chr","rsid","pos","a0","a1","n_eff", "beta_se","beta", "pvalue")

# assume all chr are combine into 1 single bed files

bed_file <- args[2]

bed_rds <- paste(bed_file,sep="")

obj.bigSNP <- snp_attach(bed_rds)

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
sample.id  <- obj.bigSNP$fam



sample_to_include_file<-args[3]

sample_to_include<-getobj(sample_to_include_file)

ind.sample<- which(sample.id$sample.ID %in% sample_to_include)

map <- obj.bigSNP$map[-3]
head(map)
names(map) <- c("chr", "rsid", "pos", "a0", "a1")

#corr_dir<-"/data/sofer/Projects/2020_RNA_seq_SOL/CorMat"

corr_dir<-args[4]

tmpdir<-args[5]
tmp <- tempfile(tmpdir = tmpdir)


ncores<-as.integer(args[6])

for (chr in 1:22) {
  
  print(chr)
  # perform SNP matching
  tmp_snp<- snp_match(summary_stat_df[summary_stat_df$chr==chr,], map)

  
  corr0 <- readRDS(paste0(corr_dir,"/corr", chr, ".rds"))
  
  dim(corr0)
  
  print(dim(corr0))
  if (chr == 1) {
    df_beta <- tmp_snp[, c("rsid", "chr", "pos", "a0","a1","beta", "beta_se", "n_eff", "_NUM_ID_")]
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    df_beta <- rbind(df_beta, tmp_snp[, c("rsid", "chr", "pos", "a0","a1","beta", "beta_se", "n_eff", "_NUM_ID_")])
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}


head(df_beta)

(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
h2_est <- ldsc[["h2"]]


message("Compute weight auto ")


multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, vec_p_init = seq_log(1e-4, 0.9, length.out = 4), ncores = ncores)

file.remove(paste0(tmp, ".sbk"))

unlink(corr_dir, recursive = TRUE) 


beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)


message("Completed beta")

G2 <- snp_fastImputeSimple(G, ncores=ncores)

message("genotyped imputed done")

pred_auto <- big_prodMat(G2, beta_auto,ind.col = df_beta[["_NUM_ID_"]], ind.row = ind.sample)


message("done pred auto")

sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])

beta_out <- cbind(df_beta, final_beta_auto)


output_file<-args[7]


write.table(beta_out, file=output_file, row.names = F, quote = F, sep="\t")


