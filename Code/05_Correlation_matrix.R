args<-commandArgs(trailingOnly = TRUE)

library(bigsnpr)
library(GWASTools)
library(dplyr)
library(data.table)
library(bigsnpr)


chr<-as.character(args[1])

message(paste0("run chr ", chr))

bed_rds <- paste(bed_files,sep="")

obj.bigSNP <- snp_attach(bed_rds)

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
sample.id  <- obj.bigSNP$fam



sample_to_include_file<-args[2]
sample_to_include<- getobj(sample_to_include_file)


ind.sample<- which(sample.id$sample.ID %in% sample_to_include)

message(paste0("run correlation matrix for ",length(ind.sample)), " people")


summary_stat_files<- args[3]

sumstats<-  bigreadr::fread2(summary_stat_files, data.table = F)
colnames(sumstats)<- tolower(colnames(sumstats))
head(sumstats)
# snp chr      pos a1 a2
colnames(sumstats)<-c("rsid", "chr", "pos", "a0","a1","beta")


map <- obj.bigSNP$map[-3]
head(map)
names(map) <- c("chr", "rsid", "pos", "a0", "a1")

# perform SNP matching

tmp_snp <- snp_match(sumstats[sumstats$chr==chr,], map)


tempdir<-args[4]
POS2 <- snp_asGeneticPos(CHR, POS, dir =tempdir )
# calculate LD
# Extract SNPs that are included in the chromosome
ind.chr <- which(tmp_snp$chr == chr)
ind.chr2 <- tmp_snp$`_NUM_ID_`[ind.chr]


message("Run correlation matrix")

ncores<-as.integer(config["ncores"])
corr0 <- snp_cor(
  G,
  ind.col = ind.chr2,
  ind.row = ind.sample,
  ncores = ncores,
  infos.pos = POS2[ind.chr2],
  size = 3 / 1000)

message(paste0("Correlation completed, the data will save to ",output_file ))

output_file<-args[5]

saveRDS(corr0, file = paste(output_file ,sep=""))

