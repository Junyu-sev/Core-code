#需要对两份GWAS文件进行Harmonization
#需要计算Z-score，去除Z2>80的SNP
#

source("/share/home/zhangjunyu/Software/PLACO-master/PLACO_v0.1.1.R")

library(genetics.binaRies)
library(ieugwasr)
library(gwasglue)
library(gwasvcf)
library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(VariantAnnotation)
library(parallel)

gwasvcf::set_bcftools("/share/home/zhangjunyu/Software/bcftools/bin/bcftools")
gwasvcf::set_plink()

id1 <- "N14_RENFAIL"
id2 <- "ukb-d-I50"

#借助IEU的gwasvcf进行vcf文件读取，用MR包进行harmonise
GWAS_route1 <- paste("/share/home/linmiao/bio/rawdata/Finngen/phenotype_vcf/finngen_R10_", id1, ".vcf", sep="")
GWAS_route2 <- paste("/share/home/linmiao/bio/rawdata/MR/IEU_GWAS/", id2, ".vcf.gz", sep="")
GWAS1 <- gwasvcf::query_gwas(GWAS_route1, pval=1, threads = 35) #这里可能要调整一下
GWAS1 <- gwasglue::gwasvcf_to_TwoSampleMR(GWAS1, type="exposure")
GWAS1 <- GWAS1 %>% filter(eaf.exposure > 0.01)
GWAS2 <- gwasvcf::query_gwas(GWAS_route2, rsid = GWAS1$SNP, proxies="no", threads = 35)
GWAS2 <- gwasglue::gwasvcf_to_TwoSampleMR(GWAS2, "outcome")
dat <- TwoSampleMR::harmonise_data(GWAS1, GWAS2)

#计算Z统计量，制作Z.matrix和P.matrix
dat <- dat[,c("SNP", "beta.exposure", "se.exposure", "pval.exposure", "beta.outcome", "se.outcome", "pval.outcome")]
dat <- dat %>% mutate(Z1 = beta.exposure/se.exposure)
dat <- dat %>% mutate(Z2 = beta.outcome/se.outcome)
dat <- subset(dat, Z1 > -sqrt(80) & Z1 < sqrt(80))
dat <- subset(dat, Z2 > -sqrt(80) & Z2 < sqrt(80) )
p <- nrow(dat)
Z.matrix = dat[,c("Z1", "Z2")] %>% as.matrix()
P.matrix = dat[,c("pval.exposure", "pval.outcome")] %>% as.matrix()
colnames(P.matrix) <- c("P1", "P2")

#测试相关性，应该问题不大，我们的样本没有重叠
#R <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)
#	"%^%" <- function(x, pow)
		with(eigen(x), vectors %*% (values^pow * t(vectors)))
#Z.matrix.decor <- Z.matrix %*% (R %^% (-0.5))
#colnames(Z.matrix.decor) <- paste("Z",1:k,sep="")

VarZ <- var.placo(Z.matrix, P.matrix, p.threshold=1e-4)

cl <- makeCluster(35, outfile = "/share/home/zhangjunyu/Project/240925_HF/Code/PLACO/logfile.txt")
clusterEvalQ(cl, {
    source("/share/home/zhangjunyu/Software/PLACO-master/PLACO_v0.1.1.R")
})
clusterExport(cl, c("p", "Z.matrix", "VarZ"))
out_placo = parLapply(cl, 1:p, function(i) placo(Z=Z.matrix[i,], VarZ=VarZ))
stopCluster(cl)

out_placo <- out_placo %>% unlist
placo_p_value <- out_placo[seq(2, length(out_placo), 2)]

dat <- data.frame(dat, placo_p = placo_p_value)
fwrite(dat, paste("/share/home/zhangjunyu/Project/240925_HF/Result/PLACO/placo_finngen_R10_", id1, "AND", id2, ".csv", sep = ""), row.names = F)

dat_p_sig <- dat %>% filter(placo_p < 5E-8)
if(nrow(dat_p_sig) != 0){
    fwrite(dat_p_sig, paste("/share/home/zhangjunyu/Project/240925_HF/Result/PLACO/placo_sig_finngen_R10_", id1, "&", id2, ".csv", sep = ""), row.names = F)
}

#save.image(file = "/share/home/zhangjunyu/Project/240925_HF/Result/PLACO/placo_finngen_R10_N14_RENFAIL&ukb-d-HEARTFAIL.RData")


