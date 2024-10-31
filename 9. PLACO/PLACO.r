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

exposure_route <- "/share/home/linmiao/bio/rawdata/Finngen/phenotype_vcf/finngen_R10_N14_RENFAIL.vcf"
expd <- gwasvcf::query_gwas(exposure_route, pval=1, threads = 35) #这里可能要调整一下
expd <- gwasglue::gwasvcf_to_TwoSampleMR(expd, type="exposure")
expd <- expd %>% filter(eaf.exposure > 0.01)
outcome_route <- "/share/home/linmiao/bio/rawdata/MR/IEU_GWAS/ukb-d-HEARTFAIL.vcf.gz"

outd <- gwasvcf::query_gwas(outcome_route, rsid = expd$SNP, proxies="no", threads = 35)
outd <- gwasglue::gwasvcf_to_TwoSampleMR(outd, "outcome")
dat <- TwoSampleMR::harmonise_data(expd, outd)

dat <- dat[,c("SNP", "beta.exposure", "se.exposure", "pval.exposure", "beta.outcome", "se.outcome", "pval.outcome")]
dat <- dat %>% mutate(Z1 = beta.exposure/se.exposure)
dat <- dat %>% mutate(Z2 = beta.outcome/se.outcome)
dat <- subset(dat, Z1 > -sqrt(80) & Z1 < sqrt(80))
dat <- subset(dat, Z2 > -sqrt(80) & Z2 < sqrt(80) )
p <- nrow(dat)
Z.matrix = dat[,c("Z1", "Z2")] %>% as.matrix()
P.matrix = dat[,c("pval.exposure", "pval.outcome")] %>% as.matrix()
colnames(P.matrix) <- c("P1", "P2")

#R <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)

VarZ <- var.placo(Z.matrix, P.matrix, p.threshold=1e-4)
out <- sapply(1:p, function(i) placo(Z=Z.matrix[i,], VarZ=VarZ))
placo_p_value <- out[2,] %>% unlist
dat <- data.frame(dat, placo_p = placo_p_value)
fwrite(dat, "/share/home/zhangjunyu/Project/240925_HF/Result/PLACO/placo_finngen_R10_N14_RENFAIL&ukb-d-HEARTFAIL.csv", row.names = F)

save.image(file = "/share/home/zhangjunyu/Project/240925_HF/Result/PLACO/placo_finngen_R10_N14_RENFAIL&ukb-d-HEARTFAIL.RData")

dat_p_sig <- dat %>% filter(placo_p < 5E-8)
if(nrow(dat_p_sig) != 0){
    fwrite(dat_p_sig, "/share/home/zhangjunyu/Project/240925_HF/Result/PLACO/placo_sig_finngen_R10_N14_RENFAIL&ukb-d-HEARTFAIL.csv", row.names = F)
}


#out1 <- mclapply(1:p, function(i) placo(Z=Z.matrix[i,], VarZ=VarZ), mc.cores=35)
#Error in mcfork() : unable to fork, possible reason: Cannot allocate memory
#Warning message:
#In mclapply(1:p, function(i) placo(Z = Z.matrix[i, ], VarZ = VarZ),  :
#  scheduled cores 4, 5, 7, 12, 32 did not deliver results, all values of the jobs will be affected
#总之，这样的并行计算不稳定，和节点有关嘛？我在fat01节点跑的，emm...





