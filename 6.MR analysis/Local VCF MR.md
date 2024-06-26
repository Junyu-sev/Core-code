# Use Local VCF Documents to Perform Mendelian Randomization

## Download LD Reference datasets

```linux
wget http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz
tar xzvf data_maf0.01_rs_ref.tgz
```

## Local VCF MR Codes
Below example codes are originated from **[website](https://mrcieu.github.io/gwasglue/articles/mr.html#other-options-1)** here ~ Some require packages need to be installed from correspond git-hub website. Search these packages in Google and you will find how to install them. Additional `VariantAnnotation` package is needed for `gwasvcf` and you should install it first from `Bioconductor`.

> I have tried R4.2 & R4.3, and both R_versions could be used for codes below. Other R-versions may also work, not for sure.

```R
library(genetics.binaRies)
library(ieugwasr)
library(gwasglue)
library(gwasvcf)
library(tidyverse)
library(data.table)
library(TwoSampleMR)

gwasvcf::set_bcftools()
expd <- gwasvcf::query_gwas("Path/To/Exposure/data.vcf.gz OR Path/To/Exposure/data.vcf", pval=5e-8)
expd <- gwasglue::gwasvcf_to_TwoSampleMR(expd, type="exposure")
retain_snps <- expd %>% dplyr::select(rsid=SNP, pval=pval.exposure) %>%
    ieugwasr::ld_clump(., plink_bin=genetics.binaRies::get_plink_binary(), 
        bfile="/Path/To/data_maf0.01_rs_ref") %>% {.$rsid} # Don't add file extension to data_maf0.01_rs_ref
expd <- subset(expd, SNP %in% retain_snps)
gwasvcf::set_plink()
outd <- gwasvcf::query_gwas("Path/To/Outcome/data.vcf.gz OR Path/To/Outcome/data.vcf", 
    rsid = expd$SNP, proxies="yes", 
    bfile="/Path/To/data_maf0.01_rs_ref") # Don't add file extension to data_maf0.01_rs_ref
outd <- gwasglue::gwasvcf_to_TwoSampleMR(outd, "outcome")
dat <- TwoSampleMR::harmonise_data(expd, outd)
TwoSampleMR::mr(dat)
```

For FinnGen data, `sample size` is lack in GWAS data. You'd better add this into the extracted `expd` or `outd`, if you used FinnGen GWAS as exposure or outcome. For example of outcome, you can make this by:

```
outd$samplesize.outcome # Be careful. Check whether it is needed to fix sample size
outd$samplesize.outcome <- 412181 # Sample size of FinnGen is 412181
```

## Mendelian Randomization

Usefull introduction and tutorials could be found [here](https://mrcieu.github.io/TwoSampleMR/articles/gwas2020.html).