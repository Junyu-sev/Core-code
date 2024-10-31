#!/bin/bash
#BSUB -J GWAS
#BSUB -n 40
#BSUB -R "span[ptile=40]"
#BSUB -o output_%J
#BSUB -e output_%J
#BSUB -q mpi

workdir="/share/home/zhangjunyu/Project/240925_HF/Result/GWAS/Renal_injury_before_HF"
Inclusion="/share/home/zhangjunyu/Project/240925_HF/Result/GWAS/Renal_injury_before_HF/Inclusion.txt"
Pheno="/share/home/zhangjunyu/Project/240925_HF/Result/GWAS/Renal_injury_before_HF/Pheno.txt"
GWAS_cov="/share/home/zhangjunyu/Project/240925_HF/Result/GWAS/Renal_injury_before_HF/Cov.txt"

cd ${workdir}

Filefold2="02.GWAS_QC"
mkdir -p ${Filefold2}
/share/home/zhangjunyu/Software/plink/plink2/plink2 \
    --memory 2000000 --threads 35 \
    --pfile ./01.Extract_participants/merge \
    --mac 10\
    --maf 0.0001 \
    --geno 0.1 \
    --mind 0.1 \
    --hwe 1e-15 \
    --make-pgen \
    --out ./${Filefold2}/GWAS_QC

Filefold3="03.GWAS"
mkdir -p ${Filefold3}
/share/home/zhangjunyu/Software/plink/plink2/plink2 \
    --memory 500000 --threads 40 \
    --pfile ./02.GWAS_QC/GWAS_QC \
    --pheno ${Pheno} \
    --pheno-col-nums 3 \
    --covar ${GWAS_cov} \
    --covar-col-nums 3-14 \
    --no-input-missing-phenotype \
    --covar-variance-standardize \
    --logistic \
    --out ./${Filefold3}/GWAS
