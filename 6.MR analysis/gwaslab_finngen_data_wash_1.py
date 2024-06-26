import pandas as pd
import numpy as np
import gwaslab as gl


Finngen = pd.read_csv("/share/home/linmiao/bio/rawdata/Finngen/summary_stats_R10_manifest.tsv", sep='\t')

gwas = Finngen['phenocode'].to_numpy()
gwas = gwas.astype(str)
gwas = np.core.defchararray.add('finngen_R10_', gwas)
gwas = np.sort(gwas)
gwas = gwas[0:241]

for gwas_id in gwas:
    mysumstats = gl.Sumstats("/share/home/linmiao/bio/rawdata/Finngen/phenotype/"+gwas_id+".gz", 
             snpid="rsids",
             chrom="#chrom",
             pos="pos",
             ea="alt",
             nea="ref",
             eaf="af_alt",
             beta="beta",
             se="sebeta",
             p="pval",
             mlog10p="mlogp",
             build="38")
    mysumstats.basic_check()
    mysumstats.assign_rsid( n_cores = 35,
                        ref_rsid_vcf ="/share/home/zhangjunyu/Rawdata/Reference_dbSNP/GCF_000001405.40.gz",
                        chr_dict = gl.get_number_to_NC(build="38"))
    mysumstats.liftover(n_cores=35, 
                    from_build="38", 
                    to_build="19",
                    remove=True)
    mysumstats.to_format("/share/home/linmiao/bio/rawdata/Finngen/phenotype_vcf/"+gwas_id, fmt="vcf", build="19",id_use="rsID",xymt_number=True)#可以再研究一下，以什么格式输出，后面在研究一下
