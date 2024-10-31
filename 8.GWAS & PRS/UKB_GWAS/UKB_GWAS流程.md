# UKB GWAS

## 概述

- 有Genotype数据的UKB人群有47w+人，**如果进行全部流程的GWAS，会花去大量的时间**。即使设定40核运行，计算kinship需要1-2周，计算PCA还需要1-2周，不可能花这么多时间的。
- UKB RAP上有做GWAS的教程，运行起来很快，大概1-2天能跑完。类似于紫英师姐一开始教我的过程，不过多了很多人群纳排的环节。结合我看到的UKB GWAS的文献，**我感觉UKB RAP上面的这个过程可以参考，我们在HPC上面实现**。
- 相比于全流程的GWAS，我们新的流程直接借用UKB基因型分析的各项结果（包括F系数、PCA等等），**省去了一些GWAS QC环节和全部的PCA过程**，因而有更短的运行时间。但是我目前感觉**这些省略应该也是没问题的**。

## 参考资料

1. [RAP上UKB做GWAS的教程1](https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/end-to-end-target-discovery-with-gwas-and-phewas)：这是一份基于RAP平台的教程，利用RAP上面的JupyterLab和Swiss Army Knife等完成GWAS的。我们参考这个RAP教程在HPC上做GWAS。
2. [RAP上UKB做GWAS的教程2](https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/gwas-using-alzheimers-disease)：这也是一份基于RAP平台的教程，人群纳排和上面这个有些许不同。

## 操作过程

### 数据准备

1. **待提取个体文件**
    **在你自己做的人群纳排的基础上，你需要额外进行以下的人群纳排**：
      - 报告性别和遗传性别是相同的
      - 基因判别为英国白人血统
      - 排除性染色体非整倍体
      - 参与者用于计算 PCA，其实也就是去掉有亲属关系的个体。
      - 去掉杂合子和缺失率的异常值

    我在这个文件`/share/home/zhangjunyu/Rawdata/UKB_GWAS/UKB_GWAS_id_QC.txt`存放了符合这五种要求的eid清单，你可以对你的数据进行过滤，只**保留这个清单内的eid编号**即可。

    最后将进行GWAS的个体编号总结在一份.txt文件中，需要两列eid编号，并将列名改为"FID"和"IID"，格式如下：

    ```sh
    FID	IID
    1000619	1000619
    1003210	1003210
    1004545	1004545
    1004716	1004716
    1005251	1005251
    1007707	1007707
    1008517	1008517
    1009506	1009506
    1009521	1009521
    ```

    将这份txt文件上传到HPC中，并记录好它所在的路径。

2. **表型文件**
   整理表型数据到一份.txt文件中。表型文件固定前两列为个体编号，列名为"FID"和"IID"，第三列或其他列为表型。
   - 连续变量作为表型，格式如下：

    ```sh
    FID	IID	pheno
    1000619	1000619	56.1603
    1003210	1003210	76.1989
    1004545	1004545	69.6949
    1004716	1004716	66.5986
    1005251	1005251	58.7845
    1007707	1007707	55.6963
    1008517	1008517	58.7492
    1009506	1009506	23.5583
    1009521	1009521	66.5572
    ```

   - 二分类变量作为表型，**注意：二分类性状中1表示control，2表示case**，格式如下：

    ```sh
    FID	IID	pheno
    1000038	1000038	1
    1000093	1000093	1
    1000108	1000108	2
    1000115	1000115	1
    1000151	1000151	1
    1000166	1000166	2
    1000173	1000173	1
    1000182	1000182	1
    1000201	1000201	1
    ```

    - 多分类变量利用plink实现似乎有些困难，可以参考[其他软件/包](https://blog.sciencenet.cn/blog-3423233-1281079.html)来实现。

    将这份txt文件上传到HPC中，并记录好它所在的路径。

3. **协变量文件**
   我们选择age、sex和PCA的前十个主成分作为协变量，下面代码是生成协变量文件的参考代码，你可以根据实际情况增添或删减协变量。

   ```R
   library(tidyverse)
   library(data.table)

   #Inclusion文件路径，要根据实际情况替换文件位置
   Inclusion <- fread("/share/home/zhangjunyu/Project/240925_HF/Result/GWAS/No_Renal_injury/Inclusion.txt")
   covariate_age <- fread("/share/home/zhangjunyu/Rawdata/UKB/Covariates/covariate.csv")
   covariate_age <- covariate_age %>% select(eid, age)
   covariate_genetic <- fread("/share/home/zhangjunyu/Rawdata/UKB_GWAS/GWAS_cov.csv")
   covariate_genetic <- covariate_genetic[, c(1, 5, 13:22)]
   covariate_genetic <- covariate_genetic %>% mutate(Genetic_sex = ifelse(Genetic_sex == "Female", 2, 1))#GWAS编码的性别中，1代表男性，2代表女性
   covariate <- inner_join(covariate_age, covariate_genetic, by = "eid")
   covariate <- covariate %>% rename(IID="eid")
   covariate_merge <- left_join(Inclusion, covariate, by = "IID")
   
   #输出文件，要根据实际情况替换输出位置，记好输出文件的路径
   fwrite(covariate_merge, "/share/home/zhangjunyu/Project/240925_HF/Result/GWAS/No_Renal_injury/Cov.txt", sep = "\t", quote = F, row.names = F)
   ```

### 进行GWAS

创建工作路径，进入工作路径，将附件中的`GWAS_binary_variable.sh`和`GWAS_continuous_variable.sh`文件，**按照实际情况修改后**上传，提交任务。

1. **连续变量的GWAS**
   连续变量作为表型的GWAS，比如身高、血糖等等。

    ```sh
    #工作路径（你准备进行GWAS分析的路径）
    workdir="/share/home/zhangjunyu/Project/240925_HF/Result/GWAS/No_Renal_injury"
    cd ${workdir}
    #提交任务
    bsub < GWAS_continuous_variable.sh
    ```

    如果你需要增多/删减协变量，或者多个表型批次分析，请在代码中修改相应的部分~

2. **二分类变量的GWAS**
   二分类变量作为表型的GWAS，比如疾病发生、是否有某种习惯等等。

    ```sh
    #工作路径（你准备进行GWAS分析的路径）
    workdir="/share/home/zhangjunyu/Project/240925_HF/Result/GWAS/No_Renal_injury"
    cd ${workdir}
    #提交任务
    bsub < GWAS_binary_variable.sh
    ```

    如果你需要增多/删减协变量，或者多个表型批次分析，请在代码中修改相应的部分~

3. **提取最终的GWAS结果**

   - 对于连续变量，找到".glm.linear"后缀的文件（下面的文件名需要替换为自己的文件），采用如下代码：

    ```sh
    head -n 1 ./03.GWAS/GWAS.pheno.glm.linear > GWAS.txt
    grep "ADD" ./03.GWAS/GWAS.pheno.glm.linear >> GWAS.txt
    ```

   - 对于二分类变量，找到".glm.logistic.hybrid"后缀的文件（下面的文件名需要替换为自己的文件），采用如下代码：

    ```sh
    head -n 1 ./03.GWAS/GWAS.pheno.glm.logistic.hybrid > GWAS.txt
    grep "ADD" ./03.GWAS/GWAS.pheno.glm.logistic.hybrid >> GWAS.txt
    ```

    **注意：UKB最终的GWAS数据，经过核实，是GRCh37（hg19）版本的数据。**