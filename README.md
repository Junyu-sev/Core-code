# Core-code
- Core codes when studying in Guangdong Provincial People's Hospital, Guangzhou.
- Personal use & ready to share.
- Most codes are adjusted from tutorials of related packages/softwares/articles.
- Reference of each package/software/method is omitted.

# Contents

## 0.Notes (in Chinese)
- Notes during studying and code writing.

## 1.Disease definition for UKB
- **R**
- ICD10 & OPSC4
- Disease contains IHD, HF...

## 2.MICE to impute missings in covariates
- **R**
- 5 times of imputation
- Rubin rules

## 3.Cox model to get HR
- **R/Python**
- In most cases, plasma proteins are risk factors
- Covariates adjustment

## 4.Machine learning model for prediction
- **Python**
- kmeans, LGBM, RF, SVM, XGBoost, DT, LR...
- SMOTE/Downsampling to deal with imbalance problems
- Bayesian-optimization
- Accuracy, precision, recall...
- ROC AUC curve & PR curve
- Cross-validation (k- fold, External validation...)

## 5.KM analysis
- **Python**
- Youden Index to determine cutoffs for each protein (biomarker effect)
- Numbers at risk & numbers at event

## 6.MR analysis
- **Python/R**
- Gwaslab to wash data into vcf framat in GR37/hg19
- Series of IEU packages to conduct Two-sample MR analysis locally

## 7.Colocalization analysis
- **R**
- Gwasglue to match vcf files in given chrompos
- Conc for colocalization analysis

## 8.GWAS & PRS
- **sh/R**
- Plink2 to perform GWAS analysis
- PRSice-2, LDpred2, and lassosum for PRS calculation

## 9.PLACO
- **R**
- Pleiotropic effect of a genetic variant on two traits
- Reference: https://github.com/RayDebashree/PLACO
