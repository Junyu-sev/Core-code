# Calculate ROC AUC confidence interval 
#https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_the_Area_Under_an_ROC_Curve.pdf
#https://cs.nyu.edu/~mohri/pub/area.pdf
# N1 The number of subjects sampled from the "positive" group.
# N2 The number of subjects sampled from the "negative" group
AUC = 0.76
N1 = 2547+3699
N2 = 9068
Q1 = AUC/(2-AUC)
Q2 = 2*AUC*AUC/(1+AUC)
se = ((AUC*(1-AUC)+(N1-1)*(Q1-AUC*AUC)+(N2-1)*(Q2-AUC*AUC))/(N1*N2))**0.5
AUC + se*1.96
AUC - se*1.96
#但是R的包是通过bootstrap来计算置信区间的

#Meta analysis through R
#https://methods.cochrane.org/sites/methods.cochrane.org.prognosis/files/uploads/r_script_meta-analysis_c-statistic.txt