#这个代码是依据师姐给的参考代码，整理出来的代码合集
#随机拆分数据集（8：2），分出来训练集和测试集
############################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
from sklearn.ensemble import RandomForestClassifier 
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import export_graphviz 
from sklearn import metrics
from sklearn.metrics import roc_curve, auc 
from sklearn.metrics import classification_report 
from sklearn.metrics import confusion_matrix 
from sklearn.model_selection import train_test_split 
import shap 
from pdpbox import pdp, info_plots 
# 分别训练 逻辑回归（Logistic Regression）、决策树（Decision Tree）、随机森林（Random Forest）、XGBoost、LightGBM 
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
import lightgbm as lgb
#画图
from sklearn.metrics import ConfusionMatrixDisplay
import matplotlib.pyplot as plt


np.random.seed(2022) 
pd.options.mode.chained_assignment = None  

dpath = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Data/'
outpath = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/'
pltpath = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/'

pro_f_lst = [
    "ACVRL1", "ADA2", "AREG", "B4GALT1", "BCL2L11", "BSG", "BST2", "BTN3A2",
    "CCL18", "CD14", "CD4", "CD74", "COLEC12", "CST3", "CTSL", "CTSS", "CTSZ",
    "CXCL16", "DDAH1", "EGFR", "EPHA2", "ESAM", "FABP4", "FAM3C", "FAS", "GDF15",
    "GFRA1", "HAVCR2", "HGF", "HLA-E", "HSPG2", "IGFBP7", "IGSF8", "ITGA5", "JAM2",
    "LAMA4", "LGALS1", "LGALS3", "LGMN", "LTBP2", "LTBR", "MB", "MSR1",
    "NDUFS6", "NECTIN2", "OGN", "PIK3IP1", "PLA2G2A", "PLA2G7", "PLAUR", "POLR2F",
    "PTGDS", "RELT", "RNASET2", "SCARB2", "SPP1", "TNFRSF10B", "TNFRSF11A", "TNFRSF14",
    "TNFRSF1B", "TREM2", "TXNDC15", "VEGFA", "WNT9A", "XG"
]

pro_df = pd.read_csv(dpath + 'Raw_data/231008_data_olink_instance_0.csv',usecols = ['ParticipantID1'] + pro_f_lst)
target_df = pd.read_csv(dpath + 'Disease_outcomes/IHD_control_definition/4without_nonIHD/IHD_outcomes.csv', usecols = ['eid', 'target_y', 'BL2Target_yrs'])
mydf = pd.merge(target_df, pro_df, how = 'inner', left_on = ['eid'], right_on = ['ParticipantID1'])
mydf = mydf.drop('ParticipantID1', axis=1)
mydf = mydf.drop('BL2Target_yrs', axis=1)

#reg_df = pd.read_csv(dpath + 'Eid_info_data.csv', usecols = ['eid', 'Region_code'])
#mydf = pd.merge(mydf, reg_df, how = 'left', on = ['eid'])

# 检查每个值是否为空值
null_values = mydf.isnull()
# 统计每列中的空值数量
null_counts = mydf.isnull().sum()

#去掉空值
mydf = mydf.dropna()######这里对于缺失值的处理有待商榷，需要再想想再想想

# 2.1.1 返回频数、均值、标准差、四分位数
print(mydf.describe(include="all"))

## 2.1.2计算数据的偏度 
#偏度值为0，则表示数据分布对称。如果偏度值为正数，则表示数据分布向左倾斜，如果偏度值为负数，则表示数据分布向右倾斜
print(mydf.skew())

# 3.1 变量间相关性
plt.figure(figsize=(20,18))
#采用sns.heatmap()是Seaborn库中的一个函数，用于绘制热力图
sns.heatmap(mydf.corr(),annot=True,cmap='YlGnBu',fmt='.2f',linewidths=2)
plt.savefig(pltpath + 'heatmapplot.pdf')

# 3.1 变量间相关性,画不了，因为heatmap这个包是py2的，环境是3.9的
#from heatmap import heatmap, corrplot
#plt.figure(figsize=(25, 25))
#corrplot(dc.corr(), size_scale=300)
#plt.savefig('./discovery/0511_corrplot.pdf')

X = mydf.drop(['eid','target_y'], axis=1)
Y = mydf['target_y']

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.20, random_state=1234,shuffle=True, stratify=None)

#很明显，我们的数据存在类别不平衡问题！！！暂时设定的class_weight='balanced'参数，有点用，但不多
#开始构建模型
# 逻辑回归
lr_model = LogisticRegression(C=1.0,class_weight='balanced')
lr_model.fit(X_train, Y_train)
# 决策树
dt_model = DecisionTreeClassifier(max_depth=None, min_samples_split=2, min_samples_leaf=1,class_weight='balanced')
dt_model.fit(X_train, Y_train)
# 随机森林
rf_model = RandomForestClassifier(n_estimators=100, max_depth=None, min_samples_split=2, min_samples_leaf=1,class_weight='balanced')
rf_model.fit(X_train, Y_train)
# XGBoost
xgb_model = xgb.XGBClassifier(max_depth=6, learning_rate=0.3, n_estimators=100, subsample=1, colsample_bytree=1,eval_metric='error',class_weight='balanced')
xgb_model.fit(X_train, Y_train)
# LightGBM
lgb_model = lgb.LGBMClassifier(max_depth=-1, learning_rate=0.1, n_estimators=100, subsample=1, colsample_bytree=1,class_weight='balanced')
lgb_model.fit(X_train, Y_train)

#my_params = {'n_estimators': 500,
#            'max_depth': 15,
#             'num_leaves': 10,
#             'subsample': 0.7,
#             'learning_rate': 0.01,
#             'colsample_bytree': 0.7}

#lgb_model = lgb.LGBMClassifier(max_depth=-1, learning_rate=0.1, n_estimators=100, subsample=1, colsample_bytree=1,class_weight='balanced')

#lgb_model = lgb.LGBMClassifier(objective = 'binary', metric = 'auc', is_unbalance = True, verbosity = 1, seed = 2023)
#lgb_model.set_params(**my_params)
#lgb_model.fit(X_train, Y_train)

# 逻辑回归混淆矩阵
lr_y_pred_bin = lr_model.predict(X_test)
lr_cm = confusion_matrix(Y_test, lr_y_pred_bin)
# 绘制逻辑回归混淆矩阵
plt.figure(figsize=(8, 6))
sns.heatmap(lr_cm, annot=True, fmt="d", cmap="Blues", xticklabels=["Control", "IHD"], yticklabels=["Control", "IHD"])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Logistic Regression Confusion Matrix')
plt.savefig(pltpath + "0511_lr_confusion.jpg", dpi=450)

# 决策树混淆矩阵
dt_y_pred_bin = dt_model.predict(X_test)
dt_cm = confusion_matrix(Y_test, dt_y_pred_bin)
# 绘制决策树混淆矩阵
plt.figure(figsize=(8, 6))
sns.heatmap(dt_cm, annot=True, fmt="d", cmap="Blues", xticklabels=["Control", "IHD"], yticklabels=["Control", "IHD"])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Decision Tree Confusion Matrix')
plt.savefig(pltpath + "0511_dt_confusion.jpg", dpi=450)

# 随机森林混淆矩阵
rf_y_pred_bin = rf_model.predict(X_test)
rf_cm = confusion_matrix(Y_test, rf_y_pred_bin)
# 绘制随机森林混淆矩阵
plt.figure(figsize=(8, 6))
sns.heatmap(rf_cm, annot=True, fmt="d", cmap="Blues", xticklabels=["Control", "IHD"], yticklabels=["Control", "IHD"])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Random Forest Confusion Matrix')
plt.savefig(pltpath + "0511_rf_confusion.jpg", dpi=450)

# XGBoost混淆矩阵
xgb_y_pred_bin = xgb_model.predict(X_test)
xgb_cm = confusion_matrix(Y_test, xgb_y_pred_bin)
# 绘制XGBoost混淆矩阵
plt.figure(figsize=(8, 6))
sns.heatmap(xgb_cm, annot=True, fmt="d", cmap="Blues", xticklabels=["Control", "IHD"], yticklabels=["Control", "IHD"])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('XGBoost Confusion Matrix')
plt.savefig(pltpath + "0511_xgb_confusion.jpg", dpi=450)

# LightGBM混淆矩阵
lgb_y_pred_bin = lgb_model.predict(X_test)
lgb_cm = confusion_matrix(Y_test, lgb_y_pred_bin)
# 绘制LightGBM混淆矩阵
plt.figure(figsize=(8, 6))
sns.heatmap(lgb_cm, annot=True, fmt="d", cmap="Blues", xticklabels=["Control", "IHD"], yticklabels=["Control", "IHD"])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('LightGBM Confusion Matrix')
plt.savefig(pltpath + "0511_lgb_confusion.jpg", dpi=450)


from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

# 对测试集进行预测
lr_pred = lr_model.predict(X_test)
dt_pred = dt_model.predict(X_test)
rf_pred = rf_model.predict(X_test)
xgb_pred = xgb_model.predict(X_test)
lgb_pred = lgb_model.predict(X_test)

# 计算准确率
lr_accuracy = accuracy_score(Y_test, lr_pred)
dt_accuracy = accuracy_score(Y_test, dt_pred)
rf_accuracy = accuracy_score(Y_test, rf_pred)
xgb_accuracy = accuracy_score(Y_test, xgb_pred)
lgb_accuracy = accuracy_score(Y_test, lgb_pred)

# 计算精确率
lr_precision = precision_score(Y_test, lr_pred)
dt_precision = precision_score(Y_test, dt_pred)
rf_precision = precision_score(Y_test, rf_pred)
xgb_precision = precision_score(Y_test, xgb_pred)
lgb_precision = precision_score(Y_test, lgb_pred)

# 计算召回率
lr_recall = recall_score(Y_test, lr_pred)
dt_recall = recall_score(Y_test, dt_pred)
rf_recall = recall_score(Y_test, rf_pred)
xgb_recall = recall_score(Y_test, xgb_pred)
lgb_recall = recall_score(Y_test, lgb_pred)

# 计算 F1 分数
lr_f1_score = f1_score(Y_test, lr_pred)
dt_f1_score = f1_score(Y_test, dt_pred)
rf_f1_score = f1_score(Y_test, rf_pred)
xgb_f1_score = f1_score(Y_test, xgb_pred)
lgb_f1_score = f1_score(Y_test, lgb_pred)

# 计算 AUC-ROC
lr_auc_roc = roc_auc_score(Y_test, lr_pred)
dt_auc_roc = roc_auc_score(Y_test, dt_pred)
rf_auc_roc = roc_auc_score(Y_test, rf_pred)
xgb_auc_roc = roc_auc_score(Y_test, xgb_pred)
lgb_auc_roc = roc_auc_score(Y_test, lgb_pred)

# 打印评估结果
print("逻辑回归的准确率：", lr_accuracy)
print("决策树的准确率：", dt_accuracy)
print("随机森林的准确率：", rf_accuracy)
print("XGBoost 的准确率：", xgb_accuracy)
print("LightGBM 的准确率：", lgb_accuracy)
print("---")
print("逻辑回归的精确率：", lr_precision)
print("决策树的精确率：", dt_precision)
print("随机森林的精确率：", rf_precision)
print("XGBoost 的精确率：", xgb_precision)
print("LightGBM 的精确率：", lgb_precision)

# 计算评估指标的均值和标准差
lr_accuracy_mean, lr_accuracy_sd = np.mean(lr_accuracy), np.std(lr_accuracy)
dt_accuracy_mean, dt_accuracy_sd = np.mean(dt_accuracy), np.std(dt_accuracy)
rf_accuracy_mean, rf_accuracy_sd = np.mean(rf_accuracy), np.std(rf_accuracy)
xgb_accuracy_mean, xgb_accuracy_sd = np.mean(xgb_accuracy), np.std(xgb_accuracy)
lgb_accuracy_mean, lgb_accuracy_sd = np.mean(lgb_accuracy), np.std(lgb_accuracy)

lr_precision_mean, lr_precision_sd = np.mean(lr_precision), np.std(lr_precision)
dt_precision_mean, dt_precision_sd = np.mean(dt_precision), np.std(dt_precision)
rf_precision_mean, rf_precision_sd = np.mean(rf_precision), np.std(rf_precision)
xgb_precision_mean, xgb_precision_sd = np.mean(xgb_precision), np.std(xgb_precision)
lgb_precision_mean, lgb_precision_sd = np.mean(lgb_precision), np.std(lgb_precision)

lr_recall_mean, lr_recall_sd = np.mean(lr_recall), np.std(lr_recall)
dt_recall_mean, dt_recall_sd = np.mean(dt_recall), np.std(dt_recall)
rf_recall_mean, rf_recall_sd = np.mean(rf_recall), np.std(rf_recall)
xgb_recall_mean, xgb_recall_sd = np.mean(xgb_recall), np.std(xgb_recall)
lgb_recall_mean, lgb_recall_sd = np.mean(lgb_recall), np.std(lgb_recall)

lr_f1_score_mean, lr_f1_score_sd = np.mean(lr_f1_score), np.std(lr_f1_score)
dt_f1_score_mean, dt_f1_score_sd = np.mean(dt_f1_score), np.std(dt_f1_score)
rf_f1_score_mean, rf_f1_score_sd = np.mean(rf_f1_score), np.std(rf_f1_score)
xgb_f1_score_mean, xgb_f1_score_sd = np.mean(xgb_f1_score), np.std(xgb_f1_score)
lgb_f1_score_mean, lgb_f1_score_sd = np.mean(lgb_f1_score), np.std(lgb_f1_score)

lr_auc_roc_mean, lr_auc_roc_sd = np.mean(lr_auc_roc), np.std(lr_auc_roc)
dt_auc_roc_mean, dt_auc_roc_sd = np.mean(dt_auc_roc), np.std(dt_auc_roc)
rf_auc_roc_mean, rf_auc_roc_sd = np.mean(rf_auc_roc), np.std(rf_auc_roc)
xgb_auc_roc_mean, xgb_auc_roc_sd = np.mean(xgb_auc_roc), np.std(xgb_auc_roc)
lgb_auc_roc_mean, lgb_auc_roc_sd = np.mean(lgb_auc_roc), np.std(lgb_auc_roc)

# 打印评估结果
print("逻辑回归的准确率：{:.2f} ± {:.2f}".format(lr_accuracy_mean, lr_accuracy_sd))
print("决策树的准确率：{:.2f} ± {:.2f}".format(dt_accuracy_mean, dt_accuracy_sd))
print("随机森林的准确率：{:.2f} ± {:.2f}".format(rf_accuracy_mean, rf_accuracy_sd))
print("XGBoost 的准确率：{:.2f} ± {:.2f}".format(xgb_accuracy_mean, xgb_accuracy_sd))
print("LightGBM 的准确率：{:.2f} ± {:.2f}".format(lgb_accuracy_mean, lgb_accuracy_sd))
print("---")
print("逻辑回归的精确率：{:.2f} ± {:.2f}".format(lr_precision_mean, lr_precision_sd))
print("决策树的精确率：{:.2f} ± {:.2f}".format(dt_precision_mean, dt_precision_sd))
print("随机森林的精确率：{:.2f} ± {:.2f}".format(rf_precision_mean, rf_precision_sd))
print("XGBoost 的精确率：{:.2f} ± {:.2f}".format(xgb_precision_mean, xgb_precision_sd))
print("LightGBM 的精确率：{:.2f} ± {:.2f}".format(lgb_precision_mean, lgb_precision_sd))
print("---")
print("逻辑回归的召回率：{:.2f} ± {:.2f}".format(lr_recall_mean, lr_recall_sd))
print("决策树的召回率：{:.2f} ± {:.2f}".format(dt_recall_mean, dt_recall_sd))
print("随机森林的召回率：{:.2f} ± {:.2f}".format(rf_recall_mean, rf_recall_sd))
print("XGBoost 的召回率：{:.2f} ± {:.2f}".format(xgb_recall_mean, xgb_recall_sd))
print("LightGBM 的召回率：{:.2f} ± {:.2f}".format(lgb_recall_mean, lgb_recall_sd))
print("---")
print("逻辑回归的F1 分数：{:.2f} ± {:.2f}".format(lr_f1_score_mean, lr_f1_score_sd))
print("决策树的F1 分数：{:.2f} ± {:.2f}".format(dt_f1_score_mean, dt_f1_score_sd))
print("随机森林的F1 分数：{:.2f} ± {:.2f}".format(rf_f1_score_mean, rf_f1_score_sd))
print("XGBoost 的F1 分数：{:.2f} ± {:.2f}".format(xgb_f1_score_mean, xgb_f1_score_sd))
print("LightGBM 的F1 分数：{:.2f} ± {:.2f}".format(lgb_f1_score_mean, lgb_f1_score_sd))
print("---")
print("逻辑回归的AUC-ROC：{:.2f} ± {:.2f}".format(lr_auc_roc_mean, lr_auc_roc_sd))
print("决策树的AUC-ROC：{:.2f} ± {:.2f}".format(dt_auc_roc_mean, dt_auc_roc_sd))
print("随机森林的AUC-ROC：{:.2f} ± {:.2f}".format(rf_auc_roc_mean, rf_auc_roc_sd))
print("XGBoost 的AUC-ROC：{:.2f} ± {:.2f}".format(xgb_auc_roc_mean, xgb_auc_roc_sd))
print("LightGBM 的AUC-ROC：{:.2f} ± {:.2f}".format(lgb_auc_roc_mean, lgb_auc_roc_sd))

# Create DataFrame with evaluation results
data = {
    'Model': ['Logistic Regression', 'Decision Tree', 'Random Forest', 'XGBoost', 'LightGBM'],
    'Accuracy': ["{:.2f} ± {:.2f}".format(lr_accuracy_mean, lr_accuracy_sd),
                 "{:.2f} ± {:.2f}".format(dt_accuracy_mean, dt_accuracy_sd),
                 "{:.2f} ± {:.2f}".format(rf_accuracy_mean, rf_accuracy_sd),
                 "{:.2f} ± {:.2f}".format(xgb_accuracy_mean, xgb_accuracy_sd),
                 "{:.2f} ± {:.2f}".format(lgb_accuracy_mean, lgb_accuracy_sd)],
    'Precision': ["{:.2f} ± {:.2f}".format(lr_precision_mean, lr_precision_sd),
                  "{:.2f} ± {:.2f}".format(dt_precision_mean, dt_precision_sd),
                  "{:.2f} ± {:.2f}".format(rf_precision_mean, rf_precision_sd),
                  "{:.2f} ± {:.2f}".format(xgb_precision_mean, xgb_precision_sd),
                  "{:.2f} ± {:.2f}".format(lgb_precision_mean, lgb_precision_sd)],
    'Recall': ["{:.2f} ± {:.2f}".format(lr_recall_mean, lr_recall_sd),
               "{:.2f} ± {:.2f}".format(dt_recall_mean, dt_recall_sd),
               "{:.2f} ± {:.2f}".format(rf_recall_mean, rf_recall_sd),
               "{:.2f} ± {:.2f}".format(xgb_recall_mean, xgb_recall_sd),
               "{:.2f} ± {:.2f}".format(lgb_recall_mean, lgb_recall_sd)],
    'F1 Score': ["{:.2f} ± {:.2f}".format(lr_f1_score_mean, lr_f1_score_sd),
                 "{:.2f} ± {:.2f}".format(dt_f1_score_mean, dt_f1_score_sd),
                 "{:.2f} ± {:.2f}".format(rf_f1_score_mean, rf_f1_score_sd),
                 "{:.2f} ± {:.2f}".format(xgb_f1_score_mean, xgb_f1_score_sd),
                 "{:.2f} ± {:.2f}".format(lgb_f1_score_mean, lgb_f1_score_sd)],
    'AUC-ROC': ["{:.2f} ± {:.2f}".format(lr_auc_roc_mean, lr_auc_roc_sd),
                "{:.2f} ± {:.2f}".format(dt_auc_roc_mean, dt_auc_roc_sd),
                "{:.2f} ± {:.2f}".format(rf_auc_roc_mean, rf_auc_roc_sd),
                "{:.2f} ± {:.2f}".format(xgb_auc_roc_mean, xgb_auc_roc_sd),
                "{:.2f} ± {:.2f}".format(lgb_auc_roc_mean, lgb_auc_roc_sd)]
}

df = pd.DataFrame(data)

# Export DataFrame as CSV
df.to_csv(pltpath + '0511_evaluation_results.csv', index=False)

#shap选择关键基因/蛋白
import shap 
from pdpbox import pdp, info_plots 
np.random.seed(2022) 

pd.options.mode.chained_assignment = None  
from sklearn.inspection import permutation_importance

# 逻辑回归模型
# 计算每个变量对模型的重要性
perm = permutation_importance(lr_model, X_test, Y_test, random_state=10)
feature_importances = perm.importances_mean
# 打印特征排序均值
print(feature_importances)

# 可视化前10个重要特征
import matplotlib.pyplot as plt
import seaborn as sns
features = X.columns
data_tuples = list(zip(feature_importances, features))
data = pd.DataFrame(data_tuples, columns=['importances', 'features'])
data = data.sort_values('importances', ascending=False)[:20]

plt.figure(figsize=(20, 6))
ax = sns.barplot(x=data['importances'], y=data['features'], palette="Blues_r", orient='h')
plt.savefig(pltpath + "0512_barplot.pdf")
data.to_csv(pltpath + '0512_dataimportance.csv')

import shap
from sklearn.inspection import permutation_importance
import xgboost as xgb
import lightgbm as lgb

# summary_plot
# summarize the effects of all the features
# SHAP
import shap
# 创建SHAP解释器
explainer = shap.KernelExplainer(lr_model.predict_proba, X_train)
# 计算SHAP值
shap_values = explainer.shap_values(X_test)#为什么这步那么那么慢，好奇怪啊
# 显示SHAP值
print(shap_values)
#shap可视化
shap.summary_plot(shap_values[1],X_test,plot_type="bar")
plt.savefig(pltpath + '0512_shap_values1.png')

# 全局条形图
# 特征重要性的条形图还有另一种绘制方法。
# 输出shap.Explanation对象
# compute SHAP values
# 创建SHAP解释器
explainer = shap.KernelExplainer(lr_model.predict_proba, X_train)
# 计算SHAP值
shap_values2 = explainer.shap_values(X_test, check_additivity=False)
# "check_additivity" 参数用于控制是否检查加性约束。加性约束表示每个样本的所有特征的 SHAP 值之和应该等于该样本的模型输出。默认情况下，SHAP 库会检查这个加性约束是否成立，如果不成立，可能会提示警告。
# 将 "check_additivity" 参数设置为 False 可以跳过这种检查，这样可以提高计算效率，特别是当对大型数据集进行 SHAP 值计算时。
# 但是完全没有变快，真奇怪啊，好吧，快了一捏捏
# 显示SHAP值
figure = plt.figure(figsize=(16,10))
shap.summary_plot(shap_values2, X_test, plot_type='bar', max_display=20, show=False)
plt.savefig(pltpath + "0512_shap_values2.jpg")

shap.summary_plot(shap_values[1],X_test)
# 转换为DataFrame并选择前20行进行打印
# 选择特定的目标类别
target_class = 1
# 转换为DataFrame并选择前20行进行打印
shap_df = pd.DataFrame(shap_values[target_class], columns=X_test.columns)
top20_shap_df = shap_df.head(20)
top20_shap_df.to_csv(pltpath + 'top20_shap_values.csv', index=False)

shap.summary_plot(shap_values2[1],X_test)
# 转换为DataFrame并选择前20行进行打印
# 选择特定的目标类别
target_class = 1
# 转换为DataFrame并选择前20行进行打印
shap_df = pd.DataFrame(shap_values2[target_class], columns=X_test.columns)
top20_shap_df = shap_df.head(20)
top20_shap_df.to_csv(pltpath + 'top20_shap_values2.csv', index=False)





#筛选分类模型
import pandas as pd
from sklearn.model_selection import cross_validate
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, BaggingClassifier, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

# 2 筛选机器学习模型
SEED=2022
models = [
    LogisticRegression(random_state=SEED),
    KNeighborsClassifier(),
    RandomForestClassifier(random_state=SEED),
    LinearDiscriminantAnalysis(),
    MLPClassifier(random_state=SEED),
    SVC(random_state=SEED),
    DecisionTreeClassifier(random_state=SEED),
    ExtraTreesClassifier(random_state=SEED),
    GradientBoostingClassifier(random_state=SEED),
    BaggingClassifier(random_state=SEED),
    AdaBoostClassifier(random_state=SEED)
]

scoring = {
    'accuracy': 'accuracy',
    'precision': 'precision',
    'recall': 'recall',
    'f1_score': 'f1',
    'auc_roc': 'roc_auc'
}

results = []
for model in models:
    cv_results = cross_validate(model, X_train, Y_train, cv=5, scoring=scoring)#k-fold交叉验证
    result = {
        'model': type(model).__name__,
        'accuracy_mean': cv_results['test_accuracy'].mean(),
        'accuracy_std': cv_results['test_accuracy'].std(),
        'precision_mean': cv_results['test_precision'].mean(),
        'precision_std': cv_results['test_precision'].std(),
        'recall_mean': cv_results['test_recall'].mean(),
        'recall_std': cv_results['test_recall'].std(),
        'f1_score_mean': cv_results['test_f1_score'].mean(),
        'f1_score_std': cv_results['test_f1_score'].std(),
        'auc_roc_mean': cv_results['test_auc_roc'].mean(),
        'auc_roc_std': cv_results['test_auc_roc'].std()
    }
    results.append(result)

    # 将结果导出为CSV文件
results_df = pd.DataFrame(results)
results_df.to_csv(pltpath + 'model_evaluation_results.csv', index=False)

import pandas as pd
from sklearn.model_selection import cross_validate
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier
from sklearn.neural_network import MLPClassifier
import xgboost as xgb
import lightgbm as lgb
from catboost import CatBoostClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

models = [
    LogisticRegression(),
    DecisionTreeClassifier(),
    SVC(),
    GaussianNB(),
    KNeighborsClassifier(),
    RandomForestClassifier(),
    GradientBoostingClassifier(),
    MLPClassifier(),
    AdaBoostClassifier(),
    xgb.XGBClassifier(),
    lgb.LGBMClassifier(),
    CatBoostClassifier()
]

# 创建评估指标列表
scoring = {
    'accuracy': 'accuracy',
    'precision': 'precision',
    'recall': 'recall',
    'f1_score': 'f1',
    'auc_roc': 'roc_auc'
}

# 执行交叉验证和模型评估
results = []
for model in models:
    cv_results = cross_validate(model, X_train, Y_train, cv=5, scoring=scoring)
    result = {
        'model': type(model).__name__,
        'accuracy_mean': cv_results['test_accuracy'].mean(),
        'accuracy_std': cv_results['test_accuracy'].std(),
        'precision_mean': cv_results['test_precision'].mean(),
        'precision_std': cv_results['test_precision'].std(),
        'recall_mean': cv_results['test_recall'].mean(),
        'recall_std': cv_results['test_recall'].std(),
        'f1_score_mean': cv_results['test_f1_score'].mean(),
        'f1_score_std': cv_results['test_f1_score'].std(),
        'auc_roc_mean': cv_results['test_auc_roc'].mean(),
        'auc_roc_std': cv_results['test_auc_roc'].std()
    }
    results.append(result)

# 将结果导出为CSV文件
results_df = pd.DataFrame(results)
results_df.to_csv(pltpath + '0512_model_evaluation_results.csv', index=False)


#画一张图捏
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.svm import SVC

# 创建模型列表
models = [
    LogisticRegression(random_state=SEED),
    KNeighborsClassifier(),
    RandomForestClassifier(random_state=SEED),
    LinearDiscriminantAnalysis(),
    MLPClassifier(random_state=SEED),
    SVC(random_state=SEED, probability=True),  # 设置probability为True
    DecisionTreeClassifier(random_state=SEED),
    ExtraTreesClassifier(random_state=SEED),
    GradientBoostingClassifier(random_state=SEED),
    BaggingClassifier(random_state=SEED),
    AdaBoostClassifier(random_state=SEED)
]

# 创建颜色列表
colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'magenta']

# 绘制AUC曲线
plt.figure(figsize=(10, 8))
for model, color in zip(models, colors):
    # 训练模型
    model.fit(X_train, Y_train)
    
    # 获取模型在测试集上的预测概率值
    y_scores = model.predict_proba(X_test)[:, 1]
    
    # 计算FPR和TPR
    fpr, tpr, _ = roc_curve(Y_test, y_scores)
    
    # 计算AUC
    roc_auc = auc(fpr, tpr)
    
    # 绘制ROC曲线
    plt.plot(fpr, tpr, color=color, label=f'{type(model).__name__} (AUC = {roc_auc:.2f})')

# 设置图例和标签
plt.legend(loc='lower right')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.savefig(pltpath + "ROC.jpg")


#保存模型
import pickle

# Save the models
with open(pltpath + 'logistic_regression_model.pkl', 'wb') as file:
    pickle.dump(lr_model, file)

with open(pltpath + 'decision_tree_model.pkl', 'wb') as file:
    pickle.dump(dt_model, file)

with open(pltpath + 'random_forest_model.pkl', 'wb') as file:
    pickle.dump(rf_model, file)

with open(pltpath + 'xgboost_model.pkl', 'wb') as file:
    pickle.dump(xgb_model, file)

with open(pltpath + 'lightgbm_model.pkl', 'wb') as file:
    pickle.dump(lgb_model, file)

