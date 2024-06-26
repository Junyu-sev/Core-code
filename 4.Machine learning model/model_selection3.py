# 这个代码进行模型评估
# https://blog.csdn.net/LegenDavid/article/details/79063044
# https://zhuanlan.zhihu.com/p/552278753
# 似乎AUC才是更稳健，更加好的评估手段

#导入需要的包
import numpy as np
import pandas as pd
import pickle
#机器学习用于cross validation和模型效果指征的包
import sklearn.model_selection as ms
from sklearn.model_selection import cross_val_score, KFold
from sklearn import metrics
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, auc, precision_recall_curve
#机器学习classifier的包
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
# import sklearn.svm as svm

# 获取所有的数据集，用于建模
# 读取数据集
dpath = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Data/'
outpath = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/'
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
#pro_f_lst = ["GDF15","HAVCR1","MMP12","PGF","EDA2R","NTproBNP","BCAN","CXCL17","CDCP1","XG","CNTN5","PRSS8"]

pro_df = pd.read_csv(dpath + 'Raw_data/231008_data_olink_instance_0.csv',usecols = ['ParticipantID1'] + pro_f_lst)
reg_df = pd.read_csv(dpath + 'Eid_info_data.csv', usecols = ['eid', 'Region_code'])
target_df = pd.read_csv(dpath + 'Disease_outcomes/IHD_control_definition/4without_nonIHD/IHD_outcomes.csv', usecols = ['eid', 'target_y', 'BL2Target_yrs'])
mydf = pd.merge(target_df, pro_df, how = 'inner', left_on = ['eid'], right_on = ['ParticipantID1'])
mydf = pd.merge(mydf, reg_df, how = 'left', on = ['eid'])
# 缺失值选择全部去掉，这其实不太好的，之后再找找文献看看有没有处理的办法（明天yebu一定）
#mydf = mydf.dropna()
# 缺失值选择取平均值处理
column_means = mydf.mean()
mydf = mydf.fillna(column_means)
# 展示各个地区的人数情况
print(mydf['Region_code'].value_counts())
#加入一些神奇的协变量试试效果，加入年龄和性别对于precision没有特别好的改变
#cov_df = pd.read_csv('/home/linmiao/ZhangjunYu/Proteomics_analysis/Data/Raw_data/covariate_before_inputation.csv', usecols= ['eid','age','sex'])
# 获取用于预测的数据和要预测的指标数据
X = mydf.drop(['eid','BL2Target_yrs','target_y','ParticipantID1'], axis=1)
X = X.reset_index(drop=True)
Y = mydf[['target_y','Region_code']]
Y = Y.reset_index(drop=True)


#创建一个函数，根据不同的模型计算accuracy，precision，recall，F1，AUC（ROC和PR）
#内外部交叉验证，按照UKB英国地区进行划分，每个地区依次作为验证数据集，其他作为训练集，然后计算AUC（C-statistic）和overall C-statistic
#展示每个地区的人数，一共10个地区
def evaluate_model(model):
    #定义需要计算的指标
    accuracies, precisions, recalls, f1s, roc_aucs, pr_aucs = [], [], [], [], [], []
    #循环每个地区，region based corssed validation
    for i in range(10):
        X_train = X.loc[X['Region_code'] != i].drop(['Region_code'], axis=1)
        Y_train = Y.loc[Y['Region_code'] != i].drop(['Region_code'], axis=1)
        X_test =  X.loc[X['Region_code'] == i].drop(['Region_code'], axis=1)
        Y_test = Y.loc[Y['Region_code'] == i].drop(['Region_code'], axis=1)
        model.fit(X_train, Y_train.values.flatten())
        y_pred = model.predict(X_test)
        y_pred_prob = model.predict_proba(X_test)[:, 1]
        accuracies.append(accuracy_score(Y_test, y_pred))
        precisions.append(precision_score(Y_test, y_pred))
        recalls.append(recall_score(Y_test, y_pred))
        f1s.append(f1_score(Y_test, y_pred))
        roc_aucs.append(roc_auc_score(Y_test, y_pred_prob))
        precision_cutoff, recall_cutoff, _ = precision_recall_curve(Y_test, y_pred_prob)
        pr_aucs.append(auc(recall_cutoff, precision_cutoff))
    evaluation_results = pd.DataFrame({'Accuracy': accuracies, 'Precision': precisions, 'recall': recalls, 'f1': f1s, 'roc_auc' : roc_aucs, 'pr_auc' : pr_aucs})
    mean_row = evaluation_results.mean()
    sd_row = evaluation_results.std()    
    evaluation_results.loc[len(evaluation_results)] = mean_row
    evaluation_results.loc[len(evaluation_results)] = sd_row
    evaluation_results.loc[len(evaluation_results)] = mean_row - 1.96*sd_row
    evaluation_results.loc[len(evaluation_results)] = mean_row + 1.96*sd_row
    evaluation_results = evaluation_results.T
    evaluation_results.columns = ['East_Midlands','London','North_East','North_West','Scotland','South_East','South_West','Wales','West_Midlands','Yorkshire_and_Humber','MEAN','SD','lower_95ci','upper_95ci']
    evaluation_results['Confidence_Interval'] = evaluation_results.apply(lambda row: f"{row['MEAN']:.2f} ({row['lower_95ci']:.2f}-{row['upper_95ci']:.2f})", axis=1)
    return evaluation_results

#RF的计算
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/random_forest.pickle'
with open(file_path, 'rb') as f:
    params = pickle.load(f)

params = params['params']
model = RandomForestClassifier(n_estimators=int(params['n_estimators']),
            min_samples_split=int(params['min_samples_split']),
            max_features=min(params['max_features'], 0.999), # float
            max_depth=int(params['max_depth']), class_weight="balanced")

outcome = evaluate_model(model)
outcome.to_csv("/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_evaluation/Random_forest.csv")
rf_evaluation = outcome['Confidence_Interval']

#dt的计算
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/dicision_tree.pickle'
with open(file_path, 'rb') as f:
    params = pickle.load(f)

params = params['params']
model = tree.DecisionTreeClassifier(max_depth=int(params['max_depth']),
            min_samples_split=int(params['min_samples_split']),
            min_samples_leaf=int(params['min_samples_leaf']),
            min_impurity_decrease=max(params['min_impurity_decrease'], 0.0),
            class_weight="balanced")

outcome = evaluate_model(model)
outcome.to_csv("/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_evaluation/Decision_tree.csv")
dt_evaluation = outcome['Confidence_Interval']

#LR的计算
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/logistic_regression.pickle'
with open(file_path, 'rb') as f:
    params = pickle.load(f)

params = params['params']
if params['penalty'] <= 0.5:
        penalty_val = 'l1'
        solver = 'liblinear' # 因为L1不支持lbfgs求解器
else:
        penalty_val = 'l2'
        solver = 'lbfgs'

model = LogisticRegression(penalty=penalty_val,
            C=params['C'],
            solver=solver,
            class_weight="balanced"
        )

outcome = evaluate_model(model)
outcome.to_csv("/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_evaluation/Logistic_regression.csv")
lr_evaluation = outcome['Confidence_Interval']

#xgb的计算
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/XGBoost.pickle'
with open(file_path, 'rb') as f:
    params = pickle.load(f)

params = params['params']
model = XGBClassifier(max_depth=int(params['max_depth']),
            learning_rate=params['learning_rate'],
            n_estimators=int(params['n_estimators']),
            min_child_weight=params['min_child_weight'],
            subsample=max(min(params['subsample'], 1), 0),
            colsample_bytree=max(min(params['colsample_bytree'], 1), 0),
            reg_alpha=max(params['reg_alpha'], 0), gamma=params['gamma'], objective='reg:squarederror',
            booster='gbtree',
            scale_pos_weight =11
        )

outcome = evaluate_model(model)
outcome.to_csv("/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_evaluation/XGBoost.csv")
xgb_evaluation = outcome['Confidence_Interval']

#lgbm的计算
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/LGBM.pickle'
with open(file_path, 'rb') as f:
    params = pickle.load(f)

params = params['params']
model = LGBMClassifier(boosting_type='gbdt', objective='binary',
            colsample_bytree=float(params['colsample_bytree']),
            min_child_samples=int(params['min_child_samples']),
            n_estimators=int(params['n_estimators']),
            num_leaves=int(params['num_leaves']),
            reg_alpha=float(params['reg_alpha']),
            reg_lambda=float(params['reg_lambda']),
            max_depth=int(params['max_depth']),
            subsample=float(params['subsample']),
            min_gain_to_split = float(params['min_gain_to_split']),
            learning_rate=float(params['learning_rate']),
            scale_pos_weight =11,
            random_state=0,
            force_col_wise=True
        )

outcome = evaluate_model(model)
outcome.to_csv("/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_evaluation/LGBM.csv")
lgbm_evaluation = outcome['Confidence_Interval']

combined_evaluation = pd.concat([rf_evaluation, dt_evaluation, lr_evaluation, xgb_evaluation, lgbm_evaluation], axis=1, ignore_index=True)
combined_evaluation.columns = ['random_forest','decision_tree','logistic_regression','XGBoost','LightGBM']
combined_evaluation.to_csv("/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_evaluation/summary_evaluation.csv")

