#这个代码用于解决类别不平衡问题，采用过采样和欠采样的方法
#类别不平衡问题
#https://zhuanlan.zhihu.com/p/36381828

#https://blog.csdn.net/weixin_41233157/article/details/131403617
#导入需要的包
import numpy as np
import pandas as pd
import pickle
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve, auc, precision_recall_curve


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
pro_df = pd.read_csv(dpath + 'Raw_data/231008_data_olink_instance_0.csv',usecols = ['ParticipantID1'] + pro_f_lst)
reg_df = pd.read_csv(dpath + 'Eid_info_data.csv', usecols = ['eid', 'Region_code'])
target_df = pd.read_csv(dpath + 'Disease_outcomes/IHD_control_definition/4without_nonIHD/IHD_outcomes.csv', usecols = ['eid', 'target_y', 'BL2Target_yrs'])
mydf = pd.merge(target_df, pro_df, how = 'inner', left_on = ['eid'], right_on = ['ParticipantID1'])
mydf = pd.merge(mydf, reg_df, how = 'left', on = ['eid'])
# 缺失值选择全部去掉，这其实不太好的，之后再找找文献看看有没有处理的办法（明天yebu一定）
mydf = mydf.dropna()
# 展示各个地区的人数情况
print(mydf['Region_code'].value_counts())
# 获取用于预测的数据和要预测的指标数据
X = mydf.drop(['eid','BL2Target_yrs','target_y','ParticipantID1'], axis=1)
X = X.reset_index(drop=True)
Y = mydf[['target_y','Region_code']]
Y = Y.reset_index(drop=True)
print('原始正样本数：', np.sum(Y['target_y'] == 1), '原始负样本数：', np.sum(Y['target_y'] == 0),   '原始总数：', len(X))

#过采样
smote = SMOTE()
x_new, y_new = smote.fit_resample(X, Y.drop(['Region_code'], axis=1))
print('smote后正样本数：', np.sum(y_new.values.flatten() == 1), 'smote后负样本数：', np.sum(y_new.values.flatten() == 0), 'smote后总数：', len(x_new))
#测试过采样后模型怎么样
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

model.fit(x_new, y_new.values.flatten())
y_pred = model.predict(x_new)
y_pred_prob = model.predict_proba(x_new)[:, 1]

#threshold = 0.8
#y_pred = (y_pred_prob > threshold).astype(int)

accuracy = accuracy_score(y_new, y_pred)
precision = precision_score(y_new, y_pred)
recall = recall_score(y_new, y_pred)
f1 = f1_score(y_new, y_pred)
roc_auc = roc_auc_score(y_new, y_pred_prob)
precision_cutoff, recall_cutoff, _ = precision_recall_curve(y_new, y_pred_prob)
pr_auc = auc(recall_cutoff, precision_cutoff)
evaluation = pd.DataFrame([accuracy,precision,recall,f1,roc_auc,pr_auc])
evaluation.columns = ['XGBoost']
evaluation.set_index(pd.Index(['Accuracy','Precision','recall','fi','roc_auc','pr_auc']),inplace=True)


#欠采样
rus = RandomUnderSampler()
x_new2, y_new2 = rus.fit_resample(X, Y.drop(['Region_code'], axis=1))
print('随机欠采样后正样本数：', np.sum(y_new2.values.flatten() == 1), '随机欠采样后负样本数：', np.sum(y_new2.values.flatten() == 0), '随机欠采样后总数：', len(x_new2))

#测试欠采样后模型怎么样
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
model.fit(x_new2, y_new2.values.flatten())
y_pred2 = model.predict(x_new2)
accuracy = accuracy_score(y_new2, y_pred2)
precision = precision_score(y_new2, y_pred2)
recall = recall_score(y_new2, y_pred2)
f1 = f1_score(y_new2, y_pred2)
y_pred_prob2 = model.predict_proba(x_new2)[:, 1]
roc_auc = roc_auc_score(y_new2, y_pred_prob2)
precision_cutoff, recall_cutoff, _ = precision_recall_curve(y_new2, y_pred_prob2)
pr_auc = auc(recall_cutoff, precision_cutoff)
evaluation = pd.DataFrame([accuracy,precision,recall,f1,roc_auc,pr_auc])
evaluation.columns = ['XGBoost']
evaluation.set_index(pd.Index(['Accuracy','Precision','recall','fi','roc_auc','pr_auc']),inplace=True)

