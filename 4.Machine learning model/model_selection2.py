# 这个代码进行机器学习模型的贝叶斯调参

#导入需要的包
import numpy as np
import pandas as pd
import pickle
#机器学习用于cross validation和模型效果指征的包
import sklearn.model_selection as ms
from sklearn.model_selection import cross_val_score, KFold
from sklearn import metrics
from sklearn.metrics import roc_auc_score, make_scorer
#机器学习classifier的包
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
#贝叶斯调参的包
from bayes_opt import BayesianOptimization
#网格搜索调参的包
from sklearn.model_selection import GridSearchCV

#读取数据集
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
target_df = pd.read_csv(dpath + 'Disease_outcomes/IHD_control_definition/4without_nonIHD/IHD_outcomes.csv', usecols = ['eid', 'target_y', 'BL2Target_yrs'])
mydf = pd.merge(target_df, pro_df, how = 'inner', left_on = ['eid'], right_on = ['ParticipantID1'])
#缺失值选择用每列平均值插补，因为有些模型不允许缺失值
mydf = mydf.fillna(mydf.mean())

#获取用于预测的数据和要预测的指标数据
X = mydf.drop(['eid','BL2Target_yrs','target_y','ParticipantID1'], axis=1)
X = X.reset_index(drop=True)
Y = mydf[['target_y']]
Y = Y.reset_index(drop=True)




# RF贝叶斯调参
# https://zhuanlan.zhihu.com/p/131216861
# 5-fold crossed validation
rf = RandomForestClassifier()
def rf_cv(n_estimators, min_samples_split, max_features, max_depth):
    val = cross_val_score(
        RandomForestClassifier(n_estimators=int(n_estimators),
            min_samples_split=int(min_samples_split),
            max_features=min(max_features, 0.999), # float
            max_depth=int(max_depth), class_weight="balanced",
            random_state=0
        ),
        X, Y.values.flatten(), scoring='roc_auc', cv=5
    ).mean() # https://www.javaroad.cn/questions/9775解释了为什么这样计算出来是对的
    return val

rf_bo = BayesianOptimization(
        rf_cv,
        {'n_estimators': (10, 250),
        'min_samples_split': (2, 25),
        'max_features': (0.1, 0.999),
        'max_depth': (5, 15)}
    )

# 指定保存文件的路径和名称，保存优化得到的参数
rf_bo.maximize(n_iter=100)
opt_parameter = rf_bo.max
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/random_forest.pickle'
# 保存模型到 pickle 文件
with open(file_path, 'wb') as f:
    pickle.dump(opt_parameter, f)


# DT贝叶斯调参
# https://blog.csdn.net/weixin_41885239/article/details/121870123
# https://blog.csdn.net/weixin_54814385/article/details/122652761
# https://blog.csdn.net/SKIp121whats112/article/details/122265766?spm=1001.2014.3001.5502
# 5-fold crossed validation
dt = tree.DecisionTreeClassifier()
def dt_cv(max_depth,min_samples_split,min_samples_leaf,min_impurity_decrease):
    val = cross_val_score(
        tree.DecisionTreeClassifier(max_depth=int(max_depth),
            min_samples_split=int(min_samples_split),
            min_samples_leaf=int(min_samples_leaf),
            min_impurity_decrease=max(min_impurity_decrease, 0.0),
            class_weight="balanced",
            random_state=0
        ),
        X, Y, scoring='roc_auc', cv=5
    ).mean() # https://www.javaroad.cn/questions/9775解释了为什么这样计算出来是对的
    return val

dt_bo = BayesianOptimization(
        dt_cv,
        {'max_depth': (10, 150),
        'min_samples_split': (2, 15),
        'min_samples_leaf': (1,25),
        'min_impurity_decrease': (0.0, 0.5)}
    )

# 指定保存文件的路径和名称，保存优化得到的参数
dt_bo.maximize(n_iter=500)
opt_parameter = dt_bo.max
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/dicision_tree.pickle'
# 保存模型到 pickle 文件
with open(file_path, 'wb') as f:
    pickle.dump(opt_parameter, f)


# LR贝叶斯调参
# https://blog.csdn.net/evolution23/article/details/85028361
# https://blog.csdn.net/qq_51205385/article/details/127785470
# 5-fold crossed validation
lr = LogisticRegression()
def lr_cv(penalty, C):
    if penalty <= 0.5:
        penalty_val = 'l1'
        solver = 'liblinear' # 因为L1不支持lbfgs求解器
    else:
        penalty_val = 'l2'
        solver = 'lbfgs'
    val = cross_val_score(
        LogisticRegression(penalty=penalty_val,
            C=C,
            solver=solver,
            class_weight="balanced",
            random_state=0
        ),
        X, Y.values.flatten(), scoring='roc_auc', cv=5
    ).mean() # https://www.javaroad.cn/questions/9775解释了为什么这样计算出来是对的
    return val

lr_bo = BayesianOptimization(
        lr_cv,
        {'penalty': (0, 1),
        'C': (0.01, 100)}
    )

# 指定保存文件的路径和名称，保存优化得到的参数
lr_bo.maximize(n_iter=500)
opt_parameter = lr_bo.max
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/logistic_regression.pickle'
# 保存模型到 pickle 文件
with open(file_path, 'wb') as f:
    pickle.dump(opt_parameter, f)


# xgb贝叶斯调参
# https://blog.csdn.net/weixin_44839513/article/details/108699862
# 5-fold crossed validation
xgb = XGBClassifier()
def xgb_cv(max_depth, learning_rate, n_estimators, min_child_weight, subsample, colsample_bytree, reg_alpha,gamma):
    val = cross_val_score(
        XGBClassifier(max_depth=int(max_depth),
            learning_rate=learning_rate,
            n_estimators=int(n_estimators),
            min_child_weight=min_child_weight,
            subsample=max(min(subsample, 1), 0),
            colsample_bytree=max(min(colsample_bytree, 1), 0),
            reg_alpha=max(reg_alpha, 0), gamma=gamma, objective='reg:squarederror',
            booster='gbtree',
            scale_pos_weight =11,
            random_state=0
        ),
        X, Y, scoring='roc_auc', cv=5
    ).mean() # https://www.javaroad.cn/questions/9775解释了为什么这样计算出来是对的
    return val

xgb_bo = BayesianOptimization(
        xgb_cv,
        {'max_depth': (1, 10),
         'learning_rate': (0.001, 0.1),
         'n_estimators': (1, 1000),
         'min_child_weight': (0, 20),
         'subsample': (0.001, 1),
         'colsample_bytree': (0.01, 1),
         'reg_alpha': (0.001, 20),
         'gamma': (0.001, 10)}
    )

# 指定保存文件的路径和名称，保存优化得到的参数
xgb_bo.maximize(n_iter=100)
opt_parameter = xgb_bo.max
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/XGBoost.pickle'
# 保存模型到 pickle 文件
with open(file_path, 'wb') as f:
    pickle.dump(opt_parameter, f)


# LGBM贝叶斯调参
# https://zhuanlan.zhihu.com/p/361832566
# 5-fold crossed validation
lgbm = LGBMClassifier()
def lgbm_cv(n_estimators,min_gain_to_split,subsample, max_depth,colsample_bytree, min_child_samples,reg_alpha,reg_lambda,num_leaves,learning_rate):
    val = cross_val_score(
        LGBMClassifier(boosting_type='gbdt', objective='binary',
            colsample_bytree=float(colsample_bytree),
            min_child_samples=int(min_child_samples),
            n_estimators=int(n_estimators),
            num_leaves=int(num_leaves),
            reg_alpha=float(reg_alpha),
            reg_lambda=float(reg_lambda),
            max_depth=int(max_depth),
            subsample=float(subsample),
            min_gain_to_split = float(min_gain_to_split),
            learning_rate=float(learning_rate),
            scale_pos_weight =11,
            random_state=0,
            force_col_wise=True,
            verbosity=-1#避免输出太多奇怪的信息
        ),
        X, Y, scoring='roc_auc', cv=5
    ).mean() # https://www.javaroad.cn/questions/9775解释了为什么这样计算出来是对的
    return val

lgbm_bo = BayesianOptimization(
        lgbm_cv,
        {'colsample_bytree': (0.5,1),
         'min_child_samples': (2, 200),
         'num_leaves': (5, 1000),
         'subsample': (0.6, 1),
         'max_depth':(2,20),
         'n_estimators': (10, 1000),
         'reg_alpha':(0,10),
         'reg_lambda':(0,10),
         'min_gain_to_split':(0,1),
         'learning_rate':(0.001,0.1)
        }
    )

# 指定保存文件的路径和名称，保存优化得到的参数
lgbm_bo.maximize(n_iter=500)
opt_parameter = lgbm_bo.max
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/LGBM.pickle'
# 保存模型到 pickle 文件
with open(file_path, 'wb') as f:
    pickle.dump(opt_parameter, f)


# SVM贝叶斯调参
# https://blog.csdn.net/weixin_48618536/article/details/131483817
# https://blog.csdn.net/qq_41076797/article/details/101037721
# https://blog.csdn.net/ztf312/article/details/98594359
# 5-fold crossed validation
svm = SVC()
def svm_cv(C,gamma):
    val = cross_val_score(
        SVC(kernel='rbf',
            C=C,
            gamma=gamma,
            class_weight="balanced",
            random_state=0
        ),
        X, Y.values.flatten(), scoring='roc_auc', cv=5
    ).mean() # https://www.javaroad.cn/questions/9775解释了为什么这样计算出来是对的
    return val

svm_bo = BayesianOptimization(
        svm_cv,
        {'C': (0.01, 1000),
        'gamma': (0, 0.5)}
    )

# 指定保存文件的路径和名称，保存优化得到的参数
svm_bo.maximize(n_iter=100)
opt_parameter = svm_bo.max
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/support_vector_machine.pickle'
# 保存模型到 pickle 文件
with open(file_path, 'wb') as f:
    pickle.dump(opt_parameter, f)


# kNN贝叶斯调参
# https://zhuanlan.zhihu.com/p/656629354
# https://blog.csdn.net/Joker_sir5/article/details/82258368
# 5-fold crossed validation
knn = KNeighborsClassifier()
def knn_cv(n_neighbors,p):
    val = cross_val_score(
        KNeighborsClassifier(weights='distance',
            n_neighbors=int(n_neighbors),
            p=int(p)
        ),
        X, Y.values.flatten(), scoring='roc_auc', cv=5
    ).mean() # https://www.javaroad.cn/questions/9775解释了为什么这样计算出来是对的
    return val

knn_bo = BayesianOptimization(
        knn_cv,
        {'n_neighbors': (10,1500),  
         'p': (1, 4)}
    )

# 指定保存文件的路径和名称，保存优化得到的参数
knn_bo.maximize(n_iter=100)
opt_parameter = knn_bo.max
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/K_nearest_neighbor.pickle'
# 保存模型到 pickle 文件
with open(file_path, 'wb') as f:
    pickle.dump(opt_parameter, f)


# LDA随即搜索调参
# https://www.cnblogs.com/solong1989/p/9593555.html
# LDA似乎没有很多要调整的地方，感觉可以不用贝叶斯调参调整，在这里选用了随机搜索调参
# 5-fold crossed validation
# LDA的拟合效果非常差，还是别用了
lda = LinearDiscriminantAnalysis()
param_space = {
    'solver': ['svd', 'lsqr', 'eigen'],
    'shrinkage': [None, 'auto'],
    'n_components': [None, 1, 2, 3, 4, 5]  # 可以根据需要扩展其他参数
}

scorer = make_scorer(roc_auc_score)
grid_search = GridSearchCV(lda, param_grid=param_space, cv=5, scoring=scorer)
grid_search.fit(X, Y.values.flatten())
# 输出最佳参数
print("Best parameters found: ", grid_search.best_params_)
# 输出最佳得分
print("Best ROC-AUC score found: ", grid_search.best_score_)

# 指定保存文件的路径和名称，保存优化得到的参数
opt_parameter = grid_search.best_params_
file_path = '/home/linmiao/ZhangjunYu/Proteomics_analysis/Result/model_selection/model_parameter/linear_discriminant_analysis.pickle'
# 保存模型到 pickle 文件
with open(file_path, 'wb') as f:
    pickle.dump(opt_parameter, f)

