#完全走一遍师兄推荐的代码
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
data = pd.read_csv("/share/home/zhangjunyu/Project/240918_metabolic_multimorbidity_clusters/Data/Covariates/blood_biochemistry_na_remove.csv")
data = data.drop(columns=['eid'])

#方差阈值
from sklearn.feature_selection import VarianceThreshold
# Create a VarianceThreshold instance with a threshold of 0.25
var_thr = VarianceThreshold(threshold=0.25)
# Fit the VarianceThreshold to the data
var_thr.fit(data)
# Get the support (boolean mask) for the selected features
var_thr.get_support()
# Identify the constant and quasi-constant columns
con_cols = [col for col in data.columns if col not in data.columns[var_thr.get_support()]]
# Drop the constant and quasi-constant columns from the dataset
data = data.drop(con_cols, axis=1)
#删除了['Apolipoprotein_A', 'Apolipoprotein_B', 'HDL_cholesterol']

#特征缩放
from sklearn.preprocessing import StandardScaler
# Initialize the StandardScaler
scaler = StandardScaler()
scaler.fit(data)
scaled_data = pd.DataFrame(scaler.transform(data), columns=data.columns)

#降维
from sklearn.decomposition import PCA
# Initialize PCA with the number of components equal to the number of scaled features
pca = PCA(n_components=scaled_data.shape[1])
pca.fit(scaled_data)
# Calculate the explained variance and cumulative explained variance
explained_variance = pca.explained_variance_ratio_
cumulative_explained_variance = np.cumsum(explained_variance)
# Visualize the explained variance to select the number of components
import matplotlib.pyplot as plt
plt.figure(figsize=(8,6))
plt.plot(cumulative_explained_variance, marker='o', linestyle='--')
plt.xlabel('Number of Components')
plt.ylabel('Cumulative Explained Variance')
plt.title('Explained Variance vs. Number of Components')
plt.grid(True)
plt.savefig('/share/home/zhangjunyu/Project/240918_metabolic_multimorbidity_clusters/Result/cluster_shixiong/explained_variance.pdf', format='pdf')
plt.close()
# Select the number of components based on the desired explained variance threshold
desired_explained_variance = 0.95
num_components = np.argmax(cumulative_explained_variance >= desired_explained_variance) + 1
# Output the selected number of components
print(f"Number of components to explain {desired_explained_variance:.2f} variance: {num_components}")
# Fit PCA with the selected number of components
pca = PCA(n_components=num_components)
PCA_ds = pd.DataFrame(pca.fit_transform(scaled_data), columns=[f"PC{i+1}" for i in range(num_components)])

#elbow method
from yellowbrick.cluster import KElbowVisualizer
from sklearn.cluster import KMeans
# Initialize the KElbowVisualizer with the KMeans estimator and a range of K values
Elbow_M = KElbowVisualizer(KMeans(), k=10)
# Fit the visualizer to the PCA-transformed data
Elbow_M.fit(PCA_ds)
Elbow_M.show(outpath="/share/home/zhangjunyu/Project/240918_metabolic_multimorbidity_clusters/Result/cluster_shixiong/elbow_method.pdf")
plt.close()

# Specify the number of clusters and other parameters
n_clusters = 4
max_iter = 1000
algorithm = 'elkan'
# Create a KMeans instance
kmeans_4 = KMeans(n_clusters=n_clusters, 
                max_iter=max_iter, 
                random_state=2024, 
                algorithm=algorithm)

# Fit the KMeans model to your data
yhat_KM_k4 = kmeans_4.fit_predict(PCA_ds)
# Assign cluster labels to your DataFrame
PCA_ds_k4 = PCA_ds.copy()
data_k4 = data.copy()
PCA_ds_k4["Cluster"] = yhat_KM_k4
data_k4["Cluster"] = yhat_KM_k4

#轮廓系数分析
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
# Define the range of K values to explore
k_values = range(2, 10)
silhouette_scores = []
best_k = None
best_silhouette_score = -1
# Iterate through different K values
for k in k_values:
    km = KMeans(n_clusters=k, init='k-means++', n_init=10, max_iter=100, random_state=2024)
    km.fit(PCA_ds)
    # Calculate the average silhouette score
    silhouette_avg = silhouette_score(PCA_ds, km.labels_)  
    silhouette_scores.append(silhouette_avg)
    # Update the best K value if a higher silhouette score is found
    if silhouette_avg > best_silhouette_score:
        best_k = k
        best_silhouette_score = silhouette_avg

# Output the best K value and corresponding Silhouette Score
print("Best K value:", best_k)
print("Best Silhouette Score:", best_silhouette_score)
import matplotlib.pyplot as plt
# Create a plot to visualize Silhouette Scores
plt.plot(k_values, silhouette_scores)
plt.xlabel("Number of Clusters (K)")
plt.ylabel("Silhouette Score")
plt.title("Silhouette Score vs. Number of Clusters")
plt.savefig('/share/home/zhangjunyu/Project/240918_metabolic_multimorbidity_clusters/Result/cluster_shixiong/Silhouette_Score.pdf', format='pdf')
plt.close()

#最终轮廓系数最优的类的数目是2
# Specify the number of clusters and other parameters
n_clusters = 2
max_iter = 1000
algorithm = 'elkan'
# Create a KMeans instance
kmeans_2 = KMeans(n_clusters=n_clusters, 
                max_iter=max_iter, 
                random_state=2024, 
                algorithm=algorithm)

# Fit the KMeans model to your data
yhat_KM_k2 = kmeans_2.fit_predict(PCA_ds)
# Assign cluster labels to your DataFrame
PCA_ds_k2 = PCA_ds.copy()
data_k2 = data.copy()
PCA_ds_k2["Cluster"] = yhat_KM_k2
data_k2["Cluster"] = yhat_KM_k2

#SSE评估k-means聚类有效性
sse_k4 = kmeans_4.inertia_ 
sse_k2 = kmeans_2.inertia_ 
print(f"When K=4 SSE = {sse_k4}\nWhen K=2 SSE = {sse_k2}")

#方差比准则（Calinski-Harabasz Index）
from sklearn.metrics import calinski_harabasz_score
# 范围内的聚类数目 (k值)
k_values = range(2, 11)
# 存储每个k值的CH Index
ch_scores = []
# 计算每个k值的Calinski-Harabasz Index
for k in k_values:
    kmeans = KMeans(n_clusters=k, random_state=2024)
    labels = kmeans.fit_predict(PCA_ds)
    ch_index = calinski_harabasz_score(PCA_ds, labels)
    ch_scores.append(ch_index)

# 找到CH Index最大的k值
best_k = k_values[np.argmax(ch_scores)]
# 输出最佳k值
print(f"最佳聚类数目为: {best_k}")
# 绘制CH Index随k值的变化
plt.figure(figsize=(8, 5))
plt.plot(k_values, ch_scores, marker='o')
plt.xlabel("聚类数目 (k)")
plt.ylabel("Calinski-Harabasz Index")
plt.title("不同聚类数目的Calinski-Harabasz Index")
plt.grid(True)
plt.savefig('/share/home/zhangjunyu/Project/240918_metabolic_multimorbidity_clusters/Result/cluster_shixiong/calinski_harabasz_score.pdf', format='pdf')
plt.close()

#聚类结果可视化-PCA
from sklearn.decomposition import PCA
#k=2
kmeans = kmeans_2
labels = yhat_KM_k2
pca = PCA(n_components=2)
X_pca = pca.fit_transform(PCA_ds)
plt.figure(figsize=(8, 6))
scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=labels, cmap='viridis', marker='o', s=50)
# 增加聚类中心
centers_pca = pca.transform(kmeans.cluster_centers_)
plt.scatter(centers_pca[:, 0], centers_pca[:, 1], c='red', marker='x', s=200, label='Cluster Centers')
plt.title('PCA 2D Visualization of KMeans Clustering')
plt.xlabel('PCA Component 1')
plt.ylabel('PCA Component 2')
plt.legend()
plt.grid(True)
plt.savefig('/share/home/zhangjunyu/Project/240918_metabolic_multimorbidity_clusters/Result/cluster_shixiong/k2_visual_PCA.png', format='png')
plt.close()
#k=4
kmeans = kmeans_4
labels = yhat_KM_k4
pca = PCA(n_components=2)
X_pca = pca.fit_transform(PCA_ds)
plt.figure(figsize=(8, 6))
scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=labels, cmap='viridis', marker='o', s=50)
# 增加聚类中心
centers_pca = pca.transform(kmeans.cluster_centers_)
plt.scatter(centers_pca[:, 0], centers_pca[:, 1], c='red', marker='x', s=200, label='Cluster Centers')
plt.title('PCA 2D Visualization of KMeans Clustering')
plt.xlabel('PCA Component 1')
plt.ylabel('PCA Component 2')
plt.legend()
plt.grid(True)
plt.savefig('/share/home/zhangjunyu/Project/240918_metabolic_multimorbidity_clusters/Result/cluster_shixiong/k4_visual_PCA.png', format='png')
plt.close()

#聚类结果可视化-PCA3维
#k=2
from mpl_toolkits.mplot3d import Axes3D
kmeans = kmeans_2
labels = yhat_KM_k2
pca = PCA(n_components=3)
X_reduced = pca.fit_transform(PCA_ds)
centers_pca = pca.transform(kmeans.cluster_centers_)
# 创建 3D 可视化
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# 绘制数据点
sc = ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=labels, cmap='viridis')
# 绘制聚类中心
ax.scatter(centers_pca[:, 0], centers_pca[:, 1], centers_pca[:, 2], 
           c='red', marker='X', s=500, linewidths=2, label='Cluster Centers')
ax.set_title('3D Visualization of Clustering with Cluster Centers')
ax.set_xlabel('PCA Component 1')
ax.set_ylabel('PCA Component 2')
ax.set_zlabel('PCA Component 3')
plt.savefig('/share/home/zhangjunyu/Project/240918_metabolic_multimorbidity_clusters/Result/cluster_shixiong/k2_visual_PCA3.png', format='png')
plt.close()
#k=4
from mpl_toolkits.mplot3d import Axes3D
kmeans = kmeans_4
labels = yhat_KM_k4
pca = PCA(n_components=3)
X_reduced = pca.fit_transform(PCA_ds)
centers_pca = pca.transform(kmeans.cluster_centers_)
# 创建 3D 可视化
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# 绘制数据点
sc = ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=labels, cmap='viridis')
# 绘制聚类中心
ax.scatter(centers_pca[:, 0], centers_pca[:, 1], centers_pca[:, 2], 
           c='red', marker='X', s=500, linewidths=2, label='Cluster Centers')
ax.set_title('3D Visualization of Clustering with Cluster Centers')
ax.set_xlabel('PCA Component 1')
ax.set_ylabel('PCA Component 2')
ax.set_zlabel('PCA Component 3')
plt.savefig('/share/home/zhangjunyu/Project/240918_metabolic_multimorbidity_clusters/Result/cluster_shixiong/k4_visual_PCA3.png', format='png')
plt.close()
