import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import os
import cuml
import cudf
import array as arr
# import umap.umap_ as umap
import umap
import umap.plot

# import pycuda.autoinit
# import pycuda.gpuarray as gpuarray
# from skcuda.cublas import *
from scipy.sparse import csr_matrix
from sklearn.preprocessing import LabelBinarizer
# from matplotlib.gridspec import GridSpec
from cuml.manifold.umap import UMAP as clUMAP
# from cuml.cluster import KMeans as cuKM
# from cuml.decomposition import PCA as cuPCA
# from sklearn.decomposition import PCA as skPCA
# from sklearn.cluster import KMeans as skKM
from sklearn.preprocessing import LabelEncoder
%matplotlib inline

counts = pd.read_csv('/home/lima/Project/simulation/Comparision/all/normalized/NormalizedCountsCutOffValue0.csv')
counts = counts.rename(columns = {"Unnamed: 0": "GeneID"})

with open("CGRAllRSEM0Cut.pkl", "rb") as pkl_in:
    cgrCounts = pickle.load(pkl_in)
pkl_in.close()

convertedSeq = pd.DataFrame.from_dict(cgrCounts, orient = 'index')
seqList = convertedSeq.iloc[ : , 0].tolist()
idx = [i for i in range(len(seqList[0]))]
convertedSeq[idx] = pd.DataFrame(seqList, index = convertedSeq.index)
convertedSeq.reset_index(inplace = True)
convertedSeq.rename(columns = {"index": "GeneID"}, inplace = True)

samples = list(counts.columns.values)
samples = samples[1 : ]

geneIDs = list(convertedSeq['GeneID'])
cgrMatrix = convertedSeq.drop(columns = 'GeneID').to_numpy()
result = []

for sample in samples:
    # print(sample)
    sampleCounts = counts[sample].to_numpy()
    # print(sampleCounts.shape)
    sampleCounts = sampleCounts.reshape(len(sampleCounts), 1)

    tmp = cgrMatrix * sampleCounts
    colIdx = [i for i in range(tmp.shape[1])]
    labels = [sample] * len(tmp)
    # col = colIdx + lables

    df = pd.DataFrame(tmp, index = geneIDs)
    df['Label'] = labels
    result.append(df)

result = pd.concat(result)
rlt = result.drop(columns = 'Label')
rlt.reset_index(drop = True, inplace = True)
target = np.array(result.Label, dtype = np.str)
targetDF = pd.DataFrame(target)
encoder = LabelEncoder()
encodedLabels = encoder.fit_transform(target)
encodedLabelsDF = pd.Series(encodedLabels)
# arrEncodedLabels = arr.array('h', encodedLabels)
classesTotal = samples

gdf = cudf.from_pandas(rlt)
gdf_target = cudf.from_pandas(encodedLabelsDF

embedding_supervised = clUMAP(verbose = False, n_neighbors = 165, min_dist = 0.0, spread = 150, init = "spectral", target_metric = "categorical").fit_transform(gdf, gdf_target)
embedding_supervised_numpy = embedding_supervised.to_pandas().values
target = encodedLabels.astype(np.float32)
fig, ax = plt.subplots(1, figsize = (14, 10))
plt.scatter(embedding_supervised_numpy[ : , 1], embedding_supervised_numpy[ : , 0], s = 0.3, c = target, cmap = 'Spectral', alpha=1.0)
plt.setp(ax, xticks = [], yticks = [])
cbar = plt.colorbar(boundaries = np.arange(14) - 0.5)
cbar.set_ticks(np.arange(13))
cbar.set_ticklabels(samples)
plt.title('Supervised Embedded via cumlUMAP');

embedding_supervised = None
embedding_supervised = clUMAP(verbose = False, n_neighbors = 100, min_dist = 0.0, spread = 160, init = "spectral", target_metric = "categorical").fit_transform(gdf, gdf_target)
embedding_supervised_numpy = embedding_supervised.to_pandas().values
target = encodedLabels.astype(np.float32)
fig, ax = plt.subplots(1, figsize = (14, 10))
plt.scatter(embedding_supervised_numpy[ : , 1], embedding_supervised_numpy[ : , 0], s = 0.3, c = target, cmap = 'Spectral', alpha=1.0)
plt.setp(ax, xticks = [], yticks = [])
cbar = plt.colorbar(boundaries = np.arange(14) - 0.5)
cbar.set_ticks(np.arange(13))
cbar.set_ticklabels(samples)
plt.title('Supervised Embedded via cumlUMAP');
embedding_supervised = None
embedding_supervised = clUMAP(verbose = False, n_neighbors = 100, min_dist = 0.0, spread = 200, init = "spectral", target_metric = "categorical").fit_transform(gdf, gdf_target)
embedding_supervised_numpy = embedding_supervised.to_pandas().values
target = encodedLabels.astype(np.float32)
fig, ax = plt.subplots(1, figsize = (14, 10))
plt.scatter(embedding_supervised_numpy[ : , 1], embedding_supervised_numpy[ : , 0], s = 0.3, c = target, cmap = 'Spectral', alpha=1.0)
plt.setp(ax, xticks = [], yticks = [])
cbar = plt.colorbar(boundaries = np.arange(14) - 0.5)
cbar.set_ticks(np.arange(13))
cbar.set_ticklabels(samples)
plt.title('Supervised Embedded via cumlUMAP');
embedding_supervised = None
embedding_supervised = clUMAP(verbose = False, n_neighbors = 140, min_dist = 0.0, spread = 200, init = "spectral", target_metric = "categorical").fit_transform(gdf, gdf_target)
embedding_supervised_numpy = embedding_supervised.to_pandas().values
target = encodedLabels.astype(np.float32)
fig, ax = plt.subplots(1, figsize = (14, 10))
plt.scatter(embedding_supervised_numpy[ : , 1], embedding_supervised_numpy[ : , 0], s = 0.3, c = target, cmap = 'Spectral', alpha=1.0)
plt.setp(ax, xticks = [], yticks = [])
cbar = plt.colorbar(boundaries = np.arange(14) - 0.5)
cbar.set_ticks(np.arange(13))
cbar.set_ticklabels(samples)
plt.title('Supervised Embedded via cumlUMAP');
