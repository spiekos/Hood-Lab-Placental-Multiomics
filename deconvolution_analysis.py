import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, cdist, squareform
import composition_stats as cs
import matplotlib.pyplot as plt
import seaborn as sns
from utils import functions as f
import matplotlib.colors as colors
from matplotlib import cm
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from PIL import Image
import glob
import deconvolution_functions as defun

##### load in CIBERSORTx results, transcripts and metadata #####
df = pd.read_csv('cibersortx_results_placenta.csv', index_col = 0)
data = pd.read_csv('transcripts-final-filtered.csv', index_col = 0)
md = pd.read_csv('./metadata-final.csv', index_col = 0)

cols_fetal = [x for x in df.columns if f.like(x, ".*Fetal.*")]
cols_maternal = [x for x in df.columns if f.like(x, ".*Maternal.*")]

pats = df.index.tolist()
md = md.loc[pats, :]

######## get PCA of raw CIBERSORT results ########
X = df.values[:,:-3]
pca = PCA()
X_pca = pca.fit_transform(X)

######## transform the raw CIBERSORTx results ####### 
ep = 10e-6
Y = X.copy()
Y[Y==0] = ep
Y = Y/Y.sum(axis=1)[:,None]
Y = defun.clr(Y)
Y_pca = pca.fit_transform(Y)

######## plot to compare ####### 

file_path = '/PATH/TO/FILE/pca_compare.png'

fig, axes = plt.subplots(n_cols = 2, figsize = (20, 8))

axes[0].scatter(X_pca[:,0], X_pca[:,1])
axes[0].set_xlabel('$PC_1$', fontsize = 14)
axes[0].set_ylabel('$PC_2$', fontsize = 14)

axes[1].scatter(Y_pca[:,0], Y_pca[:,1])
axes[1].set_xlabel('$PC_1$', fontsize = 14)
axes[1].set_ylabel('$PC_2$', fontsize = 14)

plt.savefig(file_path)

######## get statistical analysis ####### 

file_name = '/PATH/TO/FILE/clr_results.csv'
dat = pd.DataFrame(Y, columns = df.columns[:-3], index = df.index)

df_result = run_tests(dat, df, md, method = 'control')
df_result.to_csv(file_name)

plot_results(df_result, df, dat, md.loc[pats, :])
# -- the above will save figures in working directory

##### plot the results as histograms. whisker plot is done in another notebook. This formats them all into a single image ####


im_path = glob.glob('./figures/*clr*')
im_path.sort()
images1 = [ Image.open(i) for i in im_path[:9]]
images2 = [ Image.open(i) for i in im_path[9:]]
images2.extend([Image.new("RGBA",(1200,900)) for i in range(2)])
#images2.extend([Image.open(im_path[-1])])
#images2.extend([Image.new("RGBA",(1200,900)) for i in range(1)])
print(len(images))

pil_grid(images1, max_horiz = 3)
pil_grid(images2, max_horiz = 4)

# ~~~~ fin ~~~~ #






