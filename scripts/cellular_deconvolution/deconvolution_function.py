import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from statsmodels.stats.multitest import multipletests
from PIL import Image
import glob

######################################################################
####################function definitions #############################
######################################################################

## CoDA functions
def clr(X):
    # X is an np.array with imputed zeros
    X1 = X.copy()
    g = np.exp(np.log(X1).mean())
    X1[X1==1] = g
    return np.log(X1/g)

def ilr(X):
    # X is an np.array with imputed zero
    clr_x = clr(X)
    Psi = np.ones(shape = (X.shape[1]-1,X.shape[1]))
    for i in range(X.shape[1]-1):
        Psi[i,i] = np.exp(1)
    psi_p = np.linalg.pinv(Psi)
    return clr_x @ psi_p

### functions for running a bunch of KS tests inb the clr-transformed results
def run_tests(df,df1, md, col = 'Condition', adjust = True, method = '1vr'):
    grps = set(md[col].tolist())
    preDf =[]
    for g in grps:
        if (method =='control')&(g == 'Control'):
            continue
        pats = md.query(f'Condition == "{g}"').index
        if method == '1vr':
            comp2 = 'rest'
            no_pats = [p for p in md.index if p not in pats]
        if method == 'control':
            comp2 = 'control'
            no_pats = md.query('Condition == "Control"').index
        chrt1 = df.loc[pats,:]
        chrt2 = df.loc[no_pats,:]
        
        chrt1_per = df1.loc[pats,:]
        chrt2_per = df1.loc[no_pats,:]
        for c in df.columns:
            vals1 = chrt1[c].values
            vals2 = chrt2[c].values
            
            vals1_per = chrt1_per[c].values
            vals2_per = chrt2_per[c].values
            
            k, p = stats.ks_2samp(vals1, vals2)
            preDf.append([g,c, k, p, vals1.mean(), vals2.mean(), vals1.std(), vals2.std(), vals1_per.mean(), vals2_per.mean(), vals1_per.std(), vals2_per.std()])
    cols = columns = ['group', 'cell_type', 'KS statistic', 'pvalue', 'group_mean_clr', f'{comp2}_mean_clr', 'group_std_clr', f'{comp2}_std_clr','group_mean_per', f'{comp2}_mean_per', 'group_std_per', f'{comp2}_std_per']
    df_out =pd.DataFrame(preDf,columns =cols )
    if adjust:
        df_out['pvalue_adj'] = multipletests(df_out['pvalue'].values,method ='fdr_bh')[1]
    return df_out

### functions for plotting results

def plot_results(results, df_per, df_clr, md, alpha = .05, show = True, save = True):
    grps = results.groupby('group')
    if 'Control' not in grps.groups.keys():
        other = 'Control'
    else:
        other = 'Rest'
    for _, df in grps:
        if other == 'Control':
            pats2 = md.query('Condition == "Control"').index
        else:
            pats2 = md.query(f'Condition != "{_}"').index
        pats1 = md.query(f'Condition == "{_}"').index
        sig = df.query(f'pvalue_adj <={alpha}').loc[:, ['cell_type', 'KS statistic', 'pvalue_adj']].values
        if len(sig) ==0:
            continue
        for cell_type, k, p in sig: 
            vals1_per = df_per.loc[pats1, cell_type].values
            vals2_per = df_per.loc[pats2, cell_type].values
            
            vals1_clr = df_clr.loc[pats1, cell_type].values
            vals2_clr = df_clr.loc[pats2, cell_type].values
            
            h, bins = np.histogram(np.concatenate([vals1_per, vals2_per]), bins = 15, density = True)
            
            plt.hist(vals1_per, bins = bins, color = 'red', histtype = 'step', label = _, density = True)
            plt.hist(vals2_per, bins = bins, color = 'black', histtype = 'step', label = other, density = True)
            plt.xlabel(f'% {cell_type}', fontsize = 13)
            plt.ylabel(f'density', fontsize = 13)
            plt.title(f'KS = {np.round(k, 4)}, p_adj = {np.round(p, 4)}')
            plt.legend()
            if save == True:
                plt.savefig(f'./figures/{_}_{cell_type}_per.png', dpi = 200)
            if show == True:
                plt.show()
                
                
            h, bins = np.histogram(np.concatenate([vals1_clr, vals2_clr]), bins = 15, density = True)
            
            plt.hist(vals1_clr, bins = bins, color = 'red', histtype = 'step', label = _, density = True)
            plt.hist(vals2_clr, bins = bins, color = 'black', histtype = 'step', label = other, density = True)
            plt.xlabel(f'CLR(% {cell_type} ) ', fontsize = 13)
            plt.ylabel(f'density', fontsize = 13)
            plt.legend()
            plt.title(f'KS = {np.round(k, 4)}, p_adj = {np.round(p, 4)}')
            if save == True:
                plt.savefig(f'./figures/{_}_{cell_type}_clr.png', dpi = 200)
            if show == True:
                plt.show()

def pil_grid(images, max_horiz=np.iinfo(int).max):
    """ stolen from stack overflow:
    https://stackoverflow.com/questions/30227466/combine-several-images-horizontally-with-python
    """
    n_images = len(images)
    n_horiz = min(n_images, max_horiz)
    print(n_images,n_horiz)
    h_sizes, v_sizes = [0] * n_horiz, [0] * (n_images // n_horiz)
    print(h_sizes, v_sizes)
    for i, im in enumerate(images):
        h, v = i % n_horiz, i // n_horiz
        h_sizes[h] = max(h_sizes[h], im.size[0])
        v_sizes[v] = max(v_sizes[v], im.size[1])
    h_sizes, v_sizes = np.cumsum([0] + h_sizes), np.cumsum([0] + v_sizes)
    im_grid = Image.new('RGB', (h_sizes[-1], v_sizes[-1]), color='white')
    for i, im in enumerate(images):
        im_grid.paste(im, (h_sizes[i % n_horiz], v_sizes[i // n_horiz]))
    return im_grid



