#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:01:24 2022

@author: aurelien
"""


import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

import scipy

def plot_combined(combined,xname,measured,repeat, order = None,size = 8):
    """Plots results stored in a Dataframe using the SUperplots paradigm.
    
    Parameters:
        combined (Dataframe): dataframe containing merged results. 
        xname (str): name of conditions used in the x axis
        measured (str): name of the parameter to be measured (ex: D_[µm2/s])
        repeat (str): name of the repeat field
        order (list): list of condition names in the order the plots need to appear
        """
    if order is None:
        order = np.unique(combined[xname].values)
    # First measure
    ReplicateAverages = combined.groupby([xname,repeat], as_index=False).agg(
        {measured: "mean"})
    
    fig, ax = plt.subplots(1,1)
    
    sns.violinplot(x=xname, y=measured, data=combined, order = order,ax=ax,
                   color='gray',alpha=0.3)
    
    sns.swarmplot(x=xname, y=measured, hue=repeat, edgecolor="k", 
                       linewidth=2, data=ReplicateAverages, size=size, 
                       order = order,ax=ax)
    sns.swarmplot(x=xname, y=measured, hue=repeat, data=combined, order = order,ax=ax)
    ax.set_ylim(bottom=0)
    ax.legend_.remove()
    sns.reset_orig()
    

def superplot_files(files_list_list, conditions, nsum="nsum 3", keep_single_indices=True):
    
    if len(files_list_list)>1 and conditions is None:
        return ValueError('Please specify condition names')
    assert( conditions is None or len(conditions)==len(files_list_list))
    
    files = [w for flist in files_list_list for w in flist]
    excels = [pd.ExcelFile(w) for w in files]
    
    if conditions is not None:
        conditions_list = [conditions[j] for j, flist in 
                           enumerate(files_list_list) for w in flist]
    
    all_names = [w.sheet_names for w in excels]
    names0 = all_names[0]
    kept_names = list()
    for name in names0:
          if np.all([name in sublist for sublist in all_names]):
              kept_names.append(name)
    
    all_dfs = {}
    for name in kept_names:
        dfs = []
        maxindex = 0
        for j, xl in enumerate(excels):
            df = xl.parse(sheet_name=name, index_col = 0)
            if name=="parameters":
                fname = files[j]
                df["file"] = fname
                
            else:
                if keep_single_indices:
                    df['repeat']+=maxindex
                    maxindex = df['repeat'].values.max()+1
                    # print("maxindex", maxindex)
                else:
                    df['repeat'] = maxindex
                    maxindex += 1
            if conditions is not None:
                df["condition"] = conditions_list[j]
            dfs.append(df)
                
        dfs = pd.concat(dfs)
        all_dfs[name] = dfs
    xname = "condition"
    measured = "D [µm²/s]"
    plot_combined(all_dfs[nsum],xname,measured,'repeat')

if __name__=='__main__':
    import glob

    files_veg = glob.glob("/home/aurelienb/Documents/Projects/imFCS/Subtilis_veg/TFSM/*.xlsx")

    files_20 = [w for w in files_veg if '05_12' in w]
    files_37 = [w for w in files_veg if '05_12' not in w]
    files = [files_37, files_20]
    print(files)
    conditions = ["grown 37°C", "grown 20°C"]
    superplot_files(files,conditions)