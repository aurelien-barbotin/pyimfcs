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

def plot_combined(combined,xname,measured,repeat):
    order = np.unique(combined[xname].values)
    # First measure
    ReplicateAverages = combined.groupby([xname,repeat], as_index=False).agg({measured: "mean"})
    
    fig, ax = plt.subplots(1,1)
    
    sns.violinplot(x=xname, y=measured, data=combined, order = order,ax=ax,color='gray',alpha=0.5)
    
    sns.swarmplot(x=xname, y=measured, hue=repeat, edgecolor="k", 
                       linewidth=2, data=ReplicateAverages, size=8, order = order,ax=ax)
    sns.swarmplot(x=xname, y=measured, hue=repeat, data=combined, order = order,ax=ax)
    ax.set_ylim(bottom=0)
    ax.legend_.remove()
    sns.reset_orig()
    

