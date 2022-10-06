#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 11:42:06 2022

@author: aurelienb
"""

from pyimfcs.plotting import superplot_files
import matplotlib.pyplot as plt
plt.close('all')

path = "/home/aurelienb/Documents/Projects/collab Charlene/"
fileslist_list = [
    [path+"2022_10_05_control1.xlsx", path+"2022_10_05_control2.xlsx"],
    [path+"2022_10_05_flavo50+.xlsx"],
    [path+'2022_10_05_150ugmlflavo+.xlsx']
    ]

conditions = ["Control", "50 µg/mL", "150 µg/mL"]
superplot_files(fileslist_list,conditions,  keep_single_indices=True)
plt.savefig('swarmplot_results.png')