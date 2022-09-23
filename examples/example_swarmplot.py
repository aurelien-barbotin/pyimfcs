#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 11:42:06 2022

@author: aurelienb
"""

from pyimfcs.plotting import superplot_files
import matplotlib.pyplot as plt
plt.close('all')

path = "/home/aurelienb/Documents/Projects/imFCS/results/sporulation/SSB1002_sporulation/50k/"
fileslist_list = [
    [path+"2022_06_08_TFSM_W01_prespo_50k.xlsx", path+"2022_06_08_TFSM_W03_prespo_50k.xlsx",
     path+"2022_06_08_TFSM_W03_prespo_slide2_50k.xlsx"],
    [path+"2022_06_09_1_W01_beforeSporulation_TFSM_50k.xlsx"],
    [path+"2022_06_09_3_W01_TFSM_50k.xlsx",path+"/2022_06_09_5_W01_TFSM_17h_50k.xlsx"]
    ]
conditions = ["Before sporulation", "pre sporulation", "sporulation"]
superplot_files(fileslist_list,conditions,  keep_single_indices=False)
