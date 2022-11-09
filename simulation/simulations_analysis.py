#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 16:45:09 2022

@author: aurelienb
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import glob

path = "/home/aurelienb/Data/simulations/Radius/"
folders = glob.glob(path+"*/")

nsum = 2

rs = []
ds = []

for folder in folders:
    df_params = pd.read_csv(folder+"parameters.csv")
    R = int(df_params['R'].values[0])
    print(R)
    df_excel = pd.read_excel(folder+"FCS_results.xlsx",sheet_name = "nsum {}".format(nsum))
    diffs = df_excel['D [µm²/s]'].values
    ds.extend(list(diffs))
    rs.extend(list(np.ones_like(diffs)*R))

rvals = sorted(np.unique(rs))
ds = np.array(ds)
rs = np.array(rs)

ds_mean = [ds[rs==w].mean() for w in rvals]
ds_std = [ds[rs==w].std() for w in rvals]

plt.figure()
plt.errorbar(rvals,ds_mean,yerr=ds_std,capsize=5)