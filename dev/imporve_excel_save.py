#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 11:33:00 2022

@author: aurelienb
"""

import pandas as pd

file = "/home/aurelienb/Documents/Projects/imFCS/results/BA/20deg/2022_09_21_+BA.xlsx"
excel = pd.ExcelFile(file)
df = excel.parse(sheet_name='nsum 3')

means = df.groupby('repeat')['D [µm²/s]'].mean()
std = df.groupby('repeat')['D [µm²/s]'].std()
medians = df.groupby('repeat')['D [µm²/s]'].median()
count = df.groupby('repeat')['D [µm²/s]'].count()
repeat = df.groupby('repeat')['repeat'].last()
filename = df.groupby('repeat')['filename'].last()
binning = df.groupby('repeat')['binning'].last()

out_df = pd.DataFrame.from_dict({
                    "filename": filename,
                    "binning": binning,
                    "repeat":repeat,
                    "Mean":means,
                    "Stdev": std,
                    "Median":medians,
                    "Count":count})
