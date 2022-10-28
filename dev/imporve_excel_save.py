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

means = df.groupby('repeat').mean()
df.groupby('repeat').median()
df.groupby('repeat').std()