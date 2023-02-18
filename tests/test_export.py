#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 18:16:35 2023

@author: aurelienb
"""
from pyimfcs.export import merge_fcs_results
import glob
path="/home/aurelienb/Data/2022_12_13/2_44_coldshock/"
# path="/home/aurelienb/Data/2022_12_13/2_44_coldshock/Test/"
files=glob.glob(path+"*.h5")
merge_fcs_results("results_12_13",files,ith=0.8,chi_threshold=0.03,use_mask=False)