#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 18:16:35 2023

@author: aurelienb
"""
from pyimfcs.export import merge_fcs_results
from pyimfcs.class_imFCS import StackFCS

import glob
path="/home/aurelienb/Data/2022_12_13/2_44_coldshock/"
# path="/home/aurelienb/Data/2022_12_13/2_44_coldshock/Test/"
files=glob.glob(path+"*.h5")
merge_fcs_results("results_12_13",files,ith=0.8,chi_threshold=0.03,use_mask=True)

f="/home/aurelienb/Data/2022_12_13/2_44_coldshock/Image 32.h5"
stack=StackFCS(f,load_stack=False)
stack.load()
d1 = stack.extract_results(ith=0.8,chi_threshold=0.03,use_mask=True)[0][2]
d2 = stack.extract_results(ith=0.8,chi_threshold=0.03,use_mask=False)[0][2]
