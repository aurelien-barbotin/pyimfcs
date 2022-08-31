#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 16:06:15 2021

@author: aurelien
"""

import numpy as np
import matplotlib.pyplot as plt
from pyimfcs.class_imFCS import (StackFCS)
from pyimfcs.fitting import Fitter
from pyimfcs.plotting import multiplot_stack
import time
from pyimfcs.blcorr import blexp_double_offset
from pyimfcs.constants import datapath
t1 = time.time()
import glob
plt.close('all')

# Image 74_FCS_angle68.tif
path = datapath+"2021_06_08/Bacteria_NileRed/Image 31.tif"

load = False
save = True

def iqr_double(x):
    if len(x)<2:
        return np.array([0,0])
    return np.array([np.median(x)-np.percentile(x,25), np.percentile(x,75)-np.median(x)])

# binning 3, first_n 1000: 1 point OK. 3000: 3 pts OK , 5000: 5, 10000: 3, 800: 4
stack = StackFCS(path, background_correction = True,                     
                     first_n = 30000, last_n = 0, clipval = 0)

print('Exposure time {:.2f} ms'.format(stack.dt*1000))
print('Pixel size {} nm'.format(stack.xscale*1000))
xscale = stack.xscale
yscale = stack.yscale
# stack.registration(2000,plot=True)

try:
    if load:
        stack.load()
        save = False
except:
    print("no data could be loaded")
    
stack.set_bleaching_function(blexp_double_offset)
# stack.set_bleaching_function(bleaching_correct_sliding,wsize = 5000)
nsums=[2,3]
curves_avg = stack.binned_average_curves(nsums,n_norm=2)

sigmaxy = 0.2
parameters_dict = {"a":yscale, "sigma":sigmaxy}
ft = Fitter("2D",parameters_dict, ginf=True)

stack.fit_curves(ft,xmax=None)

chi_threshold = 0.03
multiplot_stack(stack,3, chi_threshold = chi_threshold)

if save:
    stack.save()
    print("saving stack")
    
