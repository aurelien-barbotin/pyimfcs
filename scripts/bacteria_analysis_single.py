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
from pyimfcs.plotting import interactive_fcs_plot
import time
from pyimfcs.blcorr import blexp_double_offset
from pyimfcs.constants import datapath
t1 = time.time()
import glob
plt.close('all')

path = datapath+"/2023_03_14/2_44_37/2p5X/Image 63.tif"

load = False
save = False


# binning 3, first_n 1000: 1 point OK. 3000: 3 pts OK , 5000: 5, 10000: 3, 800: 4
stack = StackFCS(path, background_correction = True,                     
                     first_n = 0, last_n = 0, clipval = 0)

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
nsums=[3]
curves_avg = stack.binned_average_curves(nsums,n_norm=2)

sigmaxy = 0.2
parameters_dict = {"a":yscale, "sigma":sigmaxy,"ginf":True, "mtype":"2D"}
ft = Fitter(parameters_dict)

stack.fit_curves(ft,xmax=None)

chi_threshold = 0.03
interactive_fcs_plot(stack,nsums[-1], chi_threshold = chi_threshold,light_version=False)

if save:
    stack.save()
    print("saving stack")
    
