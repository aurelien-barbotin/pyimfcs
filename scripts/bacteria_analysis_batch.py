#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 16:06:15 2021

@author: aurelien
"""

import matplotlib.pyplot as plt


from pyimfcs.class_imFCS import StackFCS
from pyimfcs.blcorr import blexp_double_offset

from pyimfcs.fitting import Fitter
from pyimfcs.hover_plot import multiplot_stack
import time
from pyimfcs.constants import datapath
import glob


t1 = time.time()

plt.close('all')
# fig 8: 10 good points
def batch_bacteria_process(files,first_n = 3000, last_n = 0):
    for path in files:
        print('Processing',path)
        stack = StackFCS(path, background_correction = True,                     
                             first_n = first_n, last_n = last_n, clipval = 0)
        
        print('Exposure time {:.2f} ms'.format(stack.dt*1000))
        print('Pixel size {} nm'.format(stack.xscale*1000))
        xscale = stack.xscale
        yscale = stack.yscale
        assert(xscale==yscale)
        stack.registration(4000,plot=True)
            
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
    
        stack.save()
        print("saving stack")

path = datapath+"/2022_08_30/2_20deg/*.tif"
files = glob.glob(path)
print(files)
files = list(filter(lambda x: "GP" not in x,files))
batch_bacteria_process(files, first_n=0, last_n = 0)

t2 = time.time()
dt = t2-t1
print("Elapsed time {:} min {:.2f} s".format(dt//60,dt%60))
