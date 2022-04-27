#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 11:28:44 2021

@author: aurelien
"""

from PyImFCS.class_imFCS import (StackFCS, get_image_metadata, 
                                 bleaching_correct_sliding, bleaching_correct_exp)
import os
import matplotlib.pyplot as plt
import numpy as np

from PyImFCS.fitting import Fitter

plt.close("all")

path = "/home/aurelien/Documents/Data/2021_04_23/imFCS_tetraspeck/VerticalShiftSpeed0-5us/imFCS1.tif"

try:
    metadata = get_image_metadata(path)
    dt = metadata['finterval']
    xscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1']*10**6
    yscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1']*10**6

except:
    dt = 4.7*10**-3
    xscale = 0.1
    yscale = 0.1

stack = StackFCS(path, background_correction = True,
                 blcorrf = lambda x:bleaching_correct_sliding(x,wsize=1000), 
                 first_n=0, last_n = 0)
#lambda x:bleaching_correct_sliding(x,wsize=2000)
# bleaching_correction_exp

nsums = [1,2,3,4,5,8,16]

curves_avg = stack.binned_average_curves(nsums,n_norm=2)
plt.suptitle(os.path.split(path)[1].rstrip(".tif"))
#plt.savefig(savepath+"curves_{}.png".format(j))
trace = stack.trace()
plt.suptitle(os.path.split(path)[1].rstrip(".tif"))

# -- fitting ---

sigmaxy = 0.11
sigmaz = 0.48
# sigmaz = 100000000

parameters_dict = {"a":yscale, "sigmaxy":sigmaxy,'sigmaz':sigmaz}
ft = Fitter("3D",parameters_dict, ginf=True)
stack.fit_curves(ft,xmax=None)
stack.plot_random_intensity(nSum=4)
stack.plot_D()

stack.plot_fits(1,maxcurves = 10, dz=0.4)
stack.plot_fits(2,maxcurves = 10, dz=0.4)
stack.plot_fits(4,maxcurves = 10, dz=0.4)
stack.plot_fits(8,maxcurves = 10, dz=0.4)
