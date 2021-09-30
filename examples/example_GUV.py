#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 16:06:15 2021

@author: aurelien
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import tifffile

from PyImFCS.class_imFCS import (StackFCS, get_image_metadata, 
                                 bleaching_correct_sliding,
                                 blexp_double_offset, bleaching_correct_segment)
from PyImFCS.fitting import Fitter
from PyImFCS.hover_plot import multiplot_stack
from PyImFCS.shift_correction import registration

from skimage.filters import threshold_otsu
import time

t1 = time.time()


plt.close('all')

path="/home/aurelien/Data/2021_07_30/GUV2/Image 7_0p5excitation.tif"
try:
    metadata = get_image_metadata(path)
    dt = metadata['finterval']
    xscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1']*10**6
    yscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1']*10**6
    
except:
    dt = 4.7*10**-3
    xscale = 0.1
    yscale = 0.1

def iqr_double(x):
    if len(x)<2:
        return np.array([0,0])
    return np.array([np.median(x)-np.percentile(x,25), np.percentile(x,75)-np.median(x)])

print('Exposure time {:.2f} ms'.format(dt*1000))
print('Pixel size {} nm'.format(xscale*1000))
stack = StackFCS(path, background_correction = True,
                     blcorrf = lambda x: bleaching_correct_segment(x, wsize=5000), 
                     first_n = 0, last_n = 0, clipval = 5000)

stack.stack = registration(stack.stack,1000,plot=True)

load = False
save = True
try:
    if load:
        stack.load()
        save = False
except:
    print("no data could be loaded")
    
nsums =[2,4,8]


curves_avg = stack.binned_average_curves(nsums,n_norm=2)
sigmaxy = 0.2
parameters_dict = {"a":yscale, "sigma":sigmaxy}
ft = Fitter("2D",parameters_dict, ginf=True)

stack.fit_curves(ft,xmax=None)
if save:
    stack.save()
    print("saving stack")
simg = stack.stack.mean(axis=0)

nss = 4
all_ds = []
chi_threshold = 0.04

for nss in nsums:
    chimap = stack.chisquares_dict[nss]<chi_threshold
    thmap = stack.get_threshold_map(nss, thf = lambda x:threshold_otsu(x))
    
    dmap = np.logical_and(thmap,chimap)
    ds = stack.parfit_dict[nss][dmap,1]
    all_ds.append(ds)

all_ds = [w[w>0] for w in all_ds]
median_ds = [np.median(w) for w in all_ds]
iqrs_double = np.array([iqr_double(w) for w in all_ds]).T

plt.figure()
plt.errorbar(nsums,median_ds, yerr=iqrs_double)
plt.xlabel("Binning (pixels)")
plt.ylabel(r"$\rm D\ [\mu m^2/s]$")

stack.plot_fits_ordered(4,maxcurves=10,order_dict=stack.chisquares_dict)
plt.suptitle('Ordered curves, binning 4')

stack.plot_fits_ordered(2,maxcurves=10,order_dict=stack.chisquares_dict)
plt.suptitle('Ordered curves, binning 2')
    
stack.plot_parameter_maps(nsums[:5],parn=0)
plt.suptitle("Number of molecules")

stack.plot_parameter_maps(nsums[:5],parn=1, maxval=15)
plt.suptitle("Diffusion Coefficients")

stack.plot_random_intensity()

stack.plot_intensity_correlation(2)
plt.suptitle('Binning 2')

nsum = 4
chis = stack.chisquares_dict[nsum].reshape(-1)
ds = stack.parfit_dict[nsum][:,:,1].reshape(-1)

plt.figure()

plt.scatter(chis,ds)
plt.xlabel('residuals')
plt.ylabel("D [um2/s]")

stack.plot_taus()

trace = stack.stack[:,5:9,5:9].mean(axis=(1,2))
bleaching_correct_segment(trace,plot=True,wsize=5000)

multiplot_stack(stack,8)
print('elapsed time : {}'.format(time.time()-t1))
