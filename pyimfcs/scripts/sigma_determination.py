#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:56:02 2022

@author: aurelienb
"""
import numpy as np
import matplotlib.pyplot as plt
from pyimfcs.class_imFCS import StackFCS
from pyimfcs.fitting import Fitter
from pyimfcs.metrics import intensity_threshold
import glob
from scipy.stats import iqr

plt.close('all')

def get_d_curve(stack):
    parfits = stack.fit_results_dict
    ns = sorted(list(parfits.keys()))

    all_diffs_dict={}
    for n in ns:
        diffs=parfits[n][:,:,1]
        chisquares = stack.chisquares_dict[n]
        intensities = stack.thumbnails_dict[n]
        msk = chisquares<chi_th
        if ith>0:
            intth = intensity_threshold(ith,intensities)
            msk = np.logical_and(msk,intensities>intth)
            print('intensity threshold',intth)
        diffs = diffs[msk]
        
        all_diffs_dict[n]=diffs
    return all_diffs_dict

# This is a path to a folder containing several readily processed 
path = "path_to_dataset"
path = "/home/aurelienb/Data/RCL44_forpaper/20deg/2023_04_28_3_16h15_20deg/NR/"
files = glob.glob(path+"/*.h5")
print(files)
# Parameters that need tuning
sigmas = np.linspace(0.12,0.22,8)
chi_th=0.015
ith=0
# binning values you want to test these binning values must have been precomputed in the stacks of interest.
nsums = [1,2,3,4,5]

make_all_ds = lambda: dict(zip(nsums,[[] for w in nsums]))
all_ds_persigma = dict(zip(sigmas,[make_all_ds() for w in sigmas]))

for file in files:
    stack = StackFCS(file[:-3]+".tif",load_stack=False)
    stack.load(light_version=True)
    print("Pixel size: {} nm".format(int(stack.yscale*10**3)))
    for sigmaxy in sigmas:
        psize = stack.xscale
        parameters_dict = {"a":psize, "sigma":sigmaxy, "ginf":True,"mtype":"2D"}
        ft = Fitter(parameters_dict)
        
        stack.fit_curves(ft,xmax=None)
        ds=get_d_curve(stack)
        ns = sorted(list(ds.keys()))
        all_ds = all_ds_persigma[sigmaxy]
        for nn in nsums:
            all_ds[nn].extend(ds[nn])
            
plt.figure()
plt.subplot(121)
for sigmaxy in sigmas:
    all_ds = all_ds_persigma[sigmaxy]
    md = [np.median(all_ds[k]) for k in nsums]
    erd = [iqr(all_ds[k]) for k in nsums]
    
    plt.errorbar(nsums,md,yerr=erd,label=sigmaxy,capsize=5)

plt.xlabel('Binning')
plt.ylabel('D [µm²/s]')
plt.legend()
plt.title('Diffusion law for different sigma values')

plt.subplot(122)
for nn in nsums:
    curve = stack.average_curve(nn,chi_th=chi_th)
    plt.semilogx(curve[:,0],curve[:,1]/curve[0,1],label=nn)
plt.legend()
plt.xlabel(r'$\rm \tau$')
plt.xlabel(r'$\rm G(\tau)$')
plt.title('Averaged FCS curves in one stack')
plt.savefig('sigma_determination.png')
