# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:09:25 2021

@author: abarbotin
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import gaussian
import tifffile

from pyimfcs.class_imFCS import StackFCS
from pyimfcs.fitting import Fitter
from pyimfcs.io import merge_fcs_results
from pyimfcs.blcorr import blexp_double_offset
import os
import pandas as pd
import datetime

def process_stack(path,first_n = 3000, last_n = 0, nsums=[2,3],
                           plot=False, default_dt= None, default_psize = None, 
                           fitter = None, export_summaries = True, 
                           chi_threshold = 0.03, ith=0.8):

    stack = StackFCS(path, background_correction = True,                     
                         first_n = first_n, last_n = last_n, clipval = 0)

    stack.dt = default_dt
    stack.xscale = default_psize
    stack.yscale = default_psize
    xscale = stack.xscale
    yscale = stack.yscale
    assert(xscale==yscale)
    stack.set_bleaching_function(blexp_double_offset)
    
    for nSum in nsums:
        stack.correlate_stack(nSum)
    if fitter is None:
        sigmaxy = sigma_psf*psize
        parameters_dict = {"a":yscale, "sigma":sigmaxy}
        ft = Fitter("2D",parameters_dict, ginf=True)
    else:
        ft = fitter
    
    stack.fit_curves(ft,xmax=None)
    
    stack.save()
 
def simulate_2D_diff(D,nsteps,nparts, 
                     savepath = "/home/aurelienb/Data/simulations/SLB/"):
    pos0 = np.random.uniform(size = (nparts,2))*npixels-npixels/2
    moves = np.random.normal(scale = np.sqrt(2*D*dt)/(psize),size = (nsteps,nparts,2) )
    
    positions = np.cumsum(moves,axis=0)
    positions = positions+pos0[np.newaxis,:,:]
    
    npix_img = 16
    stack = np.zeros((nsteps,npix_img*2+1, npix_img*2+1))
    
    for j in range(nsteps):
        
        positions_new = positions[j]
        # round is necessary to ensure fair distribution of parts and not concentration in the centre
        positions_new = np.round(positions_new).astype(int) + npix_img
        positions_new = positions_new[np.logical_and(positions_new>=0,positions_new<npix_img*2+1).all(axis=1),:]
        for k in range(len(positions_new)):
            stack[j, positions_new[k, 0],positions_new[k, 1]]+=1
        stack[j] = gaussian(stack[j],sigma = sigma_psf)
        
        if j%500==0:
            print("Processing frame {}".format(j))
    stack = stack * brightness*dt 
    stack = np.random.poisson(stack)
    x = np.linspace(-(npix_img*2+1)/2,(npix_img*2+1)/2,npix_img*2+1)
    
    xx, yy = np.meshgrid(x,x)

    stack = stack.astype(np.uint16)
    #------------ MSD stuff ---------------
    
    fname = datetime.datetime.now().__str__()[:19]
    savefolder = savepath+fname+"/"
    os.mkdir(savefolder)
    stack_name = savefolder+"stack.tif"
    tifffile.imsave(stack_name,stack)
    
    parameters_dict = {
        "psize": psize,
        "sigma_psf":sigma_psf, 
        "dt": dt,
        "D":D,
        "brightness": brightness,
        "nsteps": nsteps,
        "nparts": nparts
        }
    
    parameters_df = pd.DataFrame(parameters_dict, index=[0])
    parameters_df.to_csv(savefolder+"parameters.csv")
    
    print('---- Processing FCS acquisition-----')
    # processing
    process_stack(stack_name, first_n = 0,
                           last_n = 0,nsums = [1,2,3,4,8], default_dt = dt, 
                           default_psize = psize)
    
    # export 
    intensity_threshold = 0
    thr = 0.03
    merge_fcs_results([stack_name[:-4]+".h5"], savefolder+"FCS_results", 
          intensity_threshold = intensity_threshold, chi_threshold = thr)

# pixel size: 100 nm
psize = 0.16 #um
sigma_psf = 0.1/psize # pixels
dt = 10**-3 # s
D = 5 #um2/s

brightness = 80*10**3 #Hz/molecule

npixels = 200

nsteps = 10000
nparts = 40000
for dd in [2]:
    simulate_2D_diff(dd,nsteps,nparts)