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
from pyimfcs.fitting import Fitter, make_fitp0
from pyimfcs.export import merge_fcs_results
from pyimfcs.blcorr import blexp_double_offset
import os
import pandas as pd
import datetime
from scipy.interpolate import interp1d

def process_stack(path,first_n = 0, last_n = 0, nsums=[2,3],
                           plot=False, default_dt= None, default_psize = None, 
                           fitter = None, export_summaries = True, 
                           chi_threshold = 0.03, ith=0.8):

    stack = StackFCS(path, background_correction = True,                     
                         first_n = first_n, last_n = last_n, clipval = 0)
    make_fitp0("2D",[lambda x:max(10**-3,1/x[0,1]), lambda x: D])
    stack.dt = default_dt
    stack.xscale = default_psize
    stack.yscale = default_psize
    xscale = stack.xscale
    yscale = stack.yscale
    assert(xscale==yscale)
    #stack.set_bleaching_function(blexp_double_offset)
    
    for nSum in nsums:
        stack.correlate_stack(nSum)
    if fitter is None:
        sigmaxy = sigma_psf*psize
        print("sigmaxy", sigmaxy)
        parameters_dict = {"mtype":"2D","a":yscale, "sigma":sigmaxy,"ginf":True}
        ft = Fitter(parameters_dict)
    else:
        ft = fitter
    
    stack.fit_curves(ft,xmax=None)
    
    stack.save(exclude_list=['traces_dict'])

def g2d(x0,y0,sigma):
    y,x=coords
    return np.exp(-( (x-x0)**2 + (y-y0)**2)/(2*sigma**2))

def coord2counts(x,y):
    #!!! in pixel coordinates
    frame = g2d(x,y,sigma_psf)
    return np.random.poisson(frame* brightness*dt)

def coord2counts_exact(x,y):
    #!!! in pixel coordinates
    frame = g2d(x,y,sigma_psf)
    real_counts = np.random.poisson(brightness*dt)
    psf_lin = frame.reshape(-1)
    psf_lin = np.cumsum(psf_lin) # cdf
    
    psf_lin/=psf_lin[-1]
    # !! Dabger
    fi = interp1d(psf_lin,np.arange(frame.size),fill_value="extrapolate")
    
    samples = np.random.uniform(low=0,high=1,size=real_counts)
    coords = fi(samples).astype(int) 
    
    out = np.zeros_like(psf_lin)
    for coord in coords:
        out[coord]+=1
    return out.reshape(frame.shape)

def simulate_2D_diff(D,nsteps,nparts, 
                     savepath = "/home/aurelienb/Data/simulations/SLB/",crop=5,
                     delete_tif=False):
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
    positions = np.random.uniform(size = (nparts,2))*npixels
    
    stack = np.zeros((nsteps,npix_img*2+1, npix_img*2+1))
    
    for j in range(nsteps):
        # in pixel space
        moves = np.random.normal(scale = np.sqrt(2*D*dt)/(psize),size = (nparts,2) )
        positions_new = positions+moves
        # round is necessary to ensure fair distribution of parts and not concentration in the centre

        positions_new = np.mod(positions_new,npixels)
        positions = positions_new.copy()
        
        positions_new = positions_new-npixels//2+npix_img
        positions_new = positions_new[np.logical_and(positions_new>=0,positions_new<npix_img*2+1).all(axis=1),:]
        for k in range(len(positions_new)):
            xx,yy=positions_new[k]
            stack[j]+=coord2counts(xx, yy)
        
        if j%500==0:
            print("Processing frame {}".format(j))
    stack=stack[:,crop:-crop,crop:-crop]
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
        "nparts": nparts,
        "npixels":npixels
        }
    
    parameters_df = pd.DataFrame(parameters_dict, index=[0])
    parameters_df.to_csv(savefolder+"parameters.csv")
    
    print('---- Processing FCS acquisition-----')
    # processing
    process_stack(stack_name, first_n = 0,
                           last_n = 0,nsums = [1,2], default_dt = dt, 
                           default_psize = psize)
    
    # export 
    intensity_threshold = 0
    thr = 0.03
    print([stack_name[:-4]+".h5"])
    merge_fcs_results(savefolder+"FCS_results",[stack_name[:-4]+".h5"], 
          chi_threshold = thr, ith = intensity_threshold)
    if delete_tif:
        os.remove(savefolder+"stack.tif")
        
# pixel size: 100 nm
psize = 0.1 #um
sigma_psf = 0.2/psize # pixels
dt = 10**-3 # s
D = 2 #um2/s

brightness = 18*10**4 #Hz/molecule

npixels = 500

npix_img = 4
coords = np.meshgrid(np.arange(2*npix_img+1),np.arange(2*npix_img+1))
nsteps = 50000
nparts = 5000

for npix in [15]:
    npixels = npix
    # parts per pixel square
    parts_density= 5000/(500**2)
    nparts = int(parts_density*npixels**2)
    for j in range(6):
        simulate_2D_diff(D,nsteps,nparts,crop=3,
             savepath= "/home/aurelienb/Data/simulations/SLB/2023_06_23_npixels/",delete_tif=True )
    