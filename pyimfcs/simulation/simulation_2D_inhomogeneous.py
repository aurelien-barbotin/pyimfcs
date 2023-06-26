# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:09:25 2021

@author: abarbotin
"""

import numpy as np
import tifffile

from pyimfcs.class_imFCS import StackFCS
from pyimfcs.fitting import Fitter, make_fitp0
from pyimfcs.export import merge_fcs_results
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
    #stack.set_bleaching_function(blexp_double_offset)
    
    make_fitp0("2D",[lambda x:max(10**-3,1/x[0,1]), lambda x: D])
    for nSum in nsums:
        stack.correlate_stack(nSum)
    if fitter is None:
        sigmaxy = sigma_psf
        print("sigmaxy", sigmaxy)
        parameters_dict = {"a":yscale, "sigma":sigmaxy,"mtype":"2D","ginf":True}
        ft = Fitter(parameters_dict)
    else:
        ft = fitter
    
    stack.fit_curves(ft,xmax=None)
    
    stack.save()

def g2d(x0,y0,sigma):
    y,x=coords
    return np.exp(-( (x-x0)**2 + (y-y0)**2)/(2*sigma**2))

def coord2counts(x,y):
    #!!! in pixel coordinates
    cxy= -npix_img
    rr = np.sqrt((x+cxy)**2+(y+cxy)**2)
    if rr>=R/psize:
        return 0
    frame = g2d(x,y,sigma_psf/psize)
    
    z = R*(1-np.sqrt(1-rr**2/(R/psize)**2))
    frame*=np.exp(-z/dz_tirf)
    return np.random.poisson(frame* brightness*dt)

def simulate_2D_diff(D,nsteps,nparts, 
                     savepath = "/home/aurelienb/Data/simulations/SLB/",crop=5):
    
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
    positions = np.random.uniform(size = (nparts,2))*npixels
    
    stack = np.zeros((nsteps,npix_img*2+1, npix_img*2+1))
    
    for j in range(nsteps):
        # in pixel space
        moves = np.random.normal(scale = np.sqrt(2*D*dt)/(psize),size = (nparts,2) )
        positions_new = positions+moves
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
        "npixels":npixels,
        "R":R,
        "dz_tirf":dz_tirf
        }
    
    parameters_df = pd.DataFrame(parameters_dict, index=[0])
    parameters_df.to_csv(savefolder+"parameters.csv")
    
    print('---- Processing FCS acquisition-----')
    # processing
    process_stack(stack_name, first_n = 0,
                           last_n = 0,nsums = [1,2,3,4,8], default_dt = dt, 
                           default_psize = psize)
    
    # export 
    intensity_threshold = 0.8
    thr = 0.03
    merge_fcs_results(savefolder+"FCS_results",[stack_name[:-4]+".h5"], 
          chi_threshold = thr, ith = intensity_threshold)

# pixel size: 100 nm
psize = 0.1 #um/pixels
sigma_psf = 0.2 # um
dt = 10**-3 # s
D = 1 #um2/s
dz_tirf = 0.1 #um
brightness = 18*10**4 #Hz/molecule

npixels = 500

npix_img = 24
coords = np.meshgrid(np.arange(2*npix_img+1),np.arange(2*npix_img+1))
nsteps = 20000
nparts = 5000

for j in range(5):
    for R in [0.5,1,2,5,10]:
        simulate_2D_diff(D,nsteps,nparts,crop=4,
             savepath= "/home/aurelienb/Data/simulations/SLB/2023_06_23_tirf_field/" )