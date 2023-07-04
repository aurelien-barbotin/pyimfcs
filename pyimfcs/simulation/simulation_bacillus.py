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

import matplotlib.pyplot as plt
from skimage.filters import gaussian
plt.close('all')

def process_stack(path,first_n = 0, last_n = 0, nsums=[2,3],
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
    frame = g2d(x,y,sigma_psf/psize)
    frame*=intensity_profile
    # frame *= g2d(npix_img//2,npix_img//2,R/psize)
    return np.random.poisson(frame* brightness*dt)

def simulate_2D_diff(D,nsteps,nparts, 
                     savepath = "/home/aurelienb/Data/simulations/SLB/"):
    
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
    positions = np.random.uniform(size = (nparts,2)) #positions x, y
    positions[:,0] *=npx
    positions[:,1]*=npy
    stack = np.zeros((nsteps,npx, npy))
    
    for j in range(nsteps):
        # in pixel space
        moves = np.random.normal(scale = np.sqrt(2*D*dt)/(psize),size = (nparts,2) )
        positions_new = positions+moves
        positions_new[:,0] = np.mod(positions_new[:,0],npx)
        positions_new[:,1] = np.mod(positions_new[:,1],npy)
        positions = positions_new.copy()
        
        positions_new[:,0]  = positions_new[:,0] -npx//2
        positions_new[:,1]  = positions_new[:,1] -npy//2
        
        for k in range(len(positions_new)):
            xx,yy=positions_new[k]
            stack[j]+=coord2counts(xx, yy)
        
        if j%500==0:
            print("Processing frame {}".format(j))
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
        }
    
    parameters_df = pd.DataFrame(parameters_dict, index=[0])
    parameters_df.to_csv(savefolder+"parameters.csv")
    
    print('---- Processing FCS acquisition-----')
    # processing
    process_stack(stack_name, first_n = 0,
                           last_n = 0,nsums = [1,2,3,4], default_dt = dt, 
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
brightness = 18*10**4 #Hz/molecule

npx = 50
npy = 31
coords = np.meshgrid(np.arange(npy)-npy/2,np.arange(npx)-npx/2) #shape y, x
nsteps = 20000
nparts = 200

length = 4 #um

sigma_tirf = 0.29 #um
intensity_profile = np.zeros((npx,npy))
len_pix=int(length//psize)
space=int((npx-len_pix)//2)
intensity_profile[space:len_pix+space,npy//2]=1

intensity_profile = gaussian(intensity_profile,sigma_tirf/psize)
intensity_profile/=intensity_profile.max()
# plt.figure()
# plt.imshow(intensity_profile)

for j in range(3):
    simulate_2D_diff(D,nsteps,nparts,
         savepath= "/home/aurelienb/Data/simulations/SLB/2023_07_01_bacillus/radiusOK/" )