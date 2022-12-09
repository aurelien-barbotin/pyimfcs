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
    #stack.set_bleaching_function(blexp_double_offset)
    
    for nSum in nsums:
        stack.correlate_stack(nSum)
    if fitter is None:
        sigmaxy = sigma_psf*psize
        print("sigmaxy", sigmaxy)
        parameters_dict = {"a":yscale, "sigma":sigmaxy}
        ft = Fitter("2D",parameters_dict, ginf=True)
    else:
        ft = fitter
    
    stack.fit_curves(ft,xmax=None)
    
    stack.save()

def g2d(x0,y0,sigma):
    y,x=coords
    return np.exp(-( (x-x0)**2 + (y-y0)**2)/(2*sigma**2))

def coord2counts(x,y):
    #!!! in pixel coordinates
    frame = g2d(x,y,sigma_psf)
    return np.random.poisson(frame* brightness*dt)

def simulate_2D_diff(D,nsteps,nparts, 
                     savepath = "/home/aurelienb/Data/simulations/SLB/",crop=5,
                     bleach_radius=2):
    
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
    positions = np.random.uniform(size = (nparts,2))*npixels
    
    dists_pos = np.sqrt((positions[:,0]-npixels/2)**2+(positions[:,1]-npixels/2)**2)
    positions=positions[dists_pos>bleach_radius/psize,:]
    print(positions.shape)
    nparts_afterbleaching = positions.shape[0]
    x = np.linspace(-(npix_img*2+1)/2,(npix_img*2+1)/2,npix_img*2+1)
    xx, yy = np.meshgrid(x,x)
    
    dists = np.sqrt(xx**2+yy**2)*psize
    mask = dists<bleach_radius
    
    plt.figure()
    plt.imshow(mask)
    
    stack = np.zeros((nsteps,npix_img*2+1, npix_img*2+1))
    
    for j in range(nsteps):
        # in pixel space
        moves = np.random.normal(scale = np.sqrt(2*D*dt)/(psize),size = (nparts_afterbleaching,2) )
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
        "bleach_radius":bleach_radius
        }
    
    parameters_df = pd.DataFrame(parameters_dict, index=[0])
    parameters_df.to_csv(savefolder+"parameters.csv")
    

# pixel size: 100 nm
psize = 0.16 #um
sigma_psf = 0.1/psize # pixels
dt = 10*10**-3 # s
D = 2 #um2/s

brightness = 18*10**3 #Hz/molecule

npixels = 100

npix_img = 24
coords = np.meshgrid(np.arange(2*npix_img+1),np.arange(2*npix_img+1))
nsteps = 100
nparts = 50000

simulate_2D_diff(D,nsteps,nparts,crop=1,
                 savepath= "/home/aurelienb/Data/simulations/FRAP/" )