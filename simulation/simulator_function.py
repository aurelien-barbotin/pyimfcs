# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:09:25 2021
@author: abarbotin
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import gaussian
import tifffile
import datetime

from pyimfcs.class_imFCS import StackFCS
from pyimfcs.fitting import Fitter
from pyimfcs.io import merge_fcs_results
plt.close('all')

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
    from SO: https://stackoverflow.com/questions/13685386/
    matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    
def gauss(x,sig):
    return np.exp(-2*x**2/(sig**2))

def spherical2cart(R,phi,theta):
    x = R*np.cos(phi)*np.sin(theta)
    y = R*np.sin(phi)*np.sin(theta)
    z = R*np.cos(theta)
    return np.array([x, y, z])

def cart2spherical(x,y,z):
    theta = np.arccos(z)
    phi = np.arctan2(y,x)
    return np.array([phi, theta])


def coord2counts(x,y,z):
    zr = np.sqrt(2*np.log(2))*sigmaz
    omegaz = sigma_psf*np.sqrt(1+z/zr**2)
    frame = np.zeros((npix_img*2+1, npix_img*2+1))
    frame[x,y] = np.exp(-z/dz_tirf)
    
    frame = gaussian(frame,sigma = omegaz)*(sigma_psf/omegaz)**2
    return np.random.poisson(frame* brightness*dt )

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
    # stack.set_bleaching_function(bleaching_correct_sliding,wsize = 5000)
    
    for nSum in nsums:
        stack.correlate_stack(nSum)
    if fitter is None:
        sigmaxy = 0.2
        parameters_dict = {"a":yscale, "sigma":sigmaxy}
        ft = Fitter("2D",parameters_dict, ginf=True)
    else:
        ft = fitter
    
    stack.fit_curves(ft,xmax=None)
    
    stack.save()

plot = True
save = True

psize = 0.16
sigma_psf = 0.2/psize
sigmaz = 4*sigma_psf
dz_tirf = 0.2 # um

dt = 1*10**-3 # s
D = 3 #um2/s

R = 8 #um
brightness = 18*10**3 #Hz/molecule

npix_img = 16
z_cutoff = 10*dz_tirf

nsteps = 20000
nparts = 1000

def simulate_spherical_diffusion(R,D,nsteps,nparts):
    pos0 = np.random.uniform(size = (nparts,2))*2*np.pi
    
    # ---------- Calculation of positions-------------------
    moves = np.random.normal(scale = np.sqrt(2*D*dt)/R,
                             size = (nsteps,nparts,2) )
    
    moves[0] = 0
    
    positions = np.cumsum(moves,axis=0)
    positions = positions+pos0[np.newaxis,:,:]
    
    p1 = pos0[:,0]
    phis = moves[:,:,0]/np.sin(positions[:,:,1])
    positions[:,:,0] = p1[np.newaxis,:] + np.cumsum(phis,axis=0)
    
    positions = positions%(2*np.pi)
    z,y,x = spherical2cart(R, positions[:,:,0], positions[:,:,1])
    
    # ---------- making image ------------------
    stack = np.zeros((nsteps,npix_img*2+1, npix_img*2+1))
    
    for j in range(nsteps):
        if j%500==0:
            print("Processing frame {}".format(j))
            
        x1, y1, z1 = x[j], y[j], z[j]
        z1+=R
        # round is necessary to ensure fair distribution of parts and not concentration in the centre
        positions_new = np.array([x1,y1]).T/psize + npix_img
        positions_new = np.round(positions_new).astype(int) 
        
        msk0 = np.logical_and(positions_new>=0,positions_new<npix_img*2+1).all(axis=1)
        msk3 = np.logical_and(msk0,z1<z_cutoff)
        positions_new = positions_new[msk3,:]
        znew = z1[msk3]
        for k in range(positions_new.shape[0]):
            frame = coord2counts(positions_new[k,0], positions_new[k,1],znew[k])
            stack[j]+=frame
    
    """plt.figure()
    plt.subplot(221)
    plt.imshow(stack[0])
    plt.title('Frame 0')
    plt.subplot(222)
    plt.imshow(stack[2000])
    plt.title('Frame 2000')
    plt.subplot(223)
    plt.imshow(stack[9999])
    plt.title('Frame 9999')
    plt.subplot(224)
    plt.imshow(stack.sum(axis=0))
    plt.title('Sum of all frames')
    plt.suptitle('Summary of simulated acquisition')"""
    
    
    import os
    import pandas as pd
    #------------ MSD stuff ---------------
    
    savepath = "/home/aurelienb/Data/simulations/"
    fname = datetime.datetime.now().__str__()[:19]
    savefolder = savepath+fname+"/"
    os.mkdir(savefolder)
    stack_name = savefolder+"stack.tif"
    tifffile.imsave(stack_name,stack)
    
    parameters_dict = {
        "psize": psize,
        "sigma_psf":sigma_psf, 
        "sigmaz": sigmaz,
        "dz_tirf": dz_tirf,
        "dt": dt,
        "D":D,
        "R": R,
        "brightness": brightness,
        "nsteps": nsteps,
        "nparts": nparts
        }
    
    parameters_df = pd.DataFrame(parameters_dict, index=[0])
    parameters_df.to_csv(savefolder+"parameters.csv")
    
    print('---- Processing FCS acquisition-----')
    # processing
    process_stack(stack_name, first_n = 0,
                           last_n = 0,nsums = [2,3], default_dt = dt, 
                           default_psize = psize)
    
    # export 
    intensity_threshold = 0.8
    thr = 0.03
    merge_fcs_results([stack_name[:-4]+".h5"], savefolder+"FCS_results", 
          intensity_threshold = intensity_threshold, chi_threshold = thr)
    
nparts_ref = nparts
ds_to_test = np.linspace(0.01,15,20)

for j in range(len(ds_to_test)):
    R = 1
    nparts = min(10000,int(nparts_ref * (R/8)**2))
    nparts=max(50,nparts)
    simulate_spherical_diffusion(R,ds_to_test[j],nsteps,nparts)
