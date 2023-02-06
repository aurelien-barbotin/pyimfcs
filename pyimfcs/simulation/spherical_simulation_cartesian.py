# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:09:25 2021
@author: abarbotin

Coordinates definition: 
    https://en.wikipedia.org/wiki/Spherical_coordinate_system
Motion on spherical coordinates:
    https://www.e-education.psu.edu/meteo300/node/731
Uniform distribution on a sphere:
    https://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php
    
coordinates are stored in order (phi, theta)
"""
import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import gaussian

from scipy.stats import linregress
import datetime
import os
from numpy.linalg import norm

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

def g2d(x0,y0,sigma):
    y,x=coords
    return np.exp(-( (x-x0)**2 + (y-y0)**2)/(2*sigma**2))

def coord2counts(x,y,z):
    if remove_z_dependency:
        frame = g2d(x,y,sigma_psf/psize)
        return np.random.poisson(frame* brightness*dt)
        
    zr = np.sqrt(2*np.log(2))*sigmaz/psize
    omegaz = sigma_psf/psize*np.sqrt(1+z/zr**2)
    frame = g2d(x,y,omegaz)/(omegaz**2*np.pi/2)*np.exp(-z*psize/dz_tirf)

    return np.random.poisson(frame* brightness*dt)

def move_spherical(p0,mv):
    """mv: move, raw, amplitude sqrt(4dt)/r"""
    theta = p0[:,1].reshape(-1,1)
    c0 = 1/np.sqrt(1+np.pi**2*np.sin(theta)**4)
    # print(c0**2+c0**2*(np.pi*np.sin(theta))**2*np.sin(theta)**2)
    # print(mv.shape,p0.shape,c0.shape)
    mv_full= np.concatenate(( mv[:,0].reshape(-1,1)*c0, 
                             mv[:,1].reshape(-1,1)*c0*np.pi*np.sin(theta)),axis=1 )
    return (p0+mv_full)

def get_deltas(phi,theta,ampl,angle):
    denominator = np.sqrt(np.tan(angle)**2+1)
    dtheta = np.sqrt(ampl)*np.tan(angle)/denominator
    dphi = np.sqrt(ampl)/(np.sin(theta)*denominator)*np.sign(np.cos(angle))
    return dphi, dtheta

def move_spherical_upg(p0,ampl,angle):
    """mv: move, raw, amplitude sqrt(4dt)/r"""
    theta = p0[:,1]
    phi = p0[:,0]
    dphi,dtheta = get_deltas(phi,theta,ampl,angle)
    # print(c0**2+c0**2*(np.pi*np.sin(theta))**2*np.sin(theta)**2)
    # print(mv.shape,p0.shape,c0.shape)
    mv_full= np.concatenate(( 
        (phi+dphi).reshape(-1,1),
        (theta+dtheta).reshape(-1,1),
        ),axis=1 )
    return mv_full

# ---- processing helpers ------
from pyimfcs.class_imFCS import StackFCS
from pyimfcs.fitting import Fitter
from pyimfcs.io import merge_fcs_results

import tifffile
import pandas as pd
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
    # stack.set_bleaching_function(bleaching_correct_sliding,wsize = 5000)
    
    for nSum in nsums:
        stack.correlate_stack(nSum)
    if fitter is None:
        sigmaxy = sigma_psf
        print('sigmaxy',sigmaxy)
        parameters_dict = {"a":yscale, "sigma":sigmaxy}
        ft = Fitter("2D",parameters_dict, ginf=True)
    else:
        ft = fitter
    
    stack.fit_curves(ft,xmax=None)
    
    stack.save()
    

def simulate_spherical_diffusion(R,D,nsteps,nparts, 
                 savepath = "/home/aurelienb/Data/simulations/", plot=False,
                 return_coordinates=False,save=True):
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
    stack = np.zeros((nsteps-1,npix_img*2+1, npix_img*2+1))
    pos0 = np.random.uniform(size = (nparts,2))
    # phi
    pos0[:,0] = pos0[:,0]*2*np.pi
    # theta
    pos0[:,1] = np.arccos(2*pos0[:,1]-1)
    # ---------- Calculation of positions-------------------
    
    x0, y0, z0=spherical2cart(R,pos0[:,0],pos0[:,1])
    xyz = np.concatenate((x0.reshape(-1,1),
                              y0.reshape(-1,1),
                              z0.reshape(-1,1)),axis=1)
    
    if return_coordinates:
        out_coords = np.zeros((nsteps,nparts,3))
        out_coords[0] = xyz
    for j in range(1,nsteps):
        if j%500==0:
            print("Processing frame {}".format(j))
            
        bm = np.random.normal(scale=np.sqrt(2*D*dt)/R,size=(nparts,3))
        # norm xyz is R
        d_xyz = np.cross(xyz,bm)
        xyz_new = xyz+d_xyz
        xyz_new = R*xyz_new/norm(xyz_new,axis=1)[:,np.newaxis]
        xyz=xyz_new.copy()
        x1, y1, z1 = xyz_new[:,0],xyz_new[:,1],xyz_new[:,2].copy() # to not contaminate z
        z1+=R
        # pixel coordinates
        positions_new = np.array([x1,y1]).T/psize + npix_img
        
        msk0 = np.logical_and(positions_new>=0,positions_new<npix_img*2+1).all(axis=1)
        msk3 = np.logical_and(msk0,z1<z_cutoff)
        positions_new = positions_new[msk3,:]
        znew = z1[msk3]/psize
        for k in range(positions_new.shape[0]):
            frame = coord2counts(positions_new[k,0], positions_new[k,1],znew[k])
            stack[j-1]+=frame
            
        if return_coordinates:
            out_coords[j] = xyz_new
    if save:
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
                               last_n = 0,nsums = [1,2,3], default_dt = dt, 
                               default_psize = psize)
        
        # export 
        intensity_threshold = 0.8
        
        if remove_z_dependency:
            intensity_threshold = 0.5
            
        thr = 0.03
        merge_fcs_results([stack_name[:-4]+".h5"], savefolder+"FCS_results", 
              intensity_threshold = intensity_threshold, chi_threshold = thr)

    if return_coordinates:
        return out_coords
    
plot = True
save = True

psize = 0.16
sigma_psf = 0.2
sigmaz = 4*sigma_psf
dz_tirf = 0.2 # um
dt = 1*10**-3 # s
D = 1 #um2/s

R = 5 #um
brightness = 18*10**3 #Hz/molecule

remove_z_dependency = False
nsteps = 20000
nparts = 15000

nparts = 10000
npix_img = 16

coords = np.meshgrid(np.arange(2*npix_img+1),np.arange(2*npix_img+1))
z_cutoff = 3*dz_tirf

if __name__=='__main__':
    for R in [0.25,0.5,1,2]:
        # z_cutoff = R 
        # print('!!! Beware of z cutoff')
        nparts_new = int(nparts*(R/5)**2)
        simulate_spherical_diffusion(R,D,nsteps,nparts_new,
                     savepath="/home/aurelienb/Data/simulations/D1_sigma0.2_psize0p16/")
        