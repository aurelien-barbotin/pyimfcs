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

from scipy.stats import linregress
import datetime
import os
from numpy.linalg import norm

from pyimfcs.class_imFCS import StackFCS
from pyimfcs.fitting import Fitter
from pyimfcs.export import merge_fcs_results

import tifffile
import pandas as pd
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
    return np.exp(-( (x-x0)**2 + (y-y0)**2)/(2*sigma**2))/sigma**2

def coord2counts(x,y,z):
    # pixel coordinates
    frame = g2d(x,y,sigma_psf/psize)*np.exp(-z*psize/dz_tirf)
    return np.random.poisson(frame* brightness*dt)


# ---- processing helpers ------
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
        parameters_dict = {"a":yscale, "sigma":sigmaxy,"mtype":"2D","ginf":True}
        ft = Fitter(parameters_dict)
    else:
        ft = fitter
    ft.bounds=((10**-4,0.03,-1),
               (10,D*3,1))
    stack.fit_curves(ft,xmax=None)
    
    stack.save()    

def simulate_spherical_diffusion(R,D,nsteps,nparts, 
                 savepath = "/home/aurelienb/Data/simulations/", plot=False,
                 return_coordinates=False, save=True):
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
    stack = np.zeros((nsteps,npix_img*2+1, npix_img*2+1))
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
        
    for j in range(1,nsteps+1):
        if j%500==0:
            print("Processing frame {}".format(j))
            
        bm = np.random.normal(scale=np.sqrt(2*D*dt)/R,size=(nparts,3))
        
        # norm xyz is R
        d_xyz = np.cross(xyz,bm)
        xyz_new = xyz+d_xyz
        xyz_new = R*xyz_new/norm(xyz_new,axis=1)[:,np.newaxis]
        
        if return_coordinates:
            out_coords[j]=xyz_new
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
            
    if plot and return_coordinates:
        if plot:
            plt.figure()
            ax = plt.axes(projection='3d')
            # set_axes_equal(ax)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            # ax.scatter3D(x[-1], y[-1], z[-1],color="C0")
            for j in range(20):
                # ax.scatter3D(x[:,j], y[:,j], z[:,j],color="C0")
                print(out_coords[:,j].shape)
                
                ax.plot3D(out_coords[:,j,0],out_coords[:,j,1], out_coords[:,j,2])
            
    if save:
        fname = datetime.datetime.now().__str__()[:19]
        savefolder = savepath+fname+"/"
        os.mkdir(savefolder)
        stack_name = savefolder+"stack.tif"
        tifffile.imwrite(stack_name,stack)
        
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
                               last_n = 0,nsums = [1,2,3,4,5], default_dt = dt, 
                               default_psize = psize)
        
        # export 
        intensity_threshold = 0.8
        
            
        thr = 0.03
        merge_fcs_results( savefolder+"FCS_results",[stack_name[:-4]+".h5"], 
              ith = intensity_threshold, chi_threshold = thr)

    if return_coordinates:
        return out_coords
    
plot = True
save = True

psize = 0.05
dz_tirf = 0.1 # um
dt = 1*10**-3 # s
D = 1 #um2/s

sigma_psf = 0.19
sigmaz = 4*sigma_psf

R = 10 #um, RADIUS!
brightness = 18*10**3 #Hz/molecule

nsteps = 50000
nparts = 500
npix_img = 10

coords = np.meshgrid(np.arange(2*npix_img+1),np.arange(2*npix_img+1))
z_cutoff = 50*dz_tirf

if __name__=='__main__':
    import time
    t0 = time.time()
    # z_cutoff = R 
    # print('!!! Beware of z cutoff')
    for j in range(2):
        for ps in [0.05,0.1,0.16,0.3]:
            psize=ps
            simulate_spherical_diffusion(R,D,nsteps,nparts,plot=False,return_coordinates=False,
                 savepath="/home/aurelienb/Data/simulations/2023_06_05/psizes/ztirf100nm4/")
    print('elapsed: ',time.time()-t0)