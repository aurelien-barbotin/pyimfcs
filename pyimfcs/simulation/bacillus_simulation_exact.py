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

rng = np.random.default_rng()


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

def cylindrical2cart(R,theta,xcylindr):
    x = xcylindr
    y = R*np.cos(theta)
    z = R*np.sin(theta)
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
    # frame=np.random.poisson(frame* brightness*dt)
    frame=rng.poisson(lam=frame* brightness*dt,size=frame.shape)

    return frame


# ---- processing helpers ------
def process_stack(path,first_n = 0, last_n = 0, nsums=[2,3],
                           plot=False, default_dt= None, default_psize = None, 
                           fitter = None, export_summaries = True, 
                           chi_threshold = 0.03, ith=0.8,shifts=(0,0)):

    stack = StackFCS(path, background_correction = True,                     
                         first_n = first_n, last_n = last_n, clipval = 0)
    
    print(stack.stack.dtype)
    stack.dt = default_dt
    stack.xscale = default_psize
    stack.yscale = default_psize
    xscale = stack.xscale
    yscale = stack.yscale
    assert(xscale==yscale)
    if shifts is not None:
        sx,sy = shifts
        stack.stack = stack.stack[:,sx:,sy:]
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
    
    stack.save(exclude_list=['traces_dict'])    

def normal_to_bacillus(xyz):
    """Calculates a vector normal to a point in a bacillus configuration"""
    vec_u = xyz.copy()
    # positive values of x
    msk_pos = vec_u[:,0]>0
    vec_u[msk_pos,0] = np.max(
        np.concatenate(
            (vec_u[msk_pos,0].reshape(-1,1)-length/2,np.zeros(np.count_nonzero(msk_pos)).reshape(-1,1))
            ,axis=1)
        ,axis=1)
    
    msk_neg = vec_u[:,0]<=0
    vec_u[msk_neg,0] = np.min(
        np.concatenate(
            (vec_u[msk_neg,0].reshape(-1,1)+length/2,
             np.zeros(np.count_nonzero(msk_neg)).reshape(-1,1)
             )
            ,axis=1) # end concatenate
        ,axis=1)
    
    return vec_u

def simulate_spherical_diffusion(R,length,D,nsteps,nparts,
                 savepath = "/home/aurelienb/Data/simulations/", plot=False,
                 return_coordinates=False, save=True,shifts=(0,0),delete_tif=True, nsums=[2]):
    if not os.path.isdir(savepath):
        os.mkdir(savepath)
    stack = np.zeros((nsteps,npix_img*2+1, npix_img*2+1),dtype=int)
    f1 = R/(R+2*length) # fraction on sphere
    pos0 = np.random.uniform(size = (int(f1*nparts),2))
    # phi
    pos0[:,0] = pos0[:,0]*2*np.pi
    # theta
    pos0[:,1] = np.arccos(2*pos0[:,1]-1)
    
    # ---------- Calculation of positions on sphere part-------------------
    x0, y0, z0=spherical2cart(R,pos0[:,0],pos0[:,1])
    # cylinder is on x axis
    x0[x0<0]-=length/2
    x0[x0>0]+=length/2
    
    pos1 = np.random.uniform(size = (nparts-int(f1*nparts),2))
    pos1[:,0] = pos1[:,0]*2*np.pi # theta
    pos1[:,1] = pos1[:,1]*length-length/2 #x
    x1,y1,z1 = cylindrical2cart(R, pos1[:,0], pos1[:,1])
    xyz = np.concatenate((np.concatenate((x0,x1)).reshape(-1,1),
                          np.concatenate((y0,y1)).reshape(-1,1),
                          np.concatenate((z0,z1)).reshape(-1,1) ),axis=1)
    
    if plot:
        plt.figure()
        ax = plt.axes(projection='3d')
        # set_axes_equal(ax)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        # ax.scatter3D(x[-1], y[-1], z[-1],color="C0")
        # ax.plot3D(x0,y0,z0)
        # ax.scatter(x0,y0,z0)
        # ax.scatter(x1,y1,z1)
        ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2])
        set_axes_equal(ax)

        
    if return_coordinates:
        out_coords = np.zeros((nsteps+1,nparts,3))
        out_coords[0] = xyz
        
    for j in range(1,nsteps+1):
        if j%500==0:
            print("Processing frame {}".format(j))
            
        bm = np.random.normal(scale=np.sqrt(2*D*dt)/R,size=(nparts,3))

        vec_u = normal_to_bacillus(xyz) # has norm R in physical coordinates
        
        d_xyz = np.cross(vec_u,bm)
        xyz_new = xyz+d_xyz
        norm_u=norm(vec_u,axis=1)[:,np.newaxis]
        xyz_new = xyz_new + (R-norm_u)*vec_u/norm_u
        # xyz_new = R*xyz_new/norm_new[:,np.newaxis]
        
        if return_coordinates:
            out_coords[j]=xyz_new
        xyz=xyz_new.copy()
        x1, y1, z1 = xyz_new[:,0],xyz_new[:,1],xyz_new[:,2].copy() # to not contaminate z
        z1+=R
        # pixel coordinates
        positions_new = np.array([x1,y1]).T/psize
        # print(positions_new.min(),positions_new.max())
        positions_new+=npix_img
        
        msk0 = np.logical_and(positions_new>=0,positions_new<npix_img*2+1).all(axis=1)
        msk3 = np.logical_and(msk0,z1<z_cutoff)
        positions_new = positions_new[msk3,:]
        znew = z1[msk3]/psize
        for k in range(positions_new.shape[0]):
            frame = coord2counts(positions_new[k,0], positions_new[k,1],znew[k])
            stack[j-1]+=frame
            
        # print(stack.dtype)
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
                ax.plot3D(out_coords[:,j,0],out_coords[:,j,1], out_coords[:,j,2])
                set_axes_equal(ax)
            
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
            "length":length,
            "brightness": brightness,
            "nsteps": nsteps,
            "nparts": nparts
            }
        
        parameters_df = pd.DataFrame(parameters_dict, index=[0])
        parameters_df.to_csv(savefolder+"parameters.csv")
        
        print('---- Processing FCS acquisition-----')
        # processing
        process_stack(stack_name, first_n = 0,
                               last_n = 0,nsums = nsums, default_dt = dt, 
                               default_psize = psize,shifts=shifts)
        if delete_tif:
            os.remove(savefolder+"stack.tif")
        # export 
        intensity_threshold = 0.8
        
            
        thr = 0.03
        merge_fcs_results( savefolder+"FCS_results",[stack_name[:-4]+".h5"], 
              ith = intensity_threshold, chi_threshold = thr)

    if return_coordinates:
        return out_coords
    
plot = True
save = True

psize = 0.08
dz_tirf = 0.1 # um
dt = 1*10**-3 # s
D = 2 #um2/s

sigma_psf = 0.19
sigmaz = 4*sigma_psf

R = 0.5 #um, RADIUS!
length=3
brightness = 18*10**4 #Hz/molecule

nsteps = 50000
nparts = 500
npix_img = 30

coords = np.meshgrid(np.arange(2*npix_img+1),np.arange(2*npix_img+1))
z_cutoff = 10**3*dz_tirf

process_sums = [2,3,4]

if __name__=='__main__':
    import time
    t0 = time.time()
   
    simulate_spherical_diffusion(R,length,D,nsteps,nparts,plot=False,
                                 return_coordinates=False,
         savepath="/home/aurelienb/Data/simulations/2023_07_03_bacillus/",
         delete_tif=False, nsums=process_sums)


    print('elapsed: ',time.time()-t0)