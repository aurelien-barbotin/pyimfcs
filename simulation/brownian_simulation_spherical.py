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
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import gaussian
import tifffile
import datetime

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

def move_spherical(p0,mv):
    """mv: move, raw, amplitude sqrt(4dt)/r"""
    phi,theta = p0
    c0 = 1/np.sqrt(1+np.pi**2*np.sin(theta)**4)
    mv_full= np.array([mv[0]*c0, mv[1]*c0*np.pi*np.sin(theta)])
    return p0+mv_full
plot = True
save = True

psize = 0.16
sigma_psf = 0.2/psize
sigmaz = 4*sigma_psf
dz_tirf = 0.2 # um

dt = 1*10**-3 # s
D = 3 #um2/s

R = 5 #um
brightness = 18*10**3 #Hz/molecule

nsteps = 10000
nparts = 500

pos0 = np.random.uniform(size = (nparts,2))
# phi
pos0[:,0] = pos0[:,0]*2*np.pi
# theta
pos0[:,1] = np.arccos(2*pos0[:,1]-1)
# pos0[:,0] = pos0[:,0]/np.sin(pos0[:,1])
# ---------- Calculation of positions-------------------
moves = np.random.normal(scale = np.sqrt(2*D*dt)/R,
                         size = (nsteps,nparts,2))

moves[0] = 0

# positions = np.cumsum(moves,axis=0)
# positions = positions+pos0[np.newaxis,:,:]
positions = np.zeros_like(moves)
#thetas = moves[:,:,1]/np.sin(positions[:,:,0])
p1 = pos0[:,1]
thetas = moves[:,:,1]
positions[:,:,1] = p1[np.newaxis,:] + np.cumsum(thetas,axis=0)

p0 = pos0[:,0]
phis = moves[:,:,0]/np.sin(positions[:,:,1])
positions[:,:,0] = p0[np.newaxis,:] + np.cumsum(phis,axis=0)

positions = positions%(2*np.pi)
x,y,z = spherical2cart(R, positions[:,:,0], positions[:,:,1])

# cooper what are you doing?
# ---- plotting ----------
plt.figure()
ax = plt.axes(projection='3d')
# set_axes_equal(ax)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.scatter3D(x[-1], y[-1], z[-1],color="C0")
ax.scatter3D(0, 0, -R,color="C2")
1/0
for j in range(20):
    # ax.scatter3D(x[:,j], y[:,j], z[:,j],color="C0")
    ax.plot3D(z[:,j],y[:,j], x[:,j])

#------------ MSD stuff ---------------
def chord2arc(d,R):
    """d is chord size, R is radius"""
    return 2*R*np.arcsin(d/(2*R))

def get_chord(x,y,z,n0,n1):
    return np.sqrt( (x[n0]-x[n1])**2+(y[n0]-y[n1])**2+(z[n0]-z[n1])**2 )

def get_msds(ts,x,y,z):
    msds = list()
    for t in ts:
        chord = get_chord(x,y,z,0,t)
        msds.append( chord2arc(chord,R)**2 )
    return np.asarray(msds)
def calculate_msds(tmin=2,tmax=100,npts=10):
   
    all_mms = list()
    ts = np.linspace(tmin,tmax,npts).astype(int)
    ts= np.array(sorted(np.unique(ts)))
    
    plt.figure()
    plt.subplot(121)
    for j in range(nparts):
        
        mm = get_msds(ts,x[:,j],y[:,j],z[:,j])
        all_mms.append(mm)
        plt.plot(ts*dt,mm)
    all_mms = np.asarray(all_mms)
    mean_msd = all_mms.mean(axis=0)
    plt.plot(ts*dt,mean_msd,"k--")
    
    plt.subplot(122)
    from scipy.stats import linregress
    
    lr = linregress(ts*dt, mean_msd)
    
    plt.plot(ts*dt,mean_msd,label='Mean msd')
    plt.plot(ts*dt,lr.slope*ts*dt+lr.intercept,"k--", label="Fit")
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (um2/s)')
    plt.legend()
    
    print('MSD: {:.2f} um²/s'.format(lr.slope))
    print('D: {:.2f} um²/s'.format(lr.slope/4))

calculate_msds(tmin=2,tmax=110)

from scipy.stats import linregress
def calculate_msds_single(xx,yy,zz,tmin=2,tmax=100,npts=10,plot=False):
    # calculates msds for x,y,z set of coordinates. Returns diffusion coeff
    all_mms = list()
    ts = np.linspace(tmin,tmax,npts).astype(int)
    ts= np.array(sorted(np.unique(ts)))
    
    mm = get_msds(ts,xx,yy,zz)
    
    lr = linregress(ts*dt, mm)
    if plot:
        plt.figure()
        plt.plot(ts*dt,mm,label='Mean msd')
        plt.plot(ts*dt,lr.slope*ts*dt+lr.intercept,"k--", label="Fit")
        plt.xlabel('Time (s)')
        plt.ylabel('MSD (um2/s)')
        plt.legend()
        
    print('MSD: {:.2f} um²/s'.format(lr.slope))
    print('D: {:.2f} um²/s'.format(lr.slope/4))
    return lr.slope/4
    
"""def squaredispl(xx,yy,zz):
    return (xx[-1]-xx[0])**2+(yy[-1]-yy[0])**2+(zz[-1]-zz[0])"""
calculate_msds_single(x[:,0],y[:,0],z[:,0],plot=True)
all_ds = []
for j in range(nparts):
    xx,yy,zz=x[:,j],y[:,j],z[:,j]
    d0 = calculate_msds_single(xx,yy,zz,plot=False)
    all_ds.append(d0)

plt.figure()
plt.subplot(221)
plt.scatter(pos0[:,0],all_ds)
plt.xlabel('Initial phi')
plt.ylabel('D [µm²/s]')
plt.subplot(222)
plt.scatter(pos0[:,1],all_ds)
plt.xlabel('Initial theta')
plt.ylabel('D [µm²/s]')

plt.subplot(223)
plt.hist2d(pos0[:,0],all_ds,bins=20)
plt.xlabel('Initial phi')
plt.ylabel('D [µm²/s]')
plt.title('Phi')
plt.subplot(224)
hh=plt.hist2d(pos0[:,1],all_ds)
plt.xlabel('Initial theta')
plt.ylabel('D [µm²/s]')
print(hh[0].sum(axis=1))
all_ds=np.asarray(all_ds)
thetas_up = np.linspace(0,np.pi,10)
dmeans=[]
dstd=[]
for j in range(1,len(thetas_up)):
    msk = np.logical_and(pos0[:,1]<thetas_up[j],pos0[:,1]>thetas_up[j-1])
    ds = all_ds[msk]
    dmeans.append(np.mean(ds))
    dstd.append(np.std(ds))

plt.figure()
plt.errorbar(thetas_up[1:],dmeans,yerr=dstd,capsize=5)
plt.axhline(D,color='k',linestyle='--')
plt.xlabel('Initial theta value')
plt.ylabel('D [µm²/s]')

calculate_msds()

# density
coord_to_check=z[-1:]
subareas = np.linspace(-R,R,50)
nrs = []
for j in range(1,len(subareas)):
    msk = np.logical_and(coord_to_check>subareas[j-1], coord_to_check<=subareas[j])
    nr = np.count_nonzero(msk)
    nrs.append(nr)
nrs = np.asarray(nrs).astype(float)
nrs/=float(coord_to_check.shape[0])
# (np.pi*R**2*(subareas[1]-subareas[0]))
theoretical_density = nparts/subareas.size

plt.figure()
plt.plot(subareas[1:], nrs)
plt.axhline(theoretical_density,color='k',linestyle='--')
plt.xlabel('distance [µm]')
plt.ylabel('Particle density')
