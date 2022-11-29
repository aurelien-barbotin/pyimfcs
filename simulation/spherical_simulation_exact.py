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
https://math.stackexchange.com/questions/3725288/infinitesimal-generator-of-the-brownian-motion-on-a-sphere
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import gaussian

from scipy.stats import linregress
plt.close('all')


def spherical2cart(R,phi,theta):
    x = R*np.cos(phi)*np.sin(theta)
    y = R*np.sin(phi)*np.sin(theta)
    z = R*np.cos(theta)
    return np.array([x, y, z])

def cart2spherical(x,y,z):
    theta = np.arccos(z)
    phi = np.arctan2(y,x)
    return np.array([phi, theta])

plot = True
save = True

psize = 0.16
sigma_psf = 0.2/psize
sigmaz = 4*sigma_psf
dz_tirf = 0.2 # um

dt = 1*10**-3 # s
D = 5 #um2/s

R = 5 #um
brightness = 18*10**3 #Hz/molecule

nsteps = 1000
nparts = 1500

pos0 = np.random.uniform(size = (nparts,2))
# phi
pos0[:,0] = pos0[:,0]*2*np.pi
# theta
pos0[:,1] = np.arccos(2*pos0[:,1]-1)
# pos0[:,0] = pos0[:,0]/np.sin(pos0[:,1])
# ---------- Calculation of positions-------------------
moves = np.random.normal(scale = np.sqrt(4*D*dt)/R,
                         size = (nsteps,nparts,2) )

moves[0] = 0

# ---- new version
amplitudes = np.random.normal(scale = np.sqrt(4*D*dt)/R,
                         size = (nsteps,nparts) )**2
angles = np.random.uniform(low=0,high=2*np.pi,size=(nsteps,nparts))
# positions = np.cumsum(moves,axis=0)
# positions = positions+pos0[np.newaxis,:,:]
positions = np.zeros((nsteps,nparts,2))

moves = np.random.normal(scale = np.sqrt(2*D*dt)/R,
                         size = (nsteps,nparts,3))
positions[0] = pos0
x = np.zeros((nsteps,nparts))
y = np.zeros((nsteps,nparts))
z = np.zeros((nsteps,nparts))
x[0], y[0], z[0]=spherical2cart(R,pos0[:,0],pos0[:,1])
for j in range(1,nsteps):
    # positions[j] = move_spherical(positions[j-1],moves[j])
    p0 = positions[j-1]
    phi_t = p0[:,0]
    theta_t = p0[:,1]
    dBtheta = np.sin(phi_t)*moves[j,:,0]-np.cos(phi_t)*moves[j,:,1]
    dBphi = np.cos(theta_t)*(np.cos(phi_t)*moves[j,:,0]+np.sin(phi_t)*moves[j,:,1] )-np.sin(theta_t)*moves[j,:,2]
    dthetat = dBtheta +0.5/np.tan(theta_t)*dt/2
    dphit = dBphi/np.sin(theta_t)

    positions[j,:,0] = phi_t + dphit
    positions[j,:,1] = theta_t + dthetat
    
x,y,z = spherical2cart(R, positions[:,:,0], positions[:,:,1])

# ---- plotting ----------è
plt.figure()
ax = plt.axes(projection='3d')
# set_axes_equal(ax)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# ax.scatter3D(x[0], y[0], z[0],color="C0")
ax.scatter3D(x[-1], y[-1], z[-1],color="C1")

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
    # calculates msds all over the map
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
all_ds=np.asarray(all_ds)

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
thetas_up = np.linspace(0,np.pi,10)
dmeans=[]
dstd=[]
for j in range(1,len(thetas_up)):
    msk = np.logical_and(pos0[:,1]<thetas_up[j],pos0[:,1]>thetas_up[j-1])
    ds = all_ds[msk]
    dmeans.append(np.mean(ds))
    dstd.append(np.std(ds))
    
plt.figure()
plt.subplot(121)
plt.errorbar(thetas_up[1:],dmeans,yerr=dstd,capsize=5)
plt.axhline(D,color='k',linestyle='--')
plt.xlabel('Initial theta value')
plt.ylabel('D [µm²/s]')

plt.subplot(122)
phis_up = np.linspace(0,np.pi,10)
dmeans=[]
dstd=[]
for j in range(1,len(phis_up)):
    msk = np.logical_and(pos0[:,0]<phis_up[j],pos0[:,0]>phis_up[j-1])
    ds = all_ds[msk]
    dmeans.append(np.mean(ds))
    dstd.append(np.std(ds))
plt.errorbar(phis_up[1:],dmeans,yerr=dstd,capsize=5)
plt.axhline(D,color='k',linestyle='--')
plt.xlabel('Initial phi value')
plt.ylabel('D [µm²/s]')

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
plt.subplot(121)
plt.plot(subareas[1:], nrs)
plt.axhline(theoretical_density,color='k',linestyle='--')
plt.subplot(122)
plt.plot(subareas[1:], nrs)
plt.axhline(theoretical_density,color='k',linestyle='--')
plt.xlabel('distance [µm]')
plt.ylabel('Particle density')