#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 10:03:57 2022

@author: aurelienb
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
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
    
from spherical_simulation_cartesian import simulate_spherical_diffusion, cart2spherical

R=0.25
nparts=200
nsteps=1000
D=0.1
dt=10**-3
xyz = simulate_spherical_diffusion(R,D,nsteps,nparts,return_coordinates=True,save=False)
x,y,z = xyz[:,:,0],xyz[:,:,1],xyz[:,:,2]

calculate_msds_single(x[:,0],y[:,0],z[:,0],plot=True)
all_ds = []
for j in range(nparts):
    xx,yy,zz=x[:,j],y[:,j],z[:,j]
    d0 = calculate_msds_single(xx,yy,zz,plot=False)
    all_ds.append(d0)
all_ds=np.asarray(all_ds)

pos0 = cart2spherical(x[0]/R,y[0]/R,z[0]/R).T
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
plt.plot(subareas[1:], nrs)
plt.axhline(theoretical_density,color='k',linestyle='--')
plt.xlabel('z distance [µm]')
plt.ylabel('Particle density')

calculate_msds(tmax=20)
