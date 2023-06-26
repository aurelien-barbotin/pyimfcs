#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 09:34:59 2023

@author: aurelienb
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.close('all')

def g2d(x0,y0,sigma):
    y,x=coords
    return np.exp(-( (x-x0)**2 + (y-y0)**2)/(2*sigma**2))

npix_img = 24
coords = np.meshgrid(np.arange(2*npix_img+1),np.arange(2*npix_img+1))
psize=0.1
psf = g2d(npix_img,npix_img,0.2/psize)
sigma_psf=0.2/psize
brightness = 18*10**2
dt = 10**-3

real_counts = np.random.poisson(brightness*dt)
poisson_psf = np.random.poisson(psf* brightness*dt)

psf_lin = psf.reshape(-1)
psf_lin = np.cumsum(psf_lin) # cdf
print(psf_lin.max())

psf_lin/=psf_lin[-1]

fi = interp1d(psf_lin,np.arange(psf.size))

samples = np.random.uniform(low=0,high=1,size=real_counts)
coords_new = fi(samples).astype(int) # !!! bfof

import time

# 
def coord2counts_exact(x,y):
    #!!! in pixel coordinates
    t0 = time.time()
    frame = g2d(x,y,sigma_psf)
    t1 = time.time()
    real_counts = np.random.poisson(brightness*dt)
    psf_lin = frame.reshape(-1)
    psf_lin = np.cumsum(psf_lin) # cdf
    t2 = time.time()
    
    psf_lin/=psf_lin[-1]
    fi = interp1d(psf_lin,np.arange(psf.size))
    t3 = time.time()
    
    samples = np.random.uniform(low=0,high=1,size=real_counts)
    coords = fi(samples).astype(int) 
    
    out = np.zeros_like(psf_lin)
    for coord in coords:
        out[coord]+=1
        
    t4 = time.time()
    print('elapsed: ', t1-t0, t2-t1, t3-t2,t4-t3)
    return out.reshape(frame.shape)

out = np.zeros_like(psf_lin)
for coord in coords_new:
    out[coord]+=1

uu = coord2counts_exact(npix_img, npix_img)

plt.figure()
plt.subplot(221)
plt.imshow(psf)
plt.title('psf')
plt.subplot(222)
plt.plot(psf_lin)
plt.title('cdf')
plt.subplot(223)
plt.imshow(uu)
plt.title('exact sampling')
plt.subplot(224)
plt.imshow(poisson_psf)
plt.title('Naive sampling')

dp = poisson_psf/poisson_psf.max()-psf

plt.figure()
plt.subplot(221)
plt.imshow(psf)
plt.title('original')
plt.subplot(222)
plt.imshow(dp)
plt.title('difference')
plt.subplot(223)
plt.plot(psf[psf.shape[0]//2]/psf.max())
plt.plot(poisson_psf[poisson_psf.shape[0]//2]/poisson_psf.max())