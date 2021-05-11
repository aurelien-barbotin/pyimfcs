# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:09:25 2021

@author: abarbotin
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import gaussian
import tifffile

plt.close('all')

def gauss(x,sig):
    return np.exp(-x**2/(2*sig**2))

def make_image():
    pass

scale = 0.5 #pixels
# pixel size: 100 nm
psize = 100
sigma_psf = 110/psize
dt = 4.7*10**-3 # s
D =  #um2/s

brightness = 18*10**3 #Hz/molecule

npixels = 500

nsteps = 20000
nparts = 3000

pos0 = np.random.uniform(size = (nparts,2))*npixels-npixels/2
moves = np.random.normal(scale = np.sqrt(2*D*dt)/(psize*10**-3),size = (nsteps,nparts,2) )

positions = np.cumsum(moves,axis=0)
positions = positions+pos0[np.newaxis,:,:]

xm, ym = positions.min(), positions.max()
plot = True

npix_img = 32
stack = np.zeros((nsteps,npix_img*2+1, npix_img*2+1))

for j in range(nsteps):
    
    positions_new = positions[j]
    positions_new = positions_new[np.logical_and(np.abs(positions_new[:,0])<=npix_img+1,np.abs(positions_new[:,1])<=npix_img+1 )]
    positions_new = positions_new.astype(int) + npix_img
    for k in range(len(positions_new)):
        stack[j, positions_new[k, 0],positions_new[k, 1]]+=1
    stack[j] = gaussian(stack[j],sigma = sigma_psf)
    
    if j%500==0:
        print("Processing frame {}".format(j))
stack = stack * brightness*dt 
stack = np.random.poisson(stack)
x = np.linspace(-(npix_img*2+1)/2,(npix_img*2+1)/2,npix_img*2+1)

xx, yy = np.meshgrid(x,x)


plt.figure()
plt.imshow(stack[0])
if plot:
    plt.figure()
    plt.subplot(121)
    plt.plot(positions[:,50,0], positions[:,50,1])
    
    plt.subplot(122)
    for j in range(10):
        plt.plot(positions[:,j,0], positions[:,j,1])

save=True
if save:
    stack = stack.astype(np.uint16)
    tifffile.imwrite("simulation2_4p7msexp_D{:.2f}.tif".format(D),stack)
