#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 09:33:09 2021

@author: aurelien
"""
import tifffile
from skimage.registration import phase_cross_correlation

from skimage.registration._phase_cross_correlation import _upsampled_dft
from scipy.ndimage import fourier_shift

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
plt.close('all')

def get_shifts(stack,ns):
    """
    Parameters:
        stack (ndarray): 3D, txy stack.
        ns (int): frame pooling used for registration, e.g 1000"""
    ref = stack[:ns,].mean(axis=0)
    
    nt = stack.shape[0]
    nstep = nt//ns
    shifts  =list()
    for j in range(1,nstep):
        img = stack[j*ns:(j+1)*ns].mean(axis=0)
        shift, error, diffphase = phase_cross_correlation(ref, img,
                                                          upsample_factor=100)
        shifts.append(shift)
    shifts = np.asarray(shifts)
    return shifts

def polfit(x,a,b,c):
    return a*x**2+b*x+c*x**3

def registration(stack,ns):
    shifts = get_shifts(stack, ns)
    
    shx, shy= shifts[:,0], shifts[:,1]
    xx = np.arange(shx.size)+1
    xx*=ns
    popt1,_ = curve_fit(polfit, xx, shx)
    popt2,_ = curve_fit(polfit, xx, shy)
    
    xx2 = np.arange(stack.shape[0])
    
    xh = polfit(xx2,*popt1)
    yh = polfit(xx2,*popt2)
    plt.figure()
    plt.plot(xx,shifts[:,0],label="shift x",marker="o")
    plt.plot(xx,shifts[:,1], label="shift y",marker="v")
    plt.plot(xx2,xh,color="k",linestyle="--")
    plt.plot(xx2,yh,color="k",linestyle="--")
    plt.legend()
    plt.xlabel('Frame nr')
    plt.ylabel('shift (pixels')
    
    new_stack = np.zeros_like(stack)
    new_stack[0] = stack[0]
    for j in range(1,stack.shape[0]):
        shift = (xh[j],yh[j])
        offset_image = stack[j]
        corrected_image = fourier_shift(np.fft.fftn(offset_image), shift)
        corrected_image = np.abs(np.fft.ifftn(corrected_image))
        new_stack[j] = corrected_image
    
    return new_stack