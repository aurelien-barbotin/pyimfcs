#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 09:33:09 2021

@author: aurelien
"""
import tifffile
from skimage.registration import phase_cross_correlation
from scipy.interpolate import interp1d

from skimage.registration._phase_cross_correlation import _upsampled_dft
from scipy.ndimage import fourier_shift

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def get_shifts(stack,ns, plot=False):
    """
    Parameters:
        stack (ndarray): 3D, txy stack.
        ns (int): frame pooling used for registration, e.g 1000
    Returns:
        list: frame arrays, shift arrays"""
    ref = np.mean(stack[:ns,],axis=0)
    
    nt = stack.shape[0]
    nstep = nt//ns
    shifts =list()
    xx = list()
    for j in range(1,nstep):
        img = np.mean(stack[j*ns:(j+1)*ns], axis=0)
        shift, error, diffphase = phase_cross_correlation(ref**2, img**2,
                                                          upsample_factor=100)
        xx.append(j*ns+ns/2)
        shifts.append(shift)
    if nstep*ns<nt:
        
        img = np.mean(stack[-ns:], axis=0)
        shift, error, diffphase = phase_cross_correlation(ref**2, img**2,
                                                          upsample_factor=100)
        xx.append(nt-ns//2)
        shifts.append(shift)
        
    shifts = np.asarray(shifts)
    xx = np.asarray(xx)
    if plot:
        plt.figure()
        plt.subplot(121)
        plt.imshow(ref)
        plt.subplot(122)
        plt.imshow(img)
        plt.suptitle("shift measurement")
    return xx, shifts

def polfit(x,a,b,c):
    return a*x**2+b*x+c*x**3

def registration(stack,ns, plot = False, method='interpolation'):
    new_stack, shifts=  stackreg(stack,ns, plot = plot, method=method)
    return new_stack


def stackreg(stack,ns, plot = False, method='interpolation'):
    """Similar to the registration method,except it returns the value of the
    x and y shifts.
    Parameters:
        stack (ndarray): the data. Time should be in the first of 3 dimensions
        ns (int): stack is summed ns by ns frames at a time to increase signal
        plot (bool): if True plots the results of the registration in a separate 
            window
        method (str): interpolation method to recover the shift at every time 
            point instead of every ns.
    Returns:
        tuple: (ndarray, ndarray); (new_stack, shifts)"""
    methods = ['interpolation', 'polfit']
    if method not in methods:
        raise KeyError('Wrong method')
        
    xx, shifts = get_shifts(stack, ns)
    shx, shy = shifts[:,0], shifts[:,1]
    
    xx2 = np.arange(stack.shape[0])
    if method=='polfit':
        popt1,_ = curve_fit(polfit, xx, shx)
        popt2,_ = curve_fit(polfit, xx, shy)
        
        xh = polfit(xx2,*popt1)
        yh = polfit(xx2,*popt2)
        
    elif method=="interpolation":
        shx = np.concatenate(([0],shx))
        shy = np.concatenate(([0],shy))
        xx = np.concatenate(([0],xx))
        fx = interp1d(xx, shx, fill_value = (shx[0],shx[-1]), bounds_error=False)
        fy = interp1d(xx, shy, fill_value = (shy[0],shy[-1]), bounds_error=False)
        
        xh = fx(xx2)
        yh = fy(xx2)

    new_stack = np.zeros_like(stack)
    new_stack[0] = stack[0]
    for j in range(1,stack.shape[0]):
        shift = (xh[j],yh[j])
        offset_image = stack[j]
        corrected_image = fourier_shift(np.fft.fftn(offset_image), shift)
        corrected_image = np.abs(np.fft.ifftn(corrected_image))
        new_stack[j] = corrected_image
        
    if plot:
        plt.figure()
        plt.subplot(221)
        plt.plot(xx,shx,label="shift x",marker="o")
        plt.plot(xx,shy, label="shift y",marker="v")
        plt.plot(xx2,xh,color="k",linestyle="--")
        plt.plot(xx2,yh,color="k",linestyle="--")
        plt.legend()
        plt.xlabel('Frame nr')
        plt.ylabel('shift (pixels')
        
        plt.subplot(222)
        plt.imshow(np.mean(stack[:ns,],axis=0))
        plt.title('First frame')
        plt.axhline(stack.shape[1]//2,color="red",linestyle='--')
        plt.axvline(stack.shape[2]//2,color="red",linestyle='--')
        
        plt.subplot(223)
        plt.imshow(np.mean(stack[-ns:,],axis=0))
        plt.title('last frame')
        plt.axhline(stack.shape[1]//2,color="red",linestyle='--')
        plt.axvline(stack.shape[2]//2,color="red",linestyle='--')
        plt.tight_layout()
        
        plt.subplot(224)
        plt.imshow(np.mean(new_stack[-ns:,],axis=0))
        plt.title('last frame corrected')
        plt.axhline(stack.shape[1]//2,color="red",linestyle='--')
        plt.axvline(stack.shape[2]//2,color="red",linestyle='--')
        plt.tight_layout()
        
    return new_stack, np.array([xh,yh])

def correct_local_dips(st, plot=True):
    
    stm = st.mean(axis=(1,2))
    stdiff = np.diff(stm[::-1])[::-1]
    
    threshold = np.percentile(stdiff,99)*3
    xxs = np.arange(stdiff.size)[stdiff>threshold]+1
    
    fig, axes = plt.subplots(1,2,sharex=True)
    
    axes[0].plot(stdiff)
    axes[0].plot(xxs-1,stdiff[xxs-1],color="red",marker="x",linestyle='')
    axes[0].set_title('Localising dips')
    
    axes[1].plot(stm)
    axes[1].plot(xxs,stm[xxs],color="red",marker="x",linestyle='')
    axes[1].set_title("Dips in intensity trace")
    
    st[xxs] = st[xxs]* (stm[xxs-1]/stm[xxs])[:,None,None]
    return st