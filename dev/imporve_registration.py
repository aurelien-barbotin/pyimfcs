#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 13:39:13 2022

@author: aurelienb
"""

import tifffile
from skimage.registration import phase_cross_correlation, optical_flow_tvl1
from scipy.interpolate import interp1d

from skimage.registration._phase_cross_correlation import _upsampled_dft
from scipy.ndimage import fourier_shift

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from pyimfcs.constants import datapath
from pyimfcs.shift_correction import get_shifts
import tifffile
plt.close('all')
plot = True

# total shift measured manualy: 16, 4
path = datapath+"2021/2021_11_04/BSLB/DOPC/tmp_im6/Image 6_cropped.tif"
ns = 3000
stack = tifffile.imread(path)[15000:]
ref = np.mean(stack[:ns,],axis=0)

nt = stack.shape[0]
nstep = nt//ns
shifts =list()
global_shift = np.zeros(2)
for j in range(1,nstep):
    # ref = np.mean(stack[(j-1)*ns:j*ns],axis=0)
    img = np.mean(stack[j*ns:(j+1)*ns], axis=0)
    shift, error, diffphase = phase_cross_correlation(ref**4, img**4,
                                                      upsample_factor=100)
    ref = img
    shifts.append(shift)
    #global_shift+=shift
ref = np.mean(stack[:ns,],axis=0)
shifts = np.asarray(shifts)
shifts_sum = np.cumsum(shifts,axis=0)

if plot:
    plt.figure()
    plt.subplot(121)
    plt.imshow(ref)
    plt.subplot(122)
    plt.imshow(img)
    plt.suptitle("shift measurement")
    
plt.figure()
plt.subplot(121)
plt.plot(shifts,"-o")
plt.subplot(122)
plt.plot(shifts_sum,"-o")

j=8

img = np.mean(stack[j*ns:(j+1)*ns], axis=0)
plt.figure()
plt.subplot(121)
plt.imshow(ref**8)
plt.subplot(122)
plt.imshow(img**8)

shift = shifts[-1]
corrected_image = fourier_shift(np.fft.fftn(np.mean(stack[-ns:], axis=0)), shift)
corrected_image = np.abs(np.fft.ifftn(corrected_image))

plt.figure()
plt.subplot(131)
plt.imshow(ref)
plt.title('first image')
plt.axhline(stack.shape[1]//2,color="red",linestyle='--')
plt.axvline(stack.shape[2]//2,color="red",linestyle='--')
plt.subplot(132)
plt.imshow(stack[-ns:].mean(axis=0))
plt.title('Last image')
plt.axhline(stack.shape[1]//2,color="red",linestyle='--')
plt.axvline(stack.shape[2]//2,color="red",linestyle='--')

plt.subplot(133)
plt.imshow(corrected_image)
plt.axhline(stack.shape[1]//2,color="red",linestyle='--')
plt.axvline(stack.shape[2]//2,color="red",linestyle='--')
plt.title('Last image corrected')

shifts_ref = get_shifts(stack,ns)

lastn = np.mean(stack[-ns:], axis=0)
out=optical_flow_tvl1(ref,lastn)
plt.figure()
plt.subplot(121)
plt.imshow(out[0])
plt.subplot(122)
plt.imshow(out[1])
