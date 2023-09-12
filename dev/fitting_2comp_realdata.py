#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 10:11:14 2023

@author: aurelienb
"""

import numpy as np
import tifffile
import glob
import multipletau
import matplotlib.pyplot as plt
from pyimfcs.fitting import gim2D_2components
from scipy.special import erf
from scipy.optimize import curve_fit
from pyimfcs.blcorr import blexp_double_offset
from pyimfcs.constants import datapath
plt.close('all')

tmppath="/run/user/1001/gvfs/microscopy/Zeiss/Aurelien/"
path=tmppath+"/2023_09_06/S2_14h15/deltaMFD/Image 86.tif"
path="/run/user/1001/gvfs/microscopy/ZEISS/Aurelien/2023_09_06/S2_14h15/WT/Image 77.tif"
sigma = 0.19
a=0.16

def G_im(tau,N,D1,D2,A2):
    """"""
    k1 = a/(2*np.sqrt(D1*tau+sigma**2 ) )
    k2 = a/(2*np.sqrt(D2*tau+sigma**2 ) )
    return 1/N*( 
        A2*(erf(k1)+(np.exp(-k1**2)-1)/(k1*np.sqrt(np.pi)))**2 +
        (1-A2)*(erf(k2)+(np.exp(-k2**2)-1)/(k2*np.sqrt(np.pi)))**2
                )

def fit_corr_2comp(corr):
    x, y = corr[:,0], corr[:,1]
    
    bounds = ((10**-4,sigma**2/x.max(),sigma**2/x.max(),0),
              (10**4,sigma**2/x.min(),sigma**2/x.min(),1))
    popt,_ = curve_fit(G_im,x,y,bounds=bounds)
    #print('Measured D1 {:.3f}, D2 {:.3f}, theory: 5 and 0.1 (6 and 0.15)'.format(popt[2],popt[1]))
    yh = G_im(x,*popt)
    plt.figure()
    plt.semilogx(corr[:,0], corr[:,1])

    plt.semilogx(x,yh,'k--')
    # plt.semilogx(x,yh2,'r--')
    return popt

new_stack = tifffile.imread(path)
nSum=2
v,w=new_stack.shape[1:]
all_traces=[]
for i in range(v//nSum):
    ctmp = []
    trtmp = []
    for j in range(w//nSum):
        trace = new_stack[:, i * nSum:i * nSum + nSum,
                j * nSum:j * nSum + nSum].mean(axis=(1, 2))
        trtmp.append(trace)
    all_traces.append(trtmp)
    
new_stack=np.swapaxes(np.array(all_traces).T, 1, 2)

miniature = new_stack.mean(axis=0)
print(new_stack.shape)
j,i = 27,6
j,i = 22,1
plt.figure()
plt.imshow(miniature)
plt.plot(j,i,'rx')
print(miniature[i,j])
trace = new_stack[:,i , j]
trace_corrected = blexp_double_offset(trace)
corr = multipletau.autocorrelate(trace_corrected,m=32,normalize=True,deltat=10**-3)[1:]
x,y=corr[:,0],corr[:,1]

popt = fit_corr_2comp(corr)

yh=G_im(x,*popt)

plt.figure()
plt.subplot(121)
plt.semilogx(corr[:,0], corr[:,1])

plt.semilogx(x,yh,'k--')
# plt.semilogx(x,yh2,'r--')

plt.subplot(122)
plt.plot(trace_corrected)
plt.legend()

ds = sorted(popt[1:3])
print('Measured Ds: {:.2f}, {:.2f} µm²/s'.format(*ds))
