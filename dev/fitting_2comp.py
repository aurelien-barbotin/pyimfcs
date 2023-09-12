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

from pyimfcs.constants import datapath
plt.close('all')

path="/home/aurelienb/Data/simulations/SLB/2023_09_12/"

path=datapath+"/2023_09_06/S2_14h15/deltaMFD/Image 86.tif"

sigma = 0.19
a=0.08

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


    
    return sorted(popt[1:3])

files = glob.glob(path+"*/*.tif")

stacks = np.array([tifffile.imread(f) for f in files])

new_stack = stacks.sum(axis=0)[2:-2,2:-2]
print(new_stack.shape)

xc,yc = 3,5
xc,yc = new_stack.shape[1:]
all_ds=[]
for i in range(xc):
    print('correlating line {}'.format(i))
    for j in range(yc):
        trace = new_stack[:,i , j]
        corr = multipletau.autocorrelate(trace,m=32,normalize=True,deltat=10**-3)[1:]
        ds = fit_corr_2comp(corr)
        all_ds.append(ds)
all_ds = np.array(all_ds)
print('D_low = {}, D high = {}'.format( all_ds[:,0].mean(),all_ds[:,1].mean() ) )

"""
plt.figure()
plt.subplot(121)
plt.semilogx(corr[:,0], corr[:,1])

plt.semilogx(x,yh,'k--')
# plt.semilogx(x,yh2,'r--')

plt.subplot(122)
for d2 in [1,0.5,0.1,0.01]:
    
    yh3 = G_im(x,0.5,5,d2,0.5)
    plt.semilogx(x,yh3,label=d2)
plt.legend()
"""
