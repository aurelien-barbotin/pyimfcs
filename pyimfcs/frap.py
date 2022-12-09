#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 09:59:56 2022

@author: aurelienb
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
# bessel function of the first kind
from scipy.special import i1,i0

def make_model(f0,w):
    def f(t,C,D):
        vv=(w**2)/(2*D*t)
        return f0+C*np.exp(-1*vv)*(i0(vv)+i1(vv))
        
    return f

def make_model2(w):
    def f(t,C,D,f0):
        vv=(w**2)/(2*D*t)
        return f0+C*np.exp(-1*vv)*(i0(vv)+i1(vv))
        
    return f

def make_model3(w,C,f0):
    def f(t,D):
        vv=(w**2)/(2*D*t)
        if np.any(np.isnan(np.exp(-1*vv))):
            print("in exp")
        if np.any(np.isnan((i0(vv)+i1(vv)))):
            print("in i0")
        if np.any(np.isnan( C*np.exp(-1*vv)*(i0(vv)+i1(vv)) )):
            print("all")
            print(C*np.exp(-1*vv)*(i0(vv)+i1(vv)))
            print(i0(vv)+i1(vv))
            print(np.exp(-1*vv))
            
        return f0+C*np.exp(-1*vv)*(i0(vv)+i1(vv))
        
    return f

def model_exp(t,tau,f0,fmax):
    return f0+(fmax-f0)*(1-np.exp(-t*np.log(2)/(tau)))

if __name__=="__main__":
    from pyimfcs.constants import datapath
    from skimage.filters import threshold_otsu
    from scipy.optimize import curve_fit
    import tifffile
    import glob
    
    # D = 0.224 w2/t1/2
    f0=0
    C=1
    w=5 #µm
    D=10 #µm²/s
    model = make_model(f0,w)
    dts = np.linspace(0.1,60,500) #s
    intensity=model(dts,C,D)
    
    popt2,_ = curve_fit(model_exp,dts,intensity, p0=(1.12,intensity[0],intensity[-5:].mean()))
    yh2=model_exp(dts,*popt2)
    
    plt.figure()
    plt.plot(dts,intensity,label='model')
    plt.plot(dts,yh2,'k--',label='Exponential fit')
    
    vmax = intensity.max()
    index_half = np.argmin(np.abs(intensity-vmax/2))
    thalf = dts[index_half]
    plt.axhline(vmax,color="red")
    plt.axhline(vmax/2,color="black",label="half time")
    plt.axvline(thalf,color='k')
    plt.axvline(popt2[0],color='red',label="Half time from exp")
    plt.legend()
    thalf_theory = 0.224*w**2/D
    print('t half: theory {} - actual {}'.format(thalf_theory, thalf))
    print("ratio {:.2f}".format(thalf/thalf_theory))
    

    
    path = datapath+"../../NIKON/Aurelien/2022_12_05/FRAP_POPC/"
    files = glob.glob(path+"*.tif")
    file = files[2]
    stack = tifffile.imread(file)
    tmask = 72
    tstart=74
    # tstart=11
    len_recovery = 50
    mask = stack[tmask]>threshold_otsu(stack[tmask])
    
    # from skimage.morphology import binary_dilation
    # mask=binary_dilation(mask,footprint=np.ones((40,40)))
    trace = stack[tstart:tstart+len_recovery,mask].mean(axis=1)
    dt = 0.24 #s
    dts = np.arange(trace.size)*dt #s
    
    w = np.sqrt(np.count_nonzero(mask)/np.pi) *0.064 # um
    f0 = trace[0]
    C = trace[-5:].mean() - f0
    model = make_model3(w,C, f0)
    popt,_ = curve_fit(model,dts,trace, p0=(2))
    
    popt2,_ = curve_fit(model_exp,dts,trace, p0=(1.12,f0,trace[-5:].mean()))
    yh=model(dts,*popt)
    
    yh2=model_exp(dts,*popt2)
    
    plt.figure()
    plt.subplot(131)
    plt.imshow(stack[tmask])
    plt.subplot(132)
    plt.imshow(mask)
    plt.subplot(133)
    plt.plot(dts,trace,label="experimental")
    plt.plot(dts,yh,color="red",linestyle='--', label="Model fit")
    plt.plot(dts,yh2,color="k",linestyle='--', label="Exp fit")
    plt.axvline(popt2[0],color="k",label="Half time from exp")
    plt.legend()
    print('Measured D {}'.format(popt[0]))
    w_exponential=4.8
    print('Measured D exponential fit {}'.format(0.224*w**2/popt2[0]))
    # t1/2 = 0.47 with dt = 0.1 ms
    D0 = 0.224*w**2/(0.47*dt/0.1)
    #time series 4: D = 1.91, 2.09 