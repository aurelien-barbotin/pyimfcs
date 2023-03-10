#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:50:11 2022

@author: aurelien
"""
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.signal import fftconvolve

def cortrace(tr,fi):
    f0 = fi[0]
    new_tr = (tr)/np.sqrt(fi/f0)+f0*(1-np.sqrt(fi/f0))
    return new_tr

def bleaching_correct_exp(trace, plot = False):
    def expf(x,f0,tau):
        return f0*np.exp(-x/tau)
    subtr = trace
    xt = np.arange(subtr.size)
    bounds =((trace[0:trace.size//20].mean()*0.8,10),
             (trace[0:trace.size//20].mean()*1.2,trace.size*20))
    
    popt,_ = curve_fit(expf,xt,subtr,bounds = bounds)
    #popt=[4.7*10**6,18]
    trf = expf(xt,*popt)

    def cortrace(tr,popt):
        xt = np.arange(tr.size)
        fi = expf(xt,*popt)
        f0 = popt[0]
        new_tr = tr/np.sqrt(fi/f0)+f0*(1-np.sqrt(fi/f0))
        return new_tr
    new_tr = cortrace(trace,popt)
    if plot:
        plt.figure()
        plt.plot(xt,subtr)
        plt.plot(xt,trf,color='k',linestyle='--')
        plt.plot(xt,new_tr)
    #!!! Not clean
    return new_tr

def blexp_offset(trace, plot=False):
    #!!! Super artisanal
    def expf(x,f0,tau,c):
        return f0*np.exp(-x/tau)+c
    
    subtr = trace/trace.max()
    xt = np.arange(subtr.size)
    
    popt,_ = curve_fit(expf,xt,subtr)
    #popt=[4.7*10**6,18]
    trf = expf(xt,*popt)
    
    new_tr = cortrace(subtr,trf)*trace.max()
    
    if plot:
        plt.figure()
        plt.plot(subtr)
        plt.plot(trf,color="k",linestyle="--")
        plt.plot(new_tr)
    return new_tr

def blexp_double_offset(trace, plot=False, downsample = True, ndowns = 500):
    def expf(x,f0,tau,b,tau2,c):
        return f0*((1-b)*np.exp(-x/tau)+b*np.exp(-x/tau2))+c
    
    subtr = trace/trace.max()
    if downsample:
        nn0 = subtr.size//ndowns # take one point every nn0
        subtr = subtr[:nn0*ndowns]
        subtr = subtr.reshape(-1, ndowns).mean(axis=1)
        
        xt = (np.arange(subtr.size)+0.5)*ndowns
    
    p0 = (subtr.max(), #f0
          xt.max()//10, #tau
          0.5, #b
          xt.max()//2, #tau2
          subtr.min() #c
          )
    bounds = ((0, 1, 0, 1, 0),
              (2*subtr.max(), xt.max()*10, 1, xt.max()*10,1))
    
    try:
        popt,_ = curve_fit(expf,xt,subtr,p0=p0,bounds = bounds)
    except Exception as e:
        print("Bleaching correction failed with error:")
        print(e)
        return trace
    # redefine subtr and xt to have right number of points
    subtr = trace/trace.max()
    xt = np.arange(subtr.size)
    trf = expf(xt,*popt)
    # new_tr = cortrace(subtr,trf-popt[-1])+popt[-1]
    new_tr = cortrace(subtr,trf)
    
    if plot:
        plt.figure()
        plt.plot(subtr)
        plt.plot(trf,color="k",linestyle="--")
        plt.plot(new_tr)
        plt.axhline(np.mean(new_tr),linestyle='--',color='red')
    return new_tr*trace.max()

def bleaching_correct_sliding(trace, wsize = 5000, plot = False):
    u = trace.size//wsize
    trace = trace[:u*wsize]
    new_trace = trace.copy()
    m1 = trace[:wsize].mean()
    corrf = np.ones_like(trace).astype(float)
    for j in range(1,u):
        cf = m1/trace[j*wsize:j*wsize+wsize].mean()
        new_trace[j*wsize:j*wsize+wsize] = trace[j*wsize:j*wsize+wsize]*cf
        corrf[j*wsize:j*wsize+wsize] = cf
    
    if plot:
        plt.figure()
        plt.plot(trace)
        plt.plot(new_trace)
    return corrf*trace

def bleaching_correct_segment(trace, plot = False, wsize = 5000):
    if wsize%2==0:
        wsize+=1
    tr2 = np.pad(trace,wsize//2,mode='reflect', reflect_type='even')
    kernel = np.ones(wsize)/wsize
    # much faster
    conv = fftconvolve(tr2,kernel, mode='valid')
    new_trace = cortrace(trace,conv)
    if plot:
        plt.figure()
        plt.subplot(121)
        plt.plot(trace,label="raw")
        plt.plot(conv,'k--')
        
        ntr = new_trace+trace.max()-new_trace.min()
        ntr2 = np.pad(ntr,wsize//2,mode='symmetric')
        nconv = np.convolve(ntr2,kernel, mode='valid')
        plt.plot(ntr,label="corrected")
        plt.plot(nconv,'k--')
        plt.legend()
        plt.title('Traces')
        
        plt.subplot(122)
        plt.title('Residuals')
        r1 = trace-conv
        plt.plot(r1,label="Residuals from fit")
        r2 = ntr - nconv
        plt.plot(r2+r1.max()-r2.min())
        
    return new_trace

if __name__=="__main__":
    from pyimfcs.constants import datapath
    from tifffile import imread
    path = datapath+"/2022_08_05/2_withBA/Image 12.tif"
    stack = imread(path)
    x, y = 15, 38
    trace = stack[:,x,y]
    plt.figure()
    plt.subplot(121)
    plt.imshow(stack.mean(axis=0))
    plt.subplot(122)
    plt.plot(trace)