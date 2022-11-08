#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:20:25 2021

@author: aurelien
"""
import numpy as np
from scipy.special import erf

from pyimfcs.fitting import gim2D
import matplotlib.pyplot as plt

datapath = "/run/user/1000/gvfs/smb-share:server=data.micalis.com,share=proced/microscopy/ZEISS/Aurelien/"

def calculate_a_eff(a,w0):
    """Calculates effective observation area in imFCS from parameter a (pixel size)
    and w., 1/e2 radius of Gaussian PSF
    
    From (1) Sankaran, J.; Bag, N.; Kraut, R. S.; Wohland, T. Accuracy and 
    Precision in Camera-Based Fluorescence Correlation Spectroscopy Measurements. 
    Anal. Chem. 2013, 85 (8), 3948–3954. https://doi.org/10.1021/ac303485t.
    
    """
    a_eff = a**2/(erf(a/w0)+ (np.exp(-a**2/w0**2-1))*w0/(a*np.sqrt(np.pi)))
    return a_eff

def calculate_tauD(aeff,D):
    """From same ref as above, get diffusion coeff"""
    return aeff/(4*D)

def tauD_fromfit(D,a,w0):
    """From fitting parameters, callculates transit time"""
    aeff = calculate_a_eff(a,w0)
    return calculate_tauD(aeff,D)

def calculate_dmax(exp_time,aeff):
    tau1 = 10*exp_time
    Dmax = aeff/(4*tau1)
    return Dmax

if __name__=='__main__':
    plt.close('all')
    bins = np.arange(1,4)
    a_effs = bins*0.16
    w0 = 0.2 #um
    exptime = 4.76*10**-3
    a_effs = calculate_a_eff(a_effs,w0)
    Dinterest = 5
    
    taus=calculate_tauD(a_effs, Dinterest)*10**3
    dmax = calculate_dmax(exptime,a_effs)
    
    for j in range(len(dmax)):
        print("binning {}, Dmax {:.2f}".format(bins[j], dmax[j]))
        print("D {} µm2/s, binning {}, tauD {:.2f}".format(Dinterest,bins[j], taus[j]))
    
    g2 = gim2D(a=0.16*3,sigma=0.2)
    taus = np.logspace(-3,1)*1.26
    ds = [1,5,7,10]
    
    plt.figure()
    for dd in ds:
        y = g2(taus,1,dd)
        plt.semilogx(taus,y/y.max())