#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:20:25 2021

@author: aurelien
"""
import numpy as np
from scipy.special import erf

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
if __name__=='__main__':
    a_effs = np.arange(1,4)*0.16
    w0 = 0.2 #um
    a_effs = calculate_a_eff(a_effs,w0)