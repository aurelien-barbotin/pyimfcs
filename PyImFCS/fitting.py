#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 09:27:18 2021

@author: aurelien
"""

import numpy as np
from scipy.special import erf

def make_Gim2D(a,sigma, ginf=False):
    """Creates a fit function taking into account paramters of PSF
    
    Parameters:
        a (float): pixel side length
        sigma (float): Gaussian standard deviation """
    if ginf:
        def G_im(tau,N,D,Ginf):
            """Function to fit ACFs.
            Parameters:
                tau (ndarray): lag 
                D (float): diffusion coeff
                N (float): number of molecules
                Ginf (float): offset, should be around 0"""
            k = a/(2*np.sqrt(D*tau+sigma**2 ) )
            return 1/N*( erf(k)+(np.exp(-k**2)-1)/(k*np.sqrt(np.pi) ))**2 + Ginf
    else:
        def G_im(tau,N,D):
            """Function to fit ACFs.
            Parameters:
                tau (ndarray): lag 
                D (float): diffusion coeff
                N (float): number of molecules
                Ginf (float): offset, should be around 0"""
            k = a/(2*np.sqrt(D*tau+sigma**2 ) )
            return 1/N*( erf(k)+(np.exp(-k**2)-1)/(k*np.sqrt(np.pi) ))**2
            
    return G_im

def make_Gim3D(a,sigmaxy,sigmaz, ginf=True):
    """Creates a fit function taking into account paramters of PSF
    
    Parameters:
        a (float): pixel side length
        sigmaxy (float): Gaussian standard deviation, lateral
        sigmaz (float): Gaussian standard deviation, axial  """
    if ginf:
        def G_im(tau,N,D,Ginf):
            """Function to fit ACFs.
            Parameters:
                tau (ndarray): lag 
                D (float): diffusion coeff
                N (float): number of molecules
                Ginf (float): offset, should be around 0"""
            k = a/(2*np.sqrt(D*tau+sigmaxy**2 ) )
            gxy = 1/N*( erf(k)+(np.exp(-k**2)-1)/(k*np.sqrt(np.pi) ))**2
            gz = np.sqrt(1+D*tau/sigmaz**2)
            return gz*gxy + Ginf
    else:
        def G_im(tau,N,D):
            """Function to fit ACFs.
            Parameters:
                tau (ndarray): lag 
                D (float): diffusion coeff
                N (float): number of molecules"""
            k = a/(2*np.sqrt(D*tau+sigmaxy**2 ) )
            gxy = 1/N*( erf(k)+(np.exp(-k**2)-1)/(k*np.sqrt(np.pi) ))**2
            gz = np.sqrt(1+D*tau/sigmaz**2)
            return gz*gxy
    return G_im

def fwhm2sigma(fw):
    """let us remember it once and for all"""
    return fw/np.sqrt(8*np.log(2))

def sigma2fwhm(sig):
    """let us remember it once and for all"""
    return sig*np.sqrt(8*np.log(2))