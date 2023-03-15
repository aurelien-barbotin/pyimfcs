#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 09:27:18 2021

@author: aurelien
"""

import os
import sys
import json

import numpy as np

from inspect import signature
from scipy.special import erf
from scipy.optimize import curve_fit

BUNDLE_DIR = getattr(sys, '_MEIPASS', os.path.abspath(os.path.dirname(__file__)))

def fwhm2sigma(fw):
    """let us remember it once and for all"""
    return fw/np.sqrt(8*np.log(2))

def sigma2fwhm(sig):
    """let us remember it once and for all"""
    return sig*np.sqrt(8*np.log(2))

# --- fit functions --------------------

def gim2D(a=0.1,sigma=0.1, ginf=False,**kwargs):
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

def gim2D_2components(a=0.1,sigma=0.1, ginf=False,**kwargs):
    """Creates a fit function taking into account paramters of PSF
    
    Parameters:
        a (float): pixel side length
        sigma (float): Gaussian standard deviation """
    if ginf:
        def G_im(tau,N,D1,D2,A2,Ginf):
            """"""
            k1 = a/(2*np.sqrt(D1*tau+sigma**2 ))
            k2 = a/(2*np.sqrt(D2*tau+sigma**2 ))
            return 1/N*( 
                erf(k1)+(np.exp(-k1**2)-1)/(k1*np.sqrt(np.pi)) +
                A2* (erf(k2)+(np.exp(-k2**2)-1)/(k2*np.sqrt(np.pi)))
                        )**2 + Ginf
    else:
        def G_im(tau,N,D1,D2,A2):
            """"""
            k1 = a/(2*np.sqrt(D1*tau+sigma**2 ) )
            k2 = a/(2*np.sqrt(D2*tau+sigma**2 ) )
            return 1/N*( 
                erf(k1)+(np.exp(-k1**2)-1)/(k1*np.sqrt(np.pi)) +
                A2* (erf(k2)+(np.exp(-k2**2)-1)/(k2*np.sqrt(np.pi)))
                        )**2
            
    return G_im

def gim3D(a=0.1,sigmaxy=0.1,sigmaz=0.5, ginf=False,**kwargs):
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

fit_functions = {"2D":gim2D,
                 "3D":gim3D,
                 "2D_2c":gim2D_2components}
# microscope-dependent parameters for fitting model
fit_parameters_dict = {"2D":["sigma","ginf"],
                       "3D":["sigma","sigmaz","ginf"],
                       "2D_2c":['sigma',"ginf"]
                       }

fit_parameter_types = {'sigma': float, 'sigmaz':float,"ginf":bool}

fit_p0 = {"2D": [lambda x: 1/x[0,1]/3, lambda x: 0.23**2/4/np.median(x[:,0])],
          "2D_2c": [lambda x: 1/x[0,1]/3, lambda x: 0.23/4/np.median(x[:,0]),
                    lambda x: 0.23/2/np.median(x[:,0]),lambda x:0.5],
          "3D": [lambda x: 1/x[0,1]/3, lambda x: 0.23/4/np.median(x[:,0])]
          }

class Fitter(object):
    
    def __init__(self,parameters_dict, p0 = None,bounds=None):
        """Object used to fit curves. Creates the fitting model from the PSF parameters
        and given a desired fitting type.
        Parameters:
            parameters_dict (dict): dictionary input to one of the models
            p0 (tuple): initial guess for fitting, for scipy.optimize.curve_fit
            bounds (tuple): bounds for method curve_fit"""
        ginf = parameters_dict['ginf']
        mtype = parameters_dict['mtype']
        self.parameters_dict = parameters_dict
        self.ginf = ginf
        self.mtype = mtype
        if self.mtype not in fit_functions:
            raise KeyError('Unknown fitter')
        
        self.full_parameters_dict = parameters_dict.copy() # takes nsum in account
        self.fitter = fit_functions[self.mtype](**parameters_dict)
        self.p0 = p0
        self.p0f = fit_p0[self.mtype]
        self.bounds = bounds
        
    def set_sum(self,nsum):
        self.full_parameters_dict["a"] = self.parameters_dict["a"]*nsum
        self.fitter = fit_functions[self.mtype](**self.full_parameters_dict)
        
    def fit(self,curve):
        try:
            p0 = [f(curve) for f in self.p0f]
            if self.ginf:
                p0.append(0)
            self.p0 = tuple(p0)
            
            if self.bounds is not None:
                popt,_ = curve_fit(self.fitter,curve[:,0],curve[:,1], p0=self.p0,
                                   bounds=self.bounds)
            else:
                popt,_ = curve_fit(self.fitter,curve[:,0],curve[:,1], p0=self.p0)
            yh = self.fitter(curve[:,0],*popt)
            return popt, yh
        except Exception as e:
            # raise e
            print("Fitting error")
            print(e)
            sig = signature(self.fitter)
            popt = [-1]*(len(sig.parameters)-1)
            yh = np.zeros_like(curve[:,0])
            return popt, yh

def create_model_dict(name,parameters_dict):
    """Method used to create a parameter dictionary"""
    filename = BUNDLE_DIR+"/models/"+name+".json"
    with open(filename,"w") as f:
        json.dump(parameters_dict,f,indent=2)

if __name__=='__main__':
    pardict = {"sigma":0.20,
               "ginf":True,
               "mtype":"2D"}
    create_model_dict("2D_lens1X_legacy",pardict)