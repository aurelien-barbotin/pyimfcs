# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 20:26:38 2021

@author: abarbotin
"""

import pandas as pd
import numpy as np

import multipletau
import tifffile
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

from scipy.special import erf

def make_Gim(a,sigma):
    """Creates a fit function taking into account paramters of PSF
    
    Parameters:
        a (float): pixel side length
        sigma (float): Gaussian standard deviation """
    def G_im(tau,N,D,Ginf):
        """Function to fit ACFs.
        Parameters:
            tau (ndarray): lag time
            D (float): diffusion coeff
            N (float): number of molecules
            Ginf (float): offset, should be around 0"""
        k = a/(2*np.sqrt(D*tau+sigma**2 ) )
        return 1/N*( erf(k)+(np.exp(-k**2)-1)/(k*np.sqrt(np.pi) )) + Ginf
    return G_im

def expf(x,tau,a):
    return a*np.exp(-x/tau)

def bleach_correct_exp(trace, plot = True):
    xx = np.arange(trace.size)
    bounds = ((200,trace[0]-trace[0]/10**-3),(1000*trace.size,trace[0]+trace[0]/10**-3))
    popt, _ = curve_fit(expf,xx,trace,bounds = bounds)
    fit = expf(xx,*popt)
    trace_corrected = trace/np.sqrt(fit/fit[0])+fit[0]*(1-np.sqrt(fit/fit[0]))
    
    if plot:
        plt.figure()
        plt.subplot(121)
        plt.plot(sum_in_blocks(trace, 1))
        plt.plot(sum_in_blocks(fit, 1),color="k", linestyle = "--")
    
        plt.subplot(122)
        plt.plot(trace_corrected)
        
    return trace_corrected

plt.close("all")

path_xls = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_03_31/imFCS20percent_stack.xlsx"

def get_meanD(path):
    fd = pd.read_excel(path, sheet_name = "Fit Parameters (ACF1)")
    
    Dstp = fd[(fd["Parameter"]=="D")].values
    Ds= Dstp[0,1:].astype(float)
    
    Ds_avg = np.mean(Ds[~np.isnan(Ds)])
    print(Ds_avg)

def get_curves_xls(path):
    fd = pd.read_excel(path, sheet_name = "ACF1")
    
    return fd

def get_curves_lagtime(path):
    fd = pd.read_excel(path, sheet_name = "lagtime")
    
    return fd
def sum_in_blocks(a, c):
    # Get extent of each col for summing
    l = c*(len(a)//c)

    # Reshape to 3D considering first l rows, and "cutting" after each c rows
    # Then sum along second axis
    return a[:l].reshape(c,len(a)//c).sum(0)

def correlate_one(img,nsum,i0=0, j0=0):
    u,v,w = img.shape
    trace = img[:,i0:i0+nsum,j0:j0+nsum].sum(axis=(1,2)) 
    corr0 = multipletau.autocorrelate(trace, normalize=True, 
                                     deltat = dt)[1:]    
    return corr0

def mean_acfs(img,nsum):
    u,v,w = img.shape
    corr = 0
    for j in range(v//nsum):
        for k in range(w//nsum):
            trace = img[:,j:j+nsum,k:k+nsum].sum(axis=(1,2)) 
            corr0 = multipletau.autocorrelate(trace, normalize=True, 
                                             deltat = dt)[1:]
            corr+=corr0
    corr/=(v//nsum*w//nsum)
    return corr

load_again = False
if load_again:
    fd = get_curves_xls(path_xls)
    df_lt = get_curves_lagtime(path_xls)

psize = 0.064 #um
fwhm = 0.2 #um
dt = 4.2*10**-3 #s

path = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_03_31/100percent_exp2stack.tif"
name1 = "100% sucrose"
stack = tifffile.imread(path)

nsum = 16
i0 = 15
j0 = 8

Gfit = make_Gim(psize*nsum, fwhm/np.sqrt(8*np.log(2)))
trace = stack[:,i0:i0+nsum,j0:j0+nsum].mean(axis=(1,2))

corr = multipletau.autocorrelate(trace, normalize=True, deltat = dt)[1:]
corr = mean_acfs(stack, nsum)

bounds_fit = ((1000,10**-2,-10**-6),(10**5,10**-1,10**-6))
path2 = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_03_31/imFCS_000_t1_exp2_50percent.tif"
name2 = "50% sucrose"
stack2 = tifffile.imread(path2)


trace2 = stack2[:,i0:i0+nsum,j0:j0+nsum].mean(axis=(1,2))

corr2 = multipletau.autocorrelate(trace2, normalize=True, deltat = dt)[1:]
corr2 = mean_acfs(stack2, nsum)

#corr3 = fd["({}, {})".format(i0,j0)].values[1:]
path3 = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_03_31/imFCS20percent_stack.tif"
name3 = "20% sucrose"
stack3 = tifffile.imread(path3)
trace3 = stack3[:,i0:i0+nsum,j0:j0+nsum].mean(axis=(1,2))
corr3 = mean_acfs(stack3, nsum)

"""corr3 = 0
for k in fd.keys():
    corr3+=fd[k].values[1:]
corr3/=len(fd.keys())
lagtimes = df_lt["Lagtime"].values[1:]
name3 = "20% sucrose"
"""
#----- Fitting ------------------
popt2, _ = curve_fit(Gfit, corr2[:,0], corr2[:,1])
corr2fit = Gfit(corr2[:,0],*popt2)

popt1, _ = curve_fit(Gfit, corr[:,0], corr[:,1])
corrfit = Gfit(corr[:,0],*popt1)

popt3, _ = curve_fit(Gfit, corr3[:,0], corr3[:,1])
print(name1,popt1)
print(name2,popt2)
print(name3,popt3)
corr3fit = Gfit(corr3[:,0],*popt3)

#------ plot ---------------
normalize = True
plt.figure()
plt.subplot(121)
if normalize:
    plt.semilogx(corr[:,0],corr[:,1]/corr[:8,1].mean(), label = name1)
    plt.semilogx(corr2[:,0],corr2[:,1]/corr2[:8,1].mean(), label = name2)
    plt.semilogx(corr2[:,0],corr2fit/corr2[:8,1].mean(), color="k", linestyle="--")
    plt.semilogx(corr2[:,0],corrfit/corr[:8,1].mean(), color="gray", linestyle="--")
    plt.semilogx(corr3[:,0], corr3[:,1]/corr3[:8,1].mean(), label = name3)
    plt.semilogx(corr3[:,0], corr3fit/corr3[:8,1].mean(), color="blue",linestyle="--")
else:
    plt.semilogx(corr[:,0],corr[:,1],label = name1)
    plt.semilogx(corr2[:,0],corr2[:,1], label = name2)
    plt.semilogx(corr3[:,0], corr3[:,1], label = name3)

plt.legend()
plt.xlabel(r"$\rm \tau$")
plt.ylabel(r"$\rm G(\tau)$")
#plt.ylim(top=vm)

plt.subplot(122)
plt.plot(sum_in_blocks(trace, 50),label=name1)
plt.plot(sum_in_blocks(trace2, 50),label=name2)
plt.plot(sum_in_blocks(trace3, 50),label=name3)
plt.xlabel("Time (A.U)")
plt.ylabel("Counts")
plt.legend()

nums = [2, 4, 8, 16]
ms = [correlate_one(stack3, w) for w in nums]
fig, axes  =plt.subplots(1,2)

for j in range(len(nums)):
    corr = ms[j]
    axes[0].semilogx(corr[:,0],corr[:,1]/corr[:8,1].mean(), label = nums[j])
    axes[1].semilogx(corr[:,0],corr[:,1], label = nums[j])
axes[1].legend()
axes[0].set_title("Normalised")
axes[1].set_title("Raw")