# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 14:57:50 2021

@author: abarbotin
"""
import tifffile
import multipletau
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
    
from PIL import Image
from PIL.TiffTags import TAGS

def expf(x,tau,a):
    return a*np.exp(-x/tau)

def bleaching_correct(stack, plot = True):
    trace = stack.sum(axis=(1,2))
    x = np.arange(trace.size)
    popt, _ = curve_fit(expf, x, trace)
    
    if plot:
        yh = expf(x,*popt)
        plt.figure()
        plt.plot(x,trace)
        plt.plot(x,yh)

def bleaching_correct_sliding(stack, plot = True, wsize = 5000):
    trace = stack.sum(axis=(1,2))
    u = trace.size//wsize
    new_trace = trace.copy()
    m1 = trace[:wsize].mean()
    corrf = np.ones_like(trace).astype(float)
    for j in range(1,u+1):
        cf = m1/trace[j*wsize:j*wsize+wsize].mean()
        new_trace[j*wsize:j*wsize+wsize] = trace[j*wsize:j*wsize+wsize]*cf
        corrf[j*wsize:j*wsize+wsize] = cf
    
    if plot:
        plt.figure()
        plt.plot(trace)
        plt.plot(new_trace)
    return corrf

def get_image_metadata(path):    
    with Image.open(path) as img:
        meta_dict = {TAGS[key] : img.tag[key] for key in img.tag.keys()}
        
    description = meta_dict['ImageDescription']
    description = description[0].split('\n')
    description_dict = {}
    for d in description:
        if len(d)>0:
            k, val = d.split('=')
            try:
                description_dict[k] = float(val)
            except:
                description_dict[k] = val
                
    return description_dict
    
class StackFCS(object):
    def __init__(self, path, mfactor = 8, dt=1, background_correction = True, 
                 blcorrf = None,first_n=0):
        self.path = path
        self.stack = tifffile.imread(path)[first_n:]
        if background_correction:
            self.stack = self.stack - self.stack.min()
            print("background correction on")
        if blcorrf is not None:
            corrfactor = blcorrf(self.stack)
            self.stack = self.stack*corrfactor[:,np.newaxis,np.newaxis]
            print("bleaching correction on")
        self.correl_dicts = {}
        self.dt = dt #TODO
    
    def correlate_stack(self,nSum):
        
        if nSum>self.stack.shape[1] or nSum>self.stack.shape[2]:
            raise ValueError("Trying to sum more pixels than present in the stack")
        if nSum not in self.correl_dicts.keys():
            print("correlating stack, binning {}".format(nSum))
            correls = []
            u,v,w = self.stack.shape
            for i in range(v//nSum):
                ctmp = []
                for j in range(w//nSum):
                    trace = self.stack[:,i*nSum:i*nSum+nSum,j*nSum:j*nSum+nSum].mean(axis=(1,2))
                    corr = multipletau.autocorrelate(trace, normalize=True, deltat = self.dt)[1:]
                    ctmp.append(corr)
                correls.append(ctmp)
            correls = np.asarray(correls)
            self.correl_dicts[nSum] = correls
            
    def get_curve(self, i0 = 0, j0 =0, nSum=1):
        self.correlate_stack(nSum)
            
        correls = self.correl_dicts[nSum]
        correl = correls[i0,j0]
        return correl
    
    def plot_curve(self, i0 = 0, j0 =0, nSum=1):
        correl = self.get_curve(i0 = i0, j0=j0,nSum=nSum)
        
        plt.figure()
        plt.semilogx(correl[:,0], correl[:,1])
    
    def get_all_curves(self,nSum=1, spacing=0, npts = None):
        self.correlate_stack(nSum)
        correls = self.correl_dicts[nSum]
        u,v = correls.shape[0], correls.shape[1]
        spc = np.percentile(correls[:,:,:,1],95)*spacing
        
        plt.figure()
        for j in range(u):
            for k in range(v):
                corr = correls[j,k]
                plt.semilogx(corr[:,0], corr[:,1]+(j*u+k)*spc)
                
    def average_curve(self, nSum=1, plot = False):
        self.correlate_stack(nSum)
        
        correls = self.correl_dicts[nSum]
        print("Nsum avg curve",nSum)
        avg = correls[:,:,:,1].mean(axis=(0,1))
        if plot:
            plt.figure()
            plt.semilogx(correls[0,0,:,0],avg)
            plt.axhline(0,linestyle="--",color="k")
            plt.title("Average of binning {}".format(nSum))
        return np.array([correls[0,0,:,0],avg]).T
    
    def trace(self, plot=True):
        tr = self.stack.sum(axis=(1,2))
        if plot:
            xtr = np.arange(tr.size)*self.dt
            plt.figure()
            plt.plot(xtr,tr)
            plt.xlabel("time")
            plt.ylabel("Counts")
        return tr
    
    def binned_average_curves(self, sum_list, plot=True, n_norm = 8):
        if plot:
            fig, axes = plt.subplots(1,2)
        all_corrs = []
        for sl in sum_list:
            corr = self.average_curve(nSum=sl)
            all_corrs.append(corr)
            if plot:
                axes[0].semilogx(corr[:,0], corr[:,1], 
                                 label = "Binning {}".format(sl))
                axes[1].semilogx(corr[:,0], corr[:,1]/corr[:n_norm,1].mean(),
                                 label = "Binning {}".format(sl))
        if plot:
            axes[1].legend()
            axes[0].set_title("Raw curves")
            axes[1].set_title("Normalised curves")
            axes[0].set_xlabel(r"$\rm \tau (A.U)$")
            axes[0].set_ylabel(r"$\rm G(\tau)$")
            axes[1].set_xlabel(r"$\rm \tau (A.U)$")
            axes[1].set_ylabel(r"$\rm G(\tau)$")
        return all_corrs
    
    def plot_amplitudes(self,sum_list):
        averages = self.binned_average_curves(sum_list, plot=False)
        nvals = []
        for j in range(len(sum_list)):
            nn = averages[j][:8,1].mean()
            nvals.append(nn)
        plt.figure()
        plt.plot(sum_list, np.sqrt(nvals), label="measured")
        plt.plot(sum_list, np.sqrt(nvals)[0]*sum_list[0]/np.array(sum_list), 
                 color="k", linestyle="--",label="theory")
        plt.xlabel("Binning")
        plt.ylabel("Sqrt Curve amplitude")
        plt.xscale('log')
        plt.legend()
        
if __name__=="__main__":
    plt.close('all')
    path = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_04_07/40percent_imFCS_002_t1_stack.tif"
    # path = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_04_07/40percent_imFCS_001_t1_stack.tif"
    # path = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_04_07/test3_40percent/imFCS_000_stack.tif"
    path = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_03_31_Zeiss/50percentSucrose/ImFCS2_bleachingbefore.tif"
    path = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_03_31_Zeiss/100percentSucrose/imFCS3_40min.tif"
    path = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_03_24_Elyra7/ImFCS6.tif"
    
    path = "C:/Users/abarbotin/Documents/Python Scripts/FCS_analysis/simulation4.tif"
    path = "C:/Users/abarbotin/Desktop/analysis_tmp/2021_04_09/2_green_beads/Image 13.tif"
    path = "C:/Users/abarbotin/Desktop/analysis_tmp/2021_04_09/2_green_beads/Image 21_withAF.tif"
    path = "C:/Users/abarbotin/Desktop/analysis_tmp/2021_04_09/2_green_beads/Image 20_higherpower.tif"
    path = "C:/Users/abarbotin/Desktop/analysis_tmp/2021_04_09/postprocessed/Image 21_withAF_zoom1.tif"
    # path = "C:/Users/abarbotin/Desktop/analysis_tmp/bilayer_wohland.tif"
    #path = "C:/Users/abarbotin/Desktop/analysis_tmp/2020_04_07/60percent_imFCS_001_t1_stack.tif"
    stack = StackFCS(path, first_n=0, blcorrf=None)
    stack.plot_curve(nSum=4)
    stack.get_all_curves(nSum=8, spacing = 0.2)
    t1 = stack.average_curve(nSum=4)
    t2 = stack.average_curve(nSum=8)
    
    stack.trace()
    
    plt.figure()
    plt.subplot(121)
    plt.semilogx(t1[:,0],t1[:,1],label="4 pixel binning")
    plt.semilogx(t2[:,0],t2[:,1], label="8 pixel binning")
    plt.title("Average curves")
    plt.legend()
    plt.subplot(122)
    plt.semilogx(t1[:,0],t1[:,1]/t1[:3,1].mean(),label="4 pixel binning")
    plt.semilogx(t2[:,0],t2[:,1]/t2[:3,1].mean(), label="8 pixel binning")
    plt.title("Normalised")

    nsums = [2,4,8,12]
    corrs = stack.binned_average_curves(nsums, n_norm=3)
    stack.plot_amplitudes(nsums)
