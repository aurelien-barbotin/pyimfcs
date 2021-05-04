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

def bleaching_correct_sliding(trace, plot = False, wsize = 5000):
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
    return corrf*trace

def get_image_metadata(path):    
    img = tifffile.TiffFile(path)
    meta_dict = img.imagej_metadata
    description = meta_dict.pop('Info')
    description = description.split('\n')
    for d in description:
        if len(d)>1 and '=' in d:
            k, val = d.split('=')
            k = k.strip(' ')
            val = val.strip(' ')
            try:
                meta_dict[k] = float(val)
            except:
                
                meta_dict[k] = val
                
    return meta_dict
    
class StackFCS(object):
    def __init__(self, path, mfactor = 8, background_correction = True, 
                 blcorrf = None,first_n=0, last_n = 0, fitter = None, dt = None):
        self.path = path
        self.stack = tifffile.imread(path)
        self.stack = self.stack[first_n:self.stack.shape[0]-last_n]
        self.fitter = fitter
        
        if background_correction:
            self.stack = self.stack - self.stack.min()

        self.blcorrf = blcorrf
           
        self.correl_dicts = {}
        self.traces_dict = {}
        self.parfit_dict = {}
        self.yh_dict = {}
        if dt is None:
            try:
                metadata = get_image_metadata(path)
                dt = metadata['finterval']
                xscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1']*10**6
                yscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1']*10**6
                self.xscale = xscale
                self.yscale = yscale
                if xscale!=yscale:
                    raise ValueError('Not square pixels')
            except:
                print('error loading metadata')
                dt = 1
                xscale = 1
                yscale = 1
        self.dt = dt
    
    def correlate_stack(self,nSum):
        """Only method that correlates """
        if nSum>self.stack.shape[1] or nSum>self.stack.shape[2]:
            raise ValueError("Trying to sum more pixels than present in the stack")
        if nSum not in self.correl_dicts.keys():
            print("correlating stack, binning {}".format(nSum))
            correls = []
            traces = []
            u,v,w = self.stack.shape
            for i in range(v//nSum):
                ctmp = []
                trtmp = []
                for j in range(w//nSum):
                    trace = self.stack[:,i*nSum:i*nSum+nSum,j*nSum:j*nSum+nSum].mean(axis=(1,2))
                    if self.blcorrf is not None:
                        trace = self.blcorrf(trace)
                    corr = multipletau.autocorrelate(trace, normalize=True, deltat = self.dt)[1:]
                    ctmp.append(corr)
                    trtmp.append(trace)
                correls.append(ctmp)
                traces.append(trtmp)
            correls = np.asarray(correls)
            
            self.correl_dicts[nSum] = correls
            self.traces_dict[nSum] = np.asarray(traces)
            
    def get_curve(self, i0 = 0, j0 =0, nSum=1):
        self.correlate_stack(nSum)
            
        correls = self.correl_dicts[nSum]
        correl = correls[i0,j0]
        return correl
    
    def plot_curve(self, i0 = 0, j0 =0, nSum=1):
        correl = self.get_curve(i0 = i0, j0=j0,nSum=nSum)
        
        plt.figure()
        plt.semilogx(correl[:,0], correl[:,1])
    
    def get_all_curves(self,nSum=1, spacing=0, npts = None, plot = True):
        self.correlate_stack(nSum)
        correls = self.correl_dicts[nSum]
        u,v = correls.shape[0], correls.shape[1]
        spc = np.percentile(correls[:,:,:,1],95)*spacing
        if plot:
            plt.figure()
            for j in range(u):
                for k in range(v):
                    corr = correls[j,k]
                    plt.semilogx(corr[:,0], corr[:,1]+(j*u+k)*spc)
        return correls
    
    def get_correlation_dict(self):
        return self.correl_dicts
    
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
    
    def fit_curves(self,fitter):
        self.fitter = fitter
        nsums = self.correl_dicts.keys()
        for nsum in nsums:
            self.fitter.set_sum(nsum)
            correls = self.correl_dicts[nsum]
            popts = []
            yhs = []
            for j in range(correls.shape[0]):
                popt_tmp=[]
                yh_tmp = []
                for k in range(correls.shape[1]):
                    corr = correls[j,k]
                    popt, yh = fitter.fit(corr)
                    
                    popt_tmp.append(popt)
                    yh_tmp.append(yh)
                popts.append(popt_tmp)
                yhs.append(yh_tmp)
            self.parfit_dict[nsum] = np.array(popts)
            self.yh_dict[nsum] = np.array(yhs)
        
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
    
    def plot_random_intensity(self, nSum = None):
        if nSum is None:
            nSum = min(self.traces_dict.keys())
        traces_arr = self.traces_dict[nSum]
        trace_raw = self.stack.mean(axis=(1,2))
        u,v = traces_arr.shape[:2]
        u1 = np.random.choice(u)
        v1 = np.random.choice(v)
        trace = traces_arr[u1,v1]
        
        plt.figure()
        plt.plot(trace, label = "Corrected")
        plt.axhline(trace.mean(),color="k",linestyle="--")
        plt.plot(trace_raw,label = "Raw average intensity")
        plt.xlabel("Time (frames)")
        plt.ylabel("Intensity")
        plt.legend()
    
    def plot_fits(self,nSum,maxcurves=None,dz=0.2):
        curves = self.correl_dicts[nSum]
        sp1 = np.asarray(curves.shape)
        fits = self.yh_dict[nSum]
        fig,axes = plt.subplots(1,2,sharex=True,sharey=True)
        jj = 0
        
        indices = np.arange(sp1[0]*sp1[1])
        indices1 = indices//sp1[1]
        indices2 = indices-(indices1)*sp1[1]
        if maxcurves is not None:
            indices = np.random.choice(np.arange(sp1[0]*sp1[1]),maxcurves)
            
        for i in indices:
            j = indices1[i]
            k = indices2[i]
            corr = curves[j,k]
            yh = fits[j,k]
            a = corr[:3,1].mean()
            axes[0].semilogx(corr[:,0],corr[:,1]/a+dz*jj)
            axes[0].semilogx(corr[:,0],yh/a+dz*jj, color="k",linestyle="--")
            axes[1].semilogx(corr[:,0],yh/a-corr[:,1]/a+dz*jj, 
                             label = "curve ({},{})".format(j,k))
            jj+=1
        axes[1].legend()
        
    def plot_D(self):
        nsums = sorted(self.parfit_dict.keys())
        ds_means = list()
        ds_std = list()
        
        for ns in nsums:
            ds = self.parfit_dict[ns][:,:,1]
            ds_means.append(np.mean(ds))
            ds_std.append(np.std(ds))
        plt.figure()
        plt.errorbar(nsums, ds_means, yerr=ds_std,capsize=5)
        plt.xlabel("Binning size")
        plt.ylabel("D (um2/s)")
        return nsums, ds_means, ds_std
    
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
