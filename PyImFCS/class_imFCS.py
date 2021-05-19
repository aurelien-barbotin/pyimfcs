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
import h5py
import os

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
    
    def cortrace(tr,popt):
        xt = np.arange(tr.size)
        fi = expf(xt,*popt)
        f0 = popt[0]
        new_tr = tr/np.sqrt(fi/f0)+f0*(1-np.sqrt(fi/f0))
        return (new_tr+popt[2])
    new_tr = cortrace(subtr,popt)*trace.max()
    
    if plot:
        plt.figure()
        plt.plot(subtr)
        plt.plot(trf,color="k",linestyle="--")
        plt.plot(new_tr)
    return new_tr

def blexp_double_offset(trace, plot=False):
    def expf(x,f0,tau,b,tau2,c):
        return f0*(np.exp(-x/tau)+b*np.exp(-x/tau2))+c
    
    subtr = trace/trace.max()
    xt = np.arange(subtr.size)
    p0 = (1,trace.size//10,1,trace.size//10,subtr.min())
    bounds = ((0, 1, 0, 1, 0),
              (2, trace.size*10, 1, trace.size*10,1))
    popt,_ = curve_fit(expf,xt,subtr,p0=p0,bounds = bounds)
    #popt=[4.7*10**6,18]
    trf = expf(xt,*popt)
    
    def cortrace(tr,popt):
        xt = np.arange(tr.size)
        fi = expf(xt,*popt)
        f0 = popt[0]
        new_tr = tr/np.sqrt(fi/f0)+f0*(1-np.sqrt(fi/f0))
        return new_tr/new_tr.mean()
    new_tr = cortrace(subtr,popt)
    
    if plot:
        plt.figure()
        plt.plot(subtr)
        plt.plot(trf,color="k",linestyle="--")
        plt.plot(new_tr)
        
    return new_tr*trace.max()

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
    dic_names = ["correlations", "traces", "parameters_fits","yhat"]
    parameters_names = ["dt","xscale","yscale","path"]
    
    def __init__(self, path, mfactor = 8, background_correction = True, 
                 blcorrf = None,first_n=0, last_n = 0, fitter = None, dt = None,
                 remove_zeroes=False):
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
                self.xscale = 1
                self.yscale = 1
        self.dt = dt
    
        if remove_zeroes:
            print("Achtung! Removing zeroes in the intensity timetrace may lead to artefacts!")
            trace=self.stack.sum(axis=(1,2))
            if np.any(trace==0):
                print('frames removed')
                self.stack=self.stack[trace!=0,:,:]
                
    def save(self,name = None):
        
        if name is None:
            name = os.path.splitext(self.path)[0]+".h5"
        if not name.endswith(".h5"):
            name+=".h5"
        if os.path.isfile(name):
            os.remove(name)
            print("Removing existing file with same name")
        h5f = h5py.File(name, "w")
        
        dicts_to_save = [self.correl_dicts,self.traces_dict,self.parfit_dict, self.yh_dict]
        for j, dic in enumerate(dicts_to_save):
            dname = self.dic_names[j]
            for key, item in dic.items():
                h5f[dname+"/"+str(key)] = item
        
        for pn in self.parameters_names:
            h5f["parameters/"+pn] = getattr(self,pn)
        h5f["blcorrf"] = self.blcorrf.__name__
        
    def load(self,name = None):
        if name is None:
            name = os.path.splitext(self.path)[0]+".h5"
            
        h5f = h5py.File(name, "r")
        dicts_to_load  = {"correlations":self.correl_dicts,
                          "traces":self.traces_dict,
                          "parameters_fits":self.parfit_dict, 
                          "yhat": self.yh_dict}
        
        for j in range(len(dicts_to_load)):
            dname = self.dic_names[j]
            if dname not in h5f.keys():
                print("grougrou")
                continue
            out_dic = dicts_to_load[dname]
            ds = h5f[dname]
            for key in ds.keys():
                dd = ds[key][()]
                out_dic[int(key)] = dd
    
        for par in h5f["parameters"].keys():
            setattr(self,par,h5f["parameters"][par][()])
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
            axes[1].axvline(self.dt*self.stack.shape[0]/1000,color="k",
                            label="Max unbiased transit time")
            axes[1].axvline(self.dt*10,color="k",
                            label="Min unbiased transit time")
            axes[1].legend()
            axes[0].set_title("Raw curves")
            axes[1].set_title("Normalised curves")
            axes[0].set_xlabel(r"$\rm \tau (A.U)$")
            axes[0].set_ylabel(r"$\rm G(\tau)$")
            axes[1].set_xlabel(r"$\rm \tau (A.U)$")
            axes[1].set_ylabel(r"$\rm G(\tau)$")
        return all_corrs
    
    def fit_curves(self,fitter,xmax=None):
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
                    if xmax is None:
                        popt, yh = fitter.fit(corr)
                    else:
                        popt, yh = fitter.fit(corr[corr[:,0]<xmax,:])
                        
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
            axes[0].semilogx(corr[:yh.size,0],yh/a+dz*jj, color="k",linestyle="--")
            axes[1].semilogx(corr[:yh.size,0],yh/a-corr[:yh.size,1]/a+dz*jj, 
                             label = "curve ({},{})".format(j,k))
            jj+=1
        axes[1].legend()
        
    def plot_D(self, show_acceptable=True):
        nsums = sorted(self.parfit_dict.keys())
        nsums = np.asarray(nsums)
        ds_means = list()
        ds_std = list()
        
        if show_acceptable:
            psize = self.fitter.parameters_dict["a"]
            if "sigmaxy" in self.fitter.parameters_dict:
                sigmaxy = self.fitter.parameters_dict["sigmaxy"]
            else:
                sigmaxy = self.fitter.parameters_dict["sigma"]
            sigmaxy = sigmaxy*np.sqrt(8*np.log(2)) #fwhm
            observation_sizes = np.sqrt((nsums*psize)**2+sigmaxy**2)
            if self.fitter.name=="2D":
                factor = 4
            elif self.fitter.name=="3D":
                factor=6
            dmax = observation_sizes**2/(factor*self.dt*10)
            dmins = observation_sizes**2/(factor*self.dt*self.stack.shape[0]/1000)
        for ns in nsums:
            ds = self.parfit_dict[ns][:,:,1]
            ds_means.append(np.mean(ds))
            ds_std.append(np.std(ds))
        plt.figure()
        if show_acceptable:
            plt.plot(nsums,dmins,color="gray")
            plt.plot(nsums,dmax,color="gray")
        plt.errorbar(nsums, ds_means, yerr=ds_std,capsize=5)
        plt.xlabel("Binning size")
        plt.ylabel("D (um2/s)")
        return nsums, ds_means, ds_std
    
    def plot_taus(self, show_acceptable=True):
        nsums = sorted(self.parfit_dict.keys())
        nsums = np.asarray(nsums)
        ds_means = list()
        ds_std = list()
        
        if show_acceptable:
            psize = self.fitter.parameters_dict["a"]
            if "sigmaxy" in self.fitter.parameters_dict:
                sigmaxy = self.fitter.parameters_dict["sigmaxy"]
            else:
                sigmaxy = self.fitter.parameters_dict["sigma"]
            sigmaxy = sigmaxy*np.sqrt(8*np.log(2)) #fwhm
            observation_sizes = np.sqrt((nsums*psize)**2+sigmaxy**2)

            dmax = np.ones_like(nsums)*(self.dt*10)
            dmins = np.ones_like(nsums)*self.dt*self.stack.shape[0]/1000
        
        for j,ns in enumerate(nsums):
            ds = self.parfit_dict[ns][:,:,1]
            taus = observation_sizes[j]**2/ds
            ds_means.append(np.mean(taus))
            ds_std.append(np.std(taus))
        ds_means=np.asarray(ds_means)
        
        plt.figure()
        if show_acceptable:
            plt.plot(nsums,dmins,color="gray")
            plt.plot(nsums,dmax,color="gray")
        plt.errorbar(nsums, ds_means, yerr=ds_std,capsize=5)
        plt.xlabel("Binning size")
        plt.ylabel("tau (s)")
        return nsums, ds_means, ds_std
    
    def parameter_map(self,nsum = None,parn=1):
        if nsum is None:
            nsum = min(self.parfit_dict.keys())
        out = self.parfit_dict[nsum][:,:,parn]
        out = np.repeat(out,nsum,axis=0)
        out = np.repeat(out,nsum,axis=1)
        return out
    
    def plot_parameter_maps(self,nsums, parn=1):
        assert len(nsums)>=1
        assert len(nsums)<=5
        
        nr = int(np.sqrt(len(nsums)+1))
        nc = (len(nsums)+1)//nr
        if nr*nc<len(nsums)+1:
            nc+=1
        fig,axes = plt.subplots(nr,nc, sharex = True, sharey = True)
        axes=axes.ravel()
        ax0 = axes[0]
        ax0.imshow(self.stack.mean(axis=0),cmap="gray")
        parmaps = list()
        for j in range(len(nsums)):
            nsum = nsums[j]
            parmap = self.parameter_map(nsum=nsum,parn=parn)
            parmaps.append(parmap)
        vmin = min([w.min() for w in parmaps])
        vmax = max([w.max() for w in parmaps])
        for j in range(len(nsums)):
            im=axes[j+1].imshow(parmaps[j],cmap="hot")
            fig.colorbar(im,ax=axes[j+1])
            
    def plot_intensity_correlation(self,nsum,parn=1):
        parmap = self.parfit_dict[nsum][:,:,parn]
        parameters = list()
        intensities = list()
        stack_avg = self.stack.mean(axis=0)
        for j in range(parmap.shape[0]):
            for k in range(parmap.shape[0]):
                parameters.append(parmap[j,k])
                intensities.append(stack_avg[j*nsum:(j+1)*nsum, 
                                             k*nsum:(k+1)*nsum].mean())
        plt.figure()
        plt.scatter(parameters,intensities)
        plt.xlabel("Parameter")
        plt.ylabel("Intensity")
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
