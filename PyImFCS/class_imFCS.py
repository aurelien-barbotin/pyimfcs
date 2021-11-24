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
    
import h5py
import os

from skimage.filters import threshold_otsu
from scipy.signal import fftconvolve
from scipy.ndimage import label
from PyImFCS.shift_correction import stackreg

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

def blexp_double_offset(trace, plot=False, downsample = True):
    def expf(x,f0,tau,b,tau2,c):
        return f0*(np.exp(-x/tau)+b*np.exp(-x/tau2))+c
    
    subtr = trace/trace.max()
    if downsample:
        nn0 = subtr.size//1000 # take one point every nn0
        xt = np.arange(subtr.size)
        subtr = subtr[::nn0]
        xt = xt[::nn0]
    
    p0 = (1,trace.size//10,1,trace.size//10,subtr.min())
    bounds = ((0, 1, 0, 1, 0),
              (2, trace.size*10, 1, trace.size*10,0.8))
    
    try:
        popt,_ = curve_fit(expf,xt,subtr,p0=p0,bounds = bounds)
    except:
        print("Fitting failed arzo")
        return trace
    #popt=[4.7*10**6,18]
    
    # redefine subtr and xt to have right number of points
    subtr = trace/trace.max()
    xt = np.arange(subtr.size)
    trf = expf(xt,*popt)
    new_tr = cortrace(subtr,trf)
    
    if plot:
        plt.figure()
        plt.plot(subtr)
        plt.plot(trf,color="k",linestyle="--")
        plt.plot(new_tr)
        
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

def get_image_metadata(path):
    img = tifffile.TiffFile(path)
    meta_dict = img.imagej_metadata
    description = meta_dict.pop('Info')
    description = description.split('\n')
    for d in description:
        if len(d)>1 and '=' in d:
            oo = d.split('=')
            if len(oo)==2:
                k, val = oo
            elif len(oo)>2:
                k = oo[0]
                val = "=".join(oo[1:])
            k = k.strip(' ')
            val = val.strip(' ')
            try:
                meta_dict[k] = float(val)
            except:
                meta_dict[k] = val
    return meta_dict

def save_tiff_withmetadata(file, st, metadata):
    file = "/home/aurelien/Data/2021_06_03/test.tif"
    path = "/home/aurelien/Data/2021_06_03/imFCS1.tif"
    img = tifffile.TiffFile(path)
    meta_dict = img.imagej_metadata
    
    writer = tifffile.TiffWriter(file,imagej=True)
    writer.write(st,metadata=meta_dict)
    
def new_chi_square(y,yh):
    diff = (y-yh)/yh[0]
    diffpos = diff>0
    diffneg = diff<0
    
    diff =  diff**2
    
    chi = 0
    
    lab, num_features = label(diffpos)
    for j in range(1,num_features+1):
        msk = lab==j
        chi+=diff[msk].sum()*(np.count_nonzero(msk)-1)
        
    lab, num_features = label(diffneg)
    for j in range(1,num_features+1):
        msk = lab==j
        chi+=diff[msk].sum()*(np.count_nonzero(msk)-1)
    
    return chi/yh.size

class StackFCS(object):
    dic_names = ["correlations", "traces", "parameters_fits","yhat"]
    # parameters to save
    parameters_names = ["dt","xscale","yscale","path", "nreg", "shifts", 
                        "first_n","last_n","clipval"]
    
    def __init__(self, path, mfactor = 8, background_correction = True, 
                 blcorrf = None,first_n=0, last_n = 0, fitter = None, dt = None,
                 remove_zeroes=False, clipval = 0, load_stack = True):

        self.path = path
        self.load_stack = load_stack
        self.first_n = first_n
        self.last_n = last_n
        if load_stack:
            self.stack = tifffile.imread(path)
            self.stack = self.stack[self.first_n:self.stack.shape[0]-self.last_n]
        else:
            self.stack = np.zeros((5,5,5))
            
        self.fitter = fitter
        
        
        self.threshold_map = None
        
        # removes clipval points before and after the intensity timetrace
        # before correlation. To remove artefacts from bleaching correction
        # or image registration
        self.clipval = clipval 
        
        if background_correction:
            self.stack = self.stack - self.stack.min()

        self.blcorrf = blcorrf
        # shift correction
        self.nreg = 0
        self.shifts = np.zeros(1)
        
        # resuts dictionaries
        self.correl_dicts = {}
        self.traces_dict = {}
        self.parfit_dict = {}
        self.yh_dict = {}
        self.chisquares_dict = {}
        
        if dt is None and load_stack:
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
        else:
            self.xscale = 1
            self.yscale = 1
            dt = 1
            
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
        if self.blcorrf is not None:
            h5f["blcorrf"] = self.blcorrf.__name__
        else:
            h5f["blcorrf"] = 'None'
            
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
    
    def registration(self,nreg, plot = False):
        self.stack, shifts = stackreg(self.stack,nreg,plot=plot)
        self.nreg = nreg
        self.shifts = shifts
    
    def set_threshold_map(self,th_map):
        self.threshold_map = th_map
        
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
                    trace = self.stack[:,i*nSum:i*nSum+nSum,
                                       j*nSum:j*nSum+nSum].mean(axis=(1,2))
                    if self.blcorrf is not None:
                        trace = self.blcorrf(trace)
                        
                    if self.clipval>0:
                        trace = trace[self.clipval:-self.clipval]
                        
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
    
    def get_all_curves(self,nSum=1, spacing=0, npts = None, plot = True):
        self.correlate_stack(nSum)
        correls = self.correl_dicts[nSum]
        
        if self.threshold_map is not None:
            pass
        
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
        if self.threshold_map is not None:
            th = self.get_threshold_map(nSum)
            cs=correls[th]
            avg = cs[:,:,1].mean(axis=(0))
        else:
            avg = correls[:,:,:,1].mean(axis=(0,1))
        
        print("Nsum avg curve",nSum)
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
    
    def plot_sum_img(self):
        plt.figure()
        plt.imshow(self.stack.mean(axis=0))
        
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
            axes[1].axvline(self.dt*self.stack.shape[0]/100,color="k",
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
            axes[0].axhline(0,color="k")
            axes[1].axhline(0,color="k")
        return all_corrs
    
    def fit_curves(self,fitter,xmax=None):
        self.fitter = fitter
        nsums = self.correl_dicts.keys()
        for nsum in nsums:
            self.fitter.set_sum(nsum)
            correls = self.correl_dicts[nsum]
            popts = []
            yhs = []
            chisquares = []
            for j in range(correls.shape[0]):
                popt_tmp=[]
                yh_tmp = []
                chisquares_tmp = []
                for k in range(correls.shape[1]):
                    corr = correls[j,k]
                    if xmax is None:
                        popt, yh = fitter.fit(corr)
                    else:
                        corr = corr[corr[:,0]<xmax,:]
                        popt, yh = fitter.fit(corr)
                    # sig = signature(fitter.fitter)
                    # ddof=len(sig.parameters)-1
                    #print(corr[:,1].dtype,yh.dtype, corr[:,1].dtype,yh.dtype)
                    chi = np.square(corr[:,1]-yh).mean()/yh[0]**2
                    popt_tmp.append(popt)
                    yh_tmp.append(yh)
                    chisquares_tmp.append(chi)
                    
                popts.append(popt_tmp)
                yhs.append(yh_tmp)
                chisquares.append(chisquares_tmp)
                
            self.parfit_dict[nsum] = np.array(popts)
            self.yh_dict[nsum] = np.array(yhs)
            self.chisquares_dict[nsum] = np.array(chisquares)
            
    def parameter_map(self,nsum = None,parn=1):
        print('Caution! You are using a resampled parameter map')
        if nsum is None:
            nsum = min(self.parfit_dict.keys())
        out = self.parfit_dict[nsum][:,:,parn]
        out = np.repeat(out,nsum,axis=0)
        out = np.repeat(out,nsum,axis=1)
        return out
    
    def get_param_threshold(self,nsum,thf=None, plot = False, parn = 1):
                
        img = self.stack.mean(axis=0).astype(float)
        if thf is None:
            thf = threshold_otsu
        thresholded = img>thf(img)

        uu = thresholded.shape[0] - thresholded.shape[0]%nsum
        vv = thresholded.shape[1] - thresholded.shape[1]%nsum
        
        to_keep=self.get_threshold_map(nsum, thf=thf)
        if plot:
            plt.figure()
            plt.subplot(121)
            plt.imshow(img)
            
            plt.subplot(122)
            plt.imshow(to_keep)
            
        return self.parfit_dict[nsum][:uu,:vv,parn][to_keep]
    
    def get_param_coord(self, nsum,i0,j0,parn=1, exclude_neg = True):
        """Get the value of given parameters for all binning values below nsum"""
        sums = self.correl_dicts.keys()
        sums = sorted([w for w in sums if w<=nsum])
        ds_means=list()
        ds_std=list()
        for ns in sums:
            i00 = int(np.ceil(i0*nsum/ns))
            i01 = int(np.floor((i0+1)*nsum/ns))
            
            j00 = int(np.ceil(j0*nsum/ns))
            j01 = int(np.floor((j0+1)*nsum/ns))
            # print("ns",ns,i00,i01,"j",j00,j01)
            #print("ns {}, i {}-{}, j {}-{}".format(ns,i00,i01,j00,j01))
            
            ds = self.parfit_dict[ns][i00:i01,j00:j01, parn]
            if exclude_neg:
                ds = ds[ds>=0]
            ds_means.append(np.mean(ds))
            ds_std.append(np.std(ds))
        sums = np.asarray(sums)
        ds_means = np.asarray(ds_means)
        ds_std = np.asarray(ds_std)
        mask = ~np.isnan(ds_means)
        return sums[mask], ds_means[mask], ds_std[mask]
    
    def get_acf_coord(self, nsum,i0,j0,parn=1, average = True):
        
        sums = self.correl_dicts.keys()
        sums = sorted([w for w in sums if w<=nsum])
        all_corrs=list()
        all_yhs =list()
        all_ns = list()
        for ns in sums:
            i00 = int(np.ceil(i0*nsum/ns))
            i01 = int(np.floor((i0+1)*nsum/ns))
            
            j00 = int(np.ceil(j0*nsum/ns))
            j01 = int(np.floor((j0+1)*nsum/ns))
            corrs = self.correl_dicts[ns][i00:i01,j00:j01]
            yhs=self.yh_dict[ns][i00:i01,j00:j01]
            if average:
                corrs = corrs.mean(axis=(0,1))
                yhs=yhs.mean(axis=(0,1))
            if not np.isnan(corrs).all():
                all_corrs.append(corrs)
                all_yhs.append(yhs)
                all_ns.append(ns)
        sums = np.asarray(sums)
        all_corrs = np.asarray(all_corrs)
        all_yhs = np.asarray(all_yhs)
        return all_ns, all_corrs, all_yhs
    
    def get_threshold_map(self,nsum, thf = None):
        img = self.stack.sum(axis=0).astype(float)
        
        print('thresholding function is always otsu')
        if thf is None:
            thf = threshold_otsu
        thresholded = (img>thf(img)).astype(float)
        
        uu = thresholded.shape[0] - thresholded.shape[0]%nsum
        vv = thresholded.shape[1] - thresholded.shape[1]%nsum
        
        out = np.zeros((thresholded.shape[0]//nsum, thresholded.shape[1]//nsum))
        for j in range(nsum):
            for k in range(nsum):
                out+=thresholded[:uu,:vv][j::nsum,k::nsum]
        
        to_keep = out==nsum**2
        return to_keep
    
    def downsample_image(self, nsum):
        u,v,w=self.stack.shape
        
        img = self.stack.mean(axis=0)
        out = np.zeros((v//nsum,w//nsum))
        for i in range(v//nsum):
            for j in range(w//nsum):
                px = img[i*nsum:i*nsum+nsum,j*nsum:j*nsum+nsum].mean()
                out[i,j] = px
        return out
    
    def plot_parameter_maps(self,nsums, parn=1, cmap="jet", vmin = None, 
                            vmax = None, maxval = None):
        assert len(nsums)>=1
        assert len(nsums)<=5
        
        nr = 2
        nc = len(nsums)
        
        fig,axes = plt.subplots(nr,nc, sharex = "col", sharey = "col")
        
        parmaps = list()
        for j in range(len(nsums)):
            nsum = nsums[j]
            parmap = self.parfit_dict[nsum][:,:,parn].copy()
            if maxval is not None:
                parmap[parmap>maxval]=np.nan
            im=self.downsample_image(nsum)
            parmaps.append(parmap)
            ax0 = axes[0,j]
            ax1 = axes[1,j]
            extent = [0,im.shape[1]*nsum*self.xscale,0,
                      im.shape[0]*nsum*self.xscale]
            im0 = ax0.imshow(im,cmap="gray", extent = extent)
            ax0.set_title('Binning {}'.format(nsum))
            im1 = ax1.imshow(parmap,cmap=cmap, vmin = vmin, vmax = vmax,
                         extent = extent)
            ax1.set_xlabel("x [µm]")
            ax1.set_ylabel("y [µm]")
            
            fig.colorbar(im1,ax=ax1)
            fig.colorbar(im0,ax=ax0)
     
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
        
    def plot_curve(self, i0 = 0, j0 =0, nSum=1):
        correl = self.get_curve(i0 = i0, j0=j0,nSum=nSum)
        
        plt.figure()
        plt.semilogx(correl[:,0], correl[:,1])
        
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
            dmins = observation_sizes**2/(factor*self.dt*self.stack.shape[0]/100)
        for ns in nsums:
            ds = self.parfit_dict[ns][:,:,1]
            if self.threshold_map is not None:
                th = self.get_threshold_map(ns)
                ds = ds[th]
            ds_means.append(np.median(ds))
            ds_std.append( (np.percentile(ds,75)-np.percentile(ds,25))/2)
        ds_means = np.asarray(ds_means)
        ds_std = np.asarray(ds_std)
        
        plt.figure()
        if show_acceptable:
            plt.plot(nsums,dmins,color="gray")
            plt.plot(nsums,dmax,color="gray")
        plt.errorbar(nsums, ds_means, yerr=ds_std,capsize=5)
        plt.xlabel("Binning size")
        plt.ylabel("D (um2/s)")
        ymin = np.min(ds_means-ds_std)
        ymax = np.max(ds_means+ds_std)
        plt.ylim(bottom=ymin, top = ymax)
        return nsums, ds_means, ds_std

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
        return np.asarray(parameters),np.asarray(intensities)
    
    def plot_fits(self,nSum,maxcurves=None,dz=0.2):
        curves = self.correl_dicts[nSum]
        chisquares = self.chisquares_dict[nSum]
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
            chi = chisquares[j,k]
            yh = fits[j,k]
            a = corr[:3,1].mean()
            axes[0].semilogx(corr[:,0],corr[:,1]/a+dz*jj)
            axes[0].semilogx(corr[:yh.size,0],yh/a+dz*jj, color="k",linestyle="--")
            axes[1].semilogx(corr[:yh.size,0],yh/a-corr[:yh.size,1]/a+dz*jj, 
                             label = "curve ({},{}), chisquare {}".format(j,k, chi))
            jj+=1
        axes[1].legend()

    def plot_fits_ordered(self,nSum,maxcurves=None,dz=-0.2, order_dict = None):
        """Plots FCS curves ordered following values contained in a specific 
        dictionary"""
        assert order_dict is not None
        
        curves = self.correl_dicts[nSum].reshape(-1,*self.correl_dicts[nSum].shape[2:])
        parameters = order_dict[nSum].reshape(-1)
        fits = self.yh_dict[nSum].reshape(-1,*self.yh_dict[nSum].shape[2:])
        fig,axes = plt.subplots(1,2,sharex=True,sharey=True)
        jj = 0
        
        indices = np.arange(curves.shape[0])
        if maxcurves is not None and maxcurves<curves.shape[0]:
            indices = np.random.choice(curves.shape[0],maxcurves)
        
        indices = [x for _, x in sorted(zip(parameters[indices], indices))]
        
        for i in indices:
            corr = curves[i]
            param = parameters[i]
            yh = fits[i]
            a = yh[0]
            axes[0].semilogx(corr[:,0],corr[:,1]/a+dz*jj)
            axes[0].semilogx(corr[:yh.size,0],yh/a+dz*jj, color="k",linestyle="--")
            axes[1].semilogx(corr[:yh.size,0],yh/a-corr[:yh.size,1]/a+dz*jj, 
                             label = "corder parameter {}".format(param))
            jj+=1
        axes[1].legend()
        
    def plot_random_intensity(self, nSum = None):
        if nSum is None:
            nSum = min(self.traces_dict.keys())
        traces_arr = self.traces_dict[nSum]
        trace_raw = self.stack.mean(axis=(1,2))
        u,v = traces_arr.shape[:2]
        u1 = np.random.choice(u)
        v1 = np.random.choice(v)
        trace = traces_arr[u1,v1]
        trace_raw = self.stack[:,u1*nSum:(u1+1)*nSum,v1*nSum:(v1+1)*nSum].mean(axis=(1,2))
        
        plt.figure()
        plt.plot(trace, label = "Corrected")
        plt.axhline(trace.mean(),color="k",linestyle="--")
        plt.plot(trace_raw,label = "Raw average intensity")
        plt.xlabel("Time (frames)")
        plt.ylabel("Intensity")
        plt.legend()
        
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

        nsums = list(nsums)
        for j,ns in enumerate(nsums):
            ds = self.parfit_dict[ns][:,:,1]
            if self.threshold_map is not None:
                th = self.get_threshold_map(ns)
                ds = ds[th]
            if len(ds)==0:
                nsums.pop(j)
                continue
            taus = observation_sizes[j]**2/ds
            ds_means.append(np.median(taus))
            ds_std.append((np.percentile(taus,75)-np.percentile(taus,25)))
        ds_means=np.asarray(ds_means)
        
        plt.figure()
        if show_acceptable:
            dmax = np.ones_like(nsums)*(self.dt*10)
            dmins = np.ones_like(nsums)*self.dt*self.stack.shape[0]/100
            plt.plot(nsums,dmins,color="gray")
            plt.plot(nsums,dmax,color="gray")
        plt.errorbar(nsums, ds_means, yerr=ds_std,capsize=5)
        plt.xlabel("Binning size")
        plt.ylabel("tau (s)")
        return nsums, ds_means, ds_std
    