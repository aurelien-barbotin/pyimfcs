#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 14:18:15 2023

@author: aurelienb
"""

import matplotlib.pyplot as plt
import math
import numpy as np
from pyimfcs.methods import indices_intensity_threshold
from pyimfcs.plotting import FakeEvent

def interactive_fcs_plot_2components(stack,nsum, parn=1, normsize=1, fig = None, 
                    vmax = None, vmin = None, chi_threshold = None, 
                    ith = None, light_version = True,use_mask=False):
    """Single method for display of FCS results."""
    mutable_object = {}
    
    nlines = 2
    ncols = 3
    if light_version:
        ncols = 3
    if fig is None:
        fig,axes = plt.subplots(2,ncols,figsize = (10,7))
    else:
        axes = []
        for j in range(nlines*ncols):
            if j==1:
                axes.append(fig.add_subplot(nlines,ncols,j+1,sharex=axes[0],sharey=axes[0]))
            else:
                axes.append(fig.add_subplot(nlines,ncols,j+1))
        axes = np.asarray(axes)
    axes=axes.ravel(order='C')
    def onclick(event):
        if event.inaxes not in axes[0:4]:
            return
        X_coordinate = event.xdata
        mutable_object['click'] = X_coordinate
        
        ii0,jj0=math.floor(event.xdata+0.5), math.floor(event.ydata+0.5)
        
        line0.set_data(ii0,jj0)
        line1.set_data(ii0,jj0)
        line2.set_data(ii0,jj0)
        line3.set_data(ii0,jj0)
        axes[4].cla()
        
        barvals=np.array([dmap[jj0,ii0],dmap2[jj0,ii0]])
        barvals[np.isnan(barvals)]=0
        axes[4].bar(['slow','fast'],barvals)
        """if not light_version:
            axes[4].cla()
            axes[5].cla()
        
            ii0,jj0=math.floor(event.xdata+0.5), math.floor(event.ydata+0.5)
            trace = stack.traces_dict[nsum][jj0,ii0]
            if stack.load_stack:
                trace_raw = stack.stack[:,jj0*nsum:jj0*nsum+nsum,
                                    ii0*nsum:ii0*nsum+nsum].mean(axis=(1,2))
                trace_raw = trace_raw-trace_raw.min()+trace.max()
            else:
                trace_raw = trace
                
            line0.set_data(ii0,jj0)
            line_raw.set_data([xt2, trace_raw])
            line_corrected.set_data([xt, trace])
            
            axes[7].set_ylim(bottom=trace.min()*0.95, top=trace_raw.max()*1.05)"""
        # populates this axis with the FCS curve clicked on
        axes[5].cla()
        ns, corrs1, fits1 = stack.get_acf_coord(nsum,jj0,ii0, average=False)
        for k in range(len(ns)):
            curve = corrs1[k]
            fits = fits1[k]
            
            curve = corrs1[k].mean(axis=(0,1))
            fits = fits1[k].mean(axis=(0,1))
            a0 = fits[0]
            if np.isclose(a0,0):
                a0=1
            # shows only one curve
            if k==len(ns)-1:
                axes[5].semilogx(curve[:,0], curve[:,1],
                                 label=ns[k],color="C{}".format(ns[k]))
                axes[5].semilogx(curve[:,0], fits, color="k",linestyle='--')
            
        ns, dm, ds = stack.get_param_coord(nsum,jj0,ii0,parn=parn)
        
        axes[5].set_title("FCS curve")
        axes[5].set_xlabel(r"$\rm \tau\ [s]$")
        axes[5].set_ylabel(r"$\rm G(\tau)$")
        
        axes[4].set_title("Diffusion coefficients")
        """if not light_version:
            ns, nnm, nns = stack.get_param_coord(nsum,jj0,ii0,parn=0)
            axes[6].cla()
            axes[6].errorbar(ns,nnm,yerr=nns,capsize=5)
            axes[4].set_title("All FCS curves (normalised)")
            axes[4].set_xlabel(r"$\rm \tau$")
            axes[4].set_ylabel(r"$\rm G(\tau)$")
            
            axes[7].set_title("Intensity timetrace")
            axes[7].set_xlabel("Time (frames)")
            axes[7].set_ylabel("Counts")
        
            axes[5].set_title("Averaged FCS curves")
            axes[5].set_xlabel(r"$\rm \tau$")
            axes[5].set_ylabel(r"$\rm G(\tau)$")
            axes[5].legend()
            
            axes[6].set_title("Number of molecules")
            axes[6].set_xlabel("Binning (pixels)")
            axes[6].set_ylabel(r"$\rm N$")"""
        
        
        fig.canvas.draw_idle()
        
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    image = stack.thumbnails_dict[nsum]


    im = axes[0].imshow(image,cmap="gray")
    line0, = axes[0].plot(0,0,"x",color="red")
    axes[0].set_title("Intensity")
    fig.colorbar(im,ax=axes[0])
    
    dmap = stack.fit_results_dict[nsum][:,:,parn].copy()
    dmap=dmap.astype(float)
    dmap[dmap<0] = np.nan
    
    dmap2 = stack.fit_results_dict[nsum][:,:,2].copy()
    dmap2=dmap2.astype(float)
    dmap2[dmap2<0] = np.nan
    
    ratio_map = stack.fit_results_dict[nsum][:,:,3].copy()
    
    # !! TOUPDATE
    if vmax is not None:
        dmap2[mask_outlier_fraction(dmap2,outlier_factor=vmax)]=np.nan
    if vmin is not None:
        dmap[mask_outlier_fraction(dmap,outlier_factor=vmin)]=np.nan
        
    if chi_threshold is not None:
        if len(stack.chisquares_dict)==0:
            stack.calculate_chisquares()
        chi_map = stack.chisquares_dict[nsum]
        
        dmap[chi_map>chi_threshold] = np.nan
        dmap2[chi_map>chi_threshold] = np.nan
        ratio_map[chi_map>chi_threshold] = np.nan
        
    from pyimfcs.metrics import intensity_threshold
    from pyimfcs.methods import downsample_mask
    if ith is not None and use_mask and stack.mask is not None:
        mask = downsample_mask(stack.mask, nsum)
        labels = indices_intensity_threshold(mask,ith,image)
        for lab in np.unique(labels):
            if lab!=0:
                mask_int=labels==lab
                axes[0].contour(mask_int,levels=[0.5], colors="C{}".format(int(lab)))
    elif ith is not None:
        ithr = intensity_threshold(ith,image)
        mask_int = (image>ithr).astype(float)
        axes[0].contour(mask_int,levels=[0.5], colors="red")
    
    im2 = axes[1].imshow(dmap)
    axes[1].set_title("Diffusion Map (slow)")
    line1, = axes[1].plot(0,0,"x",color="red")
    fig.colorbar(im2,ax=axes[1])


    im3 = axes[2].imshow(dmap2)
    axes[2].set_title("Diffusion Map (fast)")
    line2, = axes[2].plot(0,0,"x",color="red")
    fig.colorbar(im3,ax=axes[2])
    
    im4 = axes[3].imshow(ratio_map,vmin=0,vmax=1,cmap="RdYlGn")
    axes[3].set_title("Ratio slow/fast")
    line3, = axes[3].plot(0,0,"x",color="red")
    fig.colorbar(im4,ax=axes[3])
    
    axes[5].legend()
    
    # intensity traces
    """if not light_version:
        trace = stack.traces_dict[nsum][0,0]
        i,j=0,0
        if stack.load_stack:
            trace_raw = stack.stack[:,i*nsum:i*nsum+nsum,
                                    j*nsum:j*nsum+nsum].mean(axis=(1,2))
        else:
            trace_raw = trace
        xt = np.arange(trace.size)
        xt2 = np.arange(trace_raw.size)
        
        line_raw, = axes[7].plot(xt,trace, label="Raw timetrace")
        line_corrected, = axes[7].plot(xt2,trace_raw/trace_raw[0]*trace[0], label = "corrected")
        axes[7].legend()"""
    onclick(FakeEvent(axes[0]))
    fig.tight_layout()
    return onclick

def mask_outlier_fraction(img,outlier_factor=5):
    msk=~np.isnan(img)
    mean = np.median(img[msk])
    print(mean,np.percentile(img[msk],75))
    print(np.count_nonzero(msk))
    msk = np.logical_or(img > (mean+outlier_factor*np.percentile(img[msk],75)),
                         img < (mean-outlier_factor*np.percentile(img[msk],25)) )
    return msk

from pyimfcs.class_imFCS import StackFCS
path = "/run/user/1001/gvfs/microscopy/ZEISS/Aurelien/2023_09_06/S2_14h15/deltaMFD/Image 95.h5"

stack=StackFCS(path,load_stack=False)
stack.load()
stack.calculate_chisquares()
nsum = 2
chi_threshold = 0.03
outlier_factor = 3

dmap1 = stack.fit_results_dict[nsum][:,:,1].copy()
dmap2 = stack.fit_results_dict[nsum][:,:,2].copy()
ratios = stack.fit_results_dict[nsum][:,:,3].copy()

chi_map = stack.chisquares_dict[nsum]

dmap1[chi_map>chi_threshold] = np.nan
dmap1[mask_outlier_fraction(dmap1)]=np.nan

dmap2[chi_map>chi_threshold] = np.nan
dmap2[mask_outlier_fraction(dmap2)]=np.nan

ratios[chi_map>chi_threshold] = np.nan

plt.figure()
plt.subplot(311)
plt.imshow(dmap1)
plt.colorbar()
plt.subplot(312)
plt.imshow(dmap2)
plt.colorbar()
plt.subplot(313)
plt.imshow(ratios,vmin=0,vmax=1,cmap="RdYlGn")
plt.colorbar()

out = interactive_fcs_plot_2components(stack, 2,vmin=3,vmax=3,ith=0.8,chi_threshold=0.03)

results=stack.extract_results()

results_ns2 = dict(zip(results.keys(),[results[w][2] for w in results.keys()]))
import pandas as pd

df = pd.DataFrame.from_dict(results_ns2)
