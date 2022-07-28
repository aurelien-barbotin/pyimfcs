#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:01:24 2022

@author: aurelien
"""

import math
import h5py


import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

from pyimfcs.io import get_dict
from pyimfcs.class_imFCS import new_chi_square

def plot_combined(combined,xname,measured,repeat, order = None,size = 8):
    """Plots results stored in a Dataframe using the SUperplots paradigm.
    
    Parameters:
        combined (Dataframe): dataframe containing merged results. 
        xname (str): name of conditions used in the x axis
        measured (str): name of the parameter to be measured (ex: D_[µm2/s])
        repeat (str): name of the repeat field
        order (list): list of condition names in the order the plots need to appear
        """
    if order is None:
        order = np.unique(combined[xname].values)
    # First measure
    ReplicateAverages = combined.groupby([xname,repeat], as_index=False).agg(
        {measured: "mean"})
    
    fig, ax = plt.subplots(1,1)
    
    sns.violinplot(x=xname, y=measured, data=combined, order = order,ax=ax,
                   color='gray',alpha=0.3)
    
    sns.swarmplot(x=xname, y=measured, hue=repeat, edgecolor="k", 
                       linewidth=2, data=ReplicateAverages, size=size, 
                       order = order,ax=ax)
    sns.swarmplot(x=xname, y=measured, hue=repeat, data=combined, order = order,ax=ax)
    ax.set_ylim(bottom=0)
    ax.legend_.remove()
    sns.reset_orig()
    

def superplot_files(files_list_list, conditions, nsum="nsum 3", keep_single_indices=True):
    
    if len(files_list_list)>1 and conditions is None:
        return ValueError('Please specify condition names')
    assert( conditions is None or len(conditions)==len(files_list_list))
    
    files = [w for flist in files_list_list for w in flist]
    excels = [pd.ExcelFile(w) for w in files]
    
    if conditions is not None:
        conditions_list = [conditions[j] for j, flist in 
                           enumerate(files_list_list) for w in flist]
    
    all_names = [w.sheet_names for w in excels]
    names0 = all_names[0]
    kept_names = list()
    for name in names0:
          if np.all([name in sublist for sublist in all_names]):
              kept_names.append(name)
    
    all_dfs = {}
    for name in kept_names:
        dfs = []
        maxindex = 0
        for j, xl in enumerate(excels):
            df = xl.parse(sheet_name=name, index_col = 0)
            if name=="parameters":
                fname = files[j]
                df["file"] = fname
                
            else:
                if keep_single_indices:
                    df['repeat']+=maxindex
                    maxindex = df['repeat'].values.max()+1
                    # print("maxindex", maxindex)
                else:
                    df['repeat'] = maxindex
                    maxindex += 1
            if conditions is not None:
                df["condition"] = conditions_list[j]
            dfs.append(df)
                
        dfs = pd.concat(dfs)
        all_dfs[name] = dfs
    xname = "condition"
    measured = "D [µm²/s]"
    plot_combined(all_dfs[nsum],xname,measured,'repeat')

plt.ion()
class FakeEvent():
    def __init__(self, ax):
        self.xdata = -0.45
        self.ydata = -0.45
        self.inaxes = ax
        
def multiplot_stack(stack,nsum, parn=1, normsize=1, fig = None, 
                    maxparval = None, chi_threshold = None):
    print('maxparval',maxparval)
    mutable_object = {}
    if fig is None:
        fig,axes = plt.subplots(2,4,figsize = (10,7))
    else:
        axes = []
        for j in range(8):
            axes.append(fig.add_subplot(2,4,j+1))
        axes = np.asarray(axes)
    axes=axes.ravel()
    axes[5:] = axes[5:][::-1]
    def onclick(event):
        if event.inaxes not in axes[0:2]:
            return
        X_coordinate = event.xdata
        Y_coordinate = event.ydata
        mutable_object['click'] = X_coordinate
        
        ii0,jj0=math.floor(event.xdata+0.5), math.floor(event.ydata+0.5)
        trace = stack.traces_dict[nsum][jj0,ii0]
        if stack.load_stack:
            trace_raw = stack.stack[:,jj0*nsum:jj0*nsum+nsum,
                                ii0*nsum:ii0*nsum+nsum].mean(axis=(1,2))
            trace_raw = trace_raw-trace_raw.min()+trace.max()
        else:
            trace_raw = trace
            
        line0.set_data(ii0,jj0)
        line1.set_data([xt2, trace_raw])
        line10.set_data([xt, trace])
        
        axes[2].set_ylim(bottom=trace.min()*0.95, top=trace_raw.max()*1.05)
        
        axes[3].cla()
        axes[4].cla()
        axes[7].cla()
        ns, corrs1, fits1 = stack.get_acf_coord(nsum,jj0,ii0, average=False)
        for k in range(len(ns)):
            curve = corrs1[k]
            fits = fits1[k]
            for u in range(curve.shape[0]):
                for v in range(curve.shape[1]):
                    cc = curve[u,v]
                    ff = fits[u,v]
                    axes[7].semilogx(cc[:,0], cc[:,1]/cc[:normsize,1].mean(),
                             label=ns[k],color="C{}".format(ns[k]))
                    axes[7].semilogx(cc[:,0], ff/cc[:normsize,1].mean(), color="k",linestyle='--')
            
            curve = corrs1[k].mean(axis=(0,1))
            fits = fits1[k].mean(axis=(0,1))
            axes[3].semilogx(curve[:,0], curve[:,1]/curve[:normsize,1].mean(),
                             label=ns[k],color="C{}".format(ns[k]))
            axes[3].semilogx(curve[:,0], fits/curve[:normsize,1].mean(), color="k",linestyle='--')
    
            axes[4].semilogx(curve[:,0], curve[:,1],
                             label=ns[k],color="C{}".format(ns[k]))
            axes[4].semilogx(curve[:,0], fits, color="k",linestyle='--')
            
        axes[3].legend()
        
        ns, dm, ds = stack.get_param_coord(nsum,jj0,ii0,parn=1)
        axes[5].cla()
        axes[5].errorbar(ns,dm,yerr=ds,capsize=5)
        
        ns, nnm, nns = stack.get_param_coord(nsum,jj0,ii0,parn=0)
        axes[6].cla()
        axes[6].errorbar(ns,nnm,yerr=nns,capsize=5)
            
        axes[2].set_title("Intensity timetrace")
        axes[2].set_xlabel("Time (frames)")
        axes[2].set_ylabel("Counts")
        
        axes[3].set_title("Averaged FCS curves (normalised)")
        axes[3].set_xlabel(r"$\rm \tau$")
        axes[3].set_ylabel(r"$\rm G(\tau)$")
        
        axes[4].set_title("Averaged FCS curves")
        axes[4].set_xlabel(r"$\rm \tau$")
        axes[4].set_ylabel(r"$\rm G(\tau)$")
        
        axes[7].set_title("All FCS curves")
        axes[7].set_xlabel(r"$\rm \tau$")
        axes[7].set_ylabel(r"$\rm G(\tau)$")
        
        axes[5].set_title("Diffusion coefficients")
        axes[5].set_xlabel("Binning (pixels)")
        axes[5].set_ylabel(r"$\rm D\ [\mu m^2/s]$")
        
        axes[6].set_title("Number of molecules")
        axes[6].set_xlabel("Binning (pixels)")
        axes[6].set_ylabel(r"$\rm N$")
        
        fig.canvas.draw_idle()
        
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    print(stack.traces_dict.keys(),nsum)
    image = stack.traces_dict[nsum].mean(axis=-1)


    im = axes[0].imshow(image)
    line0, = axes[0].plot(0,0,"x",color="red")
    axes[0].set_title("Intensity")
    fig.colorbar(im,ax=axes[0])
    
    dmap = stack.parfit_dict[nsum][:,:,parn].copy()
    dmap[dmap<0] = np.nan
    if maxparval is not None:
        dmap[dmap>maxparval] = np.nan
    
    if chi_threshold is not None:
        if len(stack.chisquares_dict)==0:
            stack.calculate_chisquares()
        chi_map = stack.chisquares_dict[nsum]
        dmap[chi_map>chi_threshold] = np.nan
        print('remove chis')
    im2 = axes[1].imshow(dmap)
    axes[1].set_title("Diffusion coeff.")
    fig.colorbar(im2,ax=axes[1])
    
    trace = stack.traces_dict[nsum][0,0]
    i,j=0,0
    if stack.load_stack:
        trace_raw = stack.stack[:,i*nsum:i*nsum+nsum,
                                j*nsum:j*nsum+nsum].mean(axis=(1,2))
    else:
        trace_raw = trace
    dt = stack.dt
    xt = np.arange(trace.size)
    xt2 = np.arange(trace_raw.size)
    
    line1, = axes[2].plot(xt,trace, label="Raw timetrace")
    line10, = axes[2].plot(xt2,trace_raw/trace_raw[0]*trace[0], label = "corrected")
    axes[2].legend()
    onclick(FakeEvent(axes[0]))
    fig.tight_layout()
    return onclick
      
def multiplot_stack_light(stack,nsum, parn=1, normsize=1, fig = None, 
                    maxparval = None, chi_threshold = None, intensity_threshold = None):
    """Light version of multiplot stack that does not load intensity traces"""
    mutable_object = {}
    if fig is None:
        fig,axes = plt.subplots(2,2,figsize = (10,7))
    else:
        axes = []
        for j in range(4):
            if j==1:
                print(axes)
                print(axes[0])
                axes.append(fig.add_subplot(2,2,j+1,sharex=axes[0],sharey=axes[0]))
            else:
                axes.append(fig.add_subplot(2,2,j+1))
        axes = np.asarray(axes)
    axes=axes.ravel()
    def onclick(event):
        if event.inaxes not in axes[0:2]:
            return
        X_coordinate = event.xdata
        Y_coordinate = event.ydata
        mutable_object['click'] = X_coordinate
        
        ii0,jj0=math.floor(event.xdata+0.5), math.floor(event.ydata+0.5)
        
        line0.set_data(ii0,jj0)
        line1.set_data(ii0,jj0)
        axes[2].cla()
        axes[3].cla()
        
        ns, corrs1, fits1 = stack.get_acf_coord(nsum,jj0,ii0, average=False)
        for k in range(len(ns)):
            curve = corrs1[k]
            fits = fits1[k]
            
            curve = corrs1[k].mean(axis=(0,1))
            fits = fits1[k].mean(axis=(0,1))
    
            axes[2].semilogx(curve[:,0], curve[:,1],
                             label=ns[k],color="C{}".format(ns[k]))
            axes[2].semilogx(curve[:,0], fits, color="k",linestyle='--')
            
        
        ns, dm, ds = stack.get_param_coord(nsum,jj0,ii0,parn=1)
        axes[3].cla()
        axes[3].errorbar(ns,dm,yerr=ds,capsize=5)
            
        
        axes[2].set_title("Averaged FCS curves")
        axes[2].set_xlabel(r"$\rm \tau$")
        axes[2].set_ylabel(r"$\rm G(\tau)$")
        
        axes[3].set_title("Diffusion coefficients")
        axes[3].set_xlabel("Binning (pixels)")
        axes[3].set_ylabel(r"$\rm D\ [\mu m^2/s]$")
        
        
        fig.canvas.draw_idle()
        
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    print(stack.traces_dict.keys(),nsum)
    image = stack.thumbnails_dict[nsum]


    im = axes[0].imshow(image)
    line0, = axes[0].plot(0,0,"x",color="red")
    axes[0].set_title("Intensity")
    fig.colorbar(im,ax=axes[0])
    
    dmap = stack.parfit_dict[nsum][:,:,parn].copy()
    dmap[dmap<0] = np.nan
    if maxparval is not None:
        dmap[dmap>maxparval] = np.nan
    
    if chi_threshold is not None:
        if len(stack.chisquares_dict)==0:
            stack.calculate_chisquares()
        chi_map = stack.chisquares_dict[nsum]
        dmap[chi_map>chi_threshold] = np.nan
    
    if intensity_threshold is not None:
        ithr = intensity_threshold*(np.percentile(image,98)- np.percentile(image,2)) + np.percentile(image,2)
        print("Intensity threshold: max, min, treshold",image.max(),image.min(),ithr)
        mask_int = (image>ithr).astype(float)
        axes[0].contour(mask_int,levels=[0.5], colors="red")
        
    im2 = axes[1].imshow(dmap)
    axes[1].set_title("Diffusion coeff.")
    line1, = axes[1].plot(0,0,"x",color="red")
    fig.colorbar(im2,ax=axes[1])

    axes[2].legend()
    onclick(FakeEvent(axes[0]))
    fig.tight_layout()
    return onclick

markers = ["o","v","s","x","1"]
def findindex(x,y,xy):
    i0 = np.abs(x-xy[0])
    i1 = np.abs(y-xy[1])
    assert(np.isclose(i0.min(),0))
    assert(np.isclose(i1.min(),0))
    
    ind0 = np.argmin(i0)
    ind1 = np.argmin(i1)
    if ind0!=ind1:
        print('values: {},{}'.format(xy[0],xy[1]))
        print(ind0,ind1)
        ind0 = np.argmin(i0+i1)
        print("new ind0: {}".format(ind0))
    # assert(ind0==ind1)
    return ind0

def hover_plot(x,y,curves_subsets, fits_subsets,labels, xlabel = 'chi', 
               ylabel = 'new chi', xlog = False, ylog = False):
    """Plots a hoverable plot that displays FCS curves.
    Parameters:
        x (list): list of value arrays for different conditions
        y (list): list of value arrays for different conditions
        curves_subsets (list)"""
    #curves = np.concatenate(curves_subsets,axis=0)
    curves = np.array([x for w in curves_subsets for x in w])
    fits = np.array([x for w in fits_subsets for x in w])
    predictions = [ [labels[w] for uu in curves_subsets[w]] for w in range(len(labels))]
    predictions = np.array(predictions)
    predictions = np.concatenate(predictions,axis=0)
    
    fig,axes = plt.subplots(1,2)
    ax = axes[0]
    xx = np.concatenate(x,axis=0)
    yy = np.concatenate(y,axis=0)
    
    msk = np.logical_and(~np.isnan(xx), ~np.isnan(yy))
    xx = xx[msk]
    yy=yy[msk]
    predictions = predictions[msk]
    curves = curves[msk]
    fits = fits[msk]
    
    sc = ax.scatter(xx, yy, s=np.ones_like(xx))
    
    for j in range(len(labels)):
        xs, ys = x[j],y[j]
        ax.scatter(xs, ys, marker = markers[j%len(markers)],label=labels[j])
        
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    
    ax1 = axes[1]
    ax1.set_xlabel("Frame")
    ax1.set_ylabel("Intensity")
    ax1.set_ylim(top=1)
    axes[1].plot(np.random.rand(50))
    
    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)
    
    def update_annot(ind):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        corrind = findindex(xx, yy, pos)
        xcurve = curves[corrind][:,0]
        ycurve = curves[corrind][:,1]
        ax1.semilogx(xcurve, ycurve/ycurve.max())
        ax1.semilogx(xcurve,fits[corrind]/ycurve.max(), color="k",linestyle="--")
        ax1.set_title(str(predictions[corrind])+" index "+str(corrind))
        ax1.set_ylim(top=1)
        
    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                ax1.cla()
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
                
            else:
                ax1.cla()
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()
    
    fig.canvas.mpl_connect("motion_notify_event", hover)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    plt.show()
    return fig,axes

def plot_fit_error(files,nsums = None, intensity_threshold = None, chi_threshold = None):
    all_curves = []
    all_fits = []
    all_chis = []
    all_chis_new = []
    all_diffs = []
    all_labels = []
    
    diffs_out = list()
    chis_out = list()
    
    # if nsums is not specified; take the nsums in every file and check that 
    # it is consistent across the dataset
    check_nsums = False
    nsums_backup = None
    if nsums is None:
        check_nsums = True
    
    for k,file in enumerate(files):
        h5f = h5py.File(file,mode='r')
        dict_diffcoeff = get_dict(h5f, "parameters_fits")
        dict_curves = get_dict(h5f, "correlations")
        dict_curves_fits = get_dict(h5f, "yhat")
        dict_traces = get_dict(h5f,'traces')
        h5f.close()
        
        if check_nsums:
            nsums = list(dict_curves.keys())
            if nsums_backup is None:
                nsums_backup = nsums
            if sorted(nsums)!=sorted(nsums_backup):
                raise KeyError('Not all files were processed with the same parameters')
            
        diffs_out.append( dict(zip(nsums,[[] for w in nsums])))
        chis_out.append(dict(zip(nsums,[[] for w in nsums])))
        
        for jj, nsum in enumerate(nsums):
            diffcoeffs = dict_diffcoeff[nsum][:,:,1]
            curves = dict_curves[nsum]
            curves_fits = dict_curves_fits[nsum]
            traces = dict_traces[nsum]
            
            intensities = traces.mean(axis=2).reshape(-1)
            # xcurves = curves[:,:,:,0]
            
            ycurves = curves[:,:,:,1]
            
            chis = np.sqrt(((ycurves-curves_fits)**2).sum(axis=2) )
            
            chis_new = np.zeros_like(chis)
            
            for i in range(chis_new.shape[0]):
                for j in range(chis_new.shape[1]):
                    
                    chis_new[i,j] = new_chi_square(ycurves[i,j], curves_fits[i,j])
        
            curves_reshaped = curves.reshape((curves.shape[0]*curves.shape[1],curves.shape[2],2))
            fits_reshaped = curves_fits.reshape((curves.shape[0]*curves.shape[1],curves.shape[2]))
            
            msk = np.ones_like(intensities, dtype = bool)
            msk[diffcoeffs.reshape(-1)<0] = 0
            if intensity_threshold is not None:
                msk = np.logical_and(msk,
                                     intensities>(intensities.max()*intensity_threshold))
                
            chis_new = chis_new.reshape(-1)[msk]
            curves_reshaped = curves_reshaped[msk]
            fits_reshaped = fits_reshaped[msk]
            diffs = diffcoeffs.reshape(-1)[msk]
            
            all_curves.append(curves_reshaped)
            all_fits.append(fits_reshaped)
            all_chis.append(chis.reshape(-1)[msk])
            all_chis_new.append(chis_new)
            all_diffs.append(diffs)
            # all_intensities.append()
            all_labels.append(file.split('/')[-1]+"_{}".format(nsum))
            
            diffs_out[k][nsum] = diffs[chis_new<chi_threshold]
            chis_out[k][nsum]= chis_new[chis_new<chi_threshold]
    
    fig,axes= hover_plot(all_chis_new,all_diffs,
               all_curves,all_fits,all_labels, xlabel = 'chinew',ylabel = "D",ylog = True, xlog = True)
    axes[0].axvline(chi_threshold,color="k")
    
    return diffs_out, chis_out

def plot_diffusion_map(file, nsum = 2, intensity_threshold = 0.4, 
                       chi_threshold = 0.02, debug_plot=True):
    h5f = h5py.File(file,mode='r')
    dict_diffcoeff = get_dict(h5f, "parameters_fits")
    dict_curves = get_dict(h5f, "correlations")
    dict_curves_fits = get_dict(h5f, "yhat")
    dict_traces = get_dict(h5f,'traces')
    xscale = h5f['parameters/xscale'][()]
    h5f.close()

    diffcoeffs = dict_diffcoeff[nsum][:,:,1]
    curves = dict_curves[nsum]
    curves_fits = dict_curves_fits[nsum]
    traces = dict_traces[nsum]
    
    intensities = traces.mean(axis=2)
    # xcurves = curves[:,:,:,0]
    
    ycurves = curves[:,:,:,1]
    
    chis = np.sqrt(((ycurves-curves_fits)**2).sum(axis=2) )
    
    chis_new = np.zeros_like(chis)
    
    for i in range(chis_new.shape[0]):
        for j in range(chis_new.shape[1]):
            
            chis_new[i,j] = new_chi_square(ycurves[i,j], curves_fits[i,j])

    
    msk0 = np.ones_like(intensities, dtype = bool)
    
    if intensity_threshold is not None:
        msk0 = intensities>(intensities.max()*intensity_threshold)
    msk = np.logical_and(msk0,chis_new<chi_threshold)
    dmap = diffcoeffs.copy()
    dmap[~msk] = np.nan
    
    if debug_plot:
        plt.figure()
        plt.subplot(221)
        plt.imshow(msk)
        plt.title('Final mask')
        
        plt.subplot(222)
        plt.imshow(msk0)
        plt.title('Intensities mask')
        
        plt.subplot(223)
        plt.imshow(chis_new)
        plt.title('Chis')
        
        plt.subplot(224)
        plt.imshow(diffcoeffs)
    
    extent = [0, intensities.shape[1]*xscale*nsum,0,intensities.shape[0]*xscale*nsum]
    plt.figure()
    plt.subplot(121)
    plt.imshow(intensities, extent=extent)
    plt.title('Intensity projection')
    plt.xlabel('x [µm]')
    plt.ylabel('y [µm]')
    plt.colorbar()
    
    plt.subplot(122)
    plt.imshow(dmap, extent=extent)
    plt.xlabel('x [µm]')
    plt.ylabel('y [µm]')
    plt.title('Diffusion map')
    cbar = plt.colorbar()
    cbar.set_label('D [µm²/s]')
    return dmap

def interactive_plot_h5(stack, fig = None, nsum = 2, vmax=  None, 
                        chi_threshold = None, light_version = False, intensity_threshold = None):
    """Wrapper function to plot imFCS results from GUI"""
    if fig is None:
        print('creating figure')
        fig  = plt.subplots(2,4,figsize = (10,7))
    if light_version:
        onclickfunction = multiplot_stack_light(stack,nsum, fig=fig, maxparval = vmax, 
                                          chi_threshold=chi_threshold, 
                                          intensity_threshold = intensity_threshold)
    else:
        onclickfunction = multiplot_stack(stack,nsum, fig=fig, maxparval = vmax, 
                                          chi_threshold=chi_threshold)
    return onclickfunction

if __name__=='__main__':
    import glob

    files_veg = glob.glob("/home/aurelienb/Documents/Projects/imFCS/Subtilis_veg/TFSM/*.xlsx")

    files_20 = [w for w in files_veg if '05_12' in w]
    files_37 = [w for w in files_veg if '05_12' not in w]
    files = [files_37, files_20]
    print(files)
    conditions = ["grown 37°C", "grown 20°C"]
    superplot_files(files,conditions)