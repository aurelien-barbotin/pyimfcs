import matplotlib.pyplot as plt
import numpy as np
import math
import h5py

from PyImFCS.io import get_dict
from PyImFCS.class_imFCS import new_chi_square

plt.ion()
class FakeEvent():
    def __init__(self, ax):
        self.xdata = -0.45
        self.ydata = -0.45
        self.inaxes = ax
        
def multiplot_stack(stack,nsum, parn=1, normsize=1, fig = None, maxparval = None):
    
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
    
    dmap = stack.parfit_dict[nsum][:,:,parn]
    if maxparval is not None:
        dmap[dmap>maxparval] = np.nan
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
    xt = np.arange(trace.size)*dt
    xt2 = np.arange(trace_raw.size)*dt
    
    line1, = axes[2].plot(xt,trace, label="Raw timetrace")
    line10, = axes[2].plot(xt2,trace_raw/trace_raw[0]*trace[0], label = "corrected")
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
    curves = [x for w in curves_subsets for x in w]
    fits = [x for w in fits_subsets for x in w]
    predictions = [ [labels[w] for uu in curves_subsets[w]] for w in range(len(labels))]
    predictions = np.array(predictions)
    predictions = np.concatenate(predictions,axis=0)
    
    fig,axes = plt.subplots(1,2)
    ax = axes[0]
    xx = np.concatenate(x,axis=0)
    yy = np.concatenate(y,axis=0)
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

def get_fit_error(files,nsums, intensity_threshold = None):
    all_curves = []
    all_fits = []
    all_chis = []
    all_chis_new = []
    all_diffs = []
    all_intensities = []
    all_labels = []
    
    for file in files:
        h5f = h5py.File(file,mode='r')
        dict_diffcoeff = get_dict(h5f, "parameters_fits")
        dict_curves = get_dict(h5f, "correlations")
        dict_curves_fits = get_dict(h5f, "yhat")
        dict_traces = get_dict(h5f,'traces')
        h5f.close()
        
        for nsum in nsums:
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
            if intensity_threshold is not None:
                msk = intensities>(intensities.max()*intensity_threshold)
                curves_reshaped = curves_reshaped[msk]
                fits_reshaped = fits_reshaped[msk]

            all_curves.append(curves_reshaped)
            all_fits.append(fits_reshaped)
            all_chis.append(chis.reshape(-1)[msk])
            all_chis_new.append(chis_new.reshape(-1)[msk])
            all_diffs.append(diffcoeffs.reshape(-1)[msk])
            # all_intensities.append()
            all_labels.append(file.split('/')[-1]+"_{}".format(nsum))
            
    
    
    """hover_plot(all_chis,all_chis_new,
               all_curves,all_fits,all_labels)"""
    
    
    hover_plot(all_chis_new,all_diffs,
               all_curves,all_fits,all_labels, xlabel = 'chinew',ylabel = "D",ylog = True, xlog = True)


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
        plt.imshow(chis_new<chi_threshold)
        plt.title('Chis mask')
        
        plt.subplot(224)
        plt.imshow(diffcoeffs)
    
    extent = [0, intensities.shape[1]*xscale*nsum,0,intensities.shape[0]*xscale*nsum]
    plt.figure()
    plt.subplot(121)
    plt.imshow(intensities, extent=extent)
    plt.title('Intensity projection')
    plt.xlabel('x [µm]')
    plt.xlabel('y [µm]')
    plt.colorbar()
    
    plt.subplot(122)
    plt.imshow(dmap, extent=extent)
    plt.xlabel('x [µm]')
    plt.xlabel('y [µm]')
    plt.title('Diffusion map')
    cbar = plt.colorbar()
    cbar.set_label('D [µm²/s]')