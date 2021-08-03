import matplotlib.pyplot as plt
import numpy as np
import math

plt.close('all')

plt.ion()
class FakeEvent():
    def __init__(self, ax):
        self.xdata = -0.45
        self.ydata = -0.45
        self.inaxes = ax
        
def multiplot_stack(stack,nsum, parn=1, normsize=2):
    
    mutable_object = {} 
    fig,axes = plt.subplots(2,4,figsize = (10,7))
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
        trace_raw = stack.stack[:,jj0*nsum:jj0*nsum+nsum,
                            ii0*nsum:ii0*nsum+nsum].mean(axis=(1,2))
        trace_raw = trace_raw-trace_raw.min()+trace.max()
        print('trace raw',trace_raw.shape,trace.shape, xt.shape, xt2.shape)
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
        
        axes[3].set_title("FCS curves (normalised)")
        axes[3].set_xlabel(r"$\rm \tau$")
        axes[3].set_ylabel(r"$\rm G(\tau)$")
        
        axes[4].set_title("FCS curves")
        axes[4].set_xlabel(r"$\rm \tau$")
        axes[4].set_ylabel(r"$\rm G(\tau)$")
        
        axes[5].set_xlabel("Binning (pixels)")
        axes[5].set_ylabel(r"$\rm D\ [\mu m^2/s]$")
        
        axes[6].set_xlabel("Binning (pixels)")
        axes[6].set_ylabel(r"$\rm N$")
        
        fig.canvas.draw_idle()
        
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    im = axes[0].imshow(stack.downsample_image(nsum))
    line0, = axes[0].plot(0,0,"x",color="red")
    axes[0].set_title("Intensity")
    im2 = axes[1].imshow(stack.parfit_dict[nsum][:,:,parn])
    axes[1].set_title("Diffusion coeff.")
    fig.colorbar(im2,ax=axes[1])
    
    trace = stack.traces_dict[nsum][0,0]
    i,j=0,0
    trace_raw = stack.stack[:,i*nsum:i*nsum+nsum,
                            j*nsum:j*nsum+nsum].mean(axis=(1,2))
    dt = stack.dt
    xt = np.arange(trace.size)*dt
    xt2 = np.arange(trace_raw.size)*dt
    print(trace_raw.shape,trace.shape,xt.shape,xt2.shape)
    line1, = axes[2].plot(xt,trace, label="Raw timetrace")
    line10, = axes[2].plot(xt2,trace_raw/trace_raw[0]*trace[0], label = "corrected")
    axes[2].legend()
    onclick(FakeEvent(axes[0]))
    fig.tight_layout()

# multiplot_stack(stack,6)
