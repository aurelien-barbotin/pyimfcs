import matplotlib.pyplot as plt
import numpy as np
import math

plt.close('all')

plt.ion()
def multiplot_stack(stack,nsum, parn=1, normsize=2):
    
    mutable_object = {} 
    fig,axes = plt.subplots(2,3)
    axes=axes.ravel()
    
    def onclick(event):
        X_coordinate = event.xdata
        Y_coordinate = event.ydata
        mutable_object['click'] = X_coordinate
        
        i,j=math.floor(event.xdata+0.5), math.floor(event.ydata+0.5)
        trace = stack.traces_dict[nsum][j,i]
        line1.set_data([xt,trace])
        axes[2].set_ylim(bottom=trace.min()*0.8, top=trace.max()*1.2)
        
        axes[3].cla()
        axes[4].cla()
        ns, corrs1, fits1 = stack.get_acf_coord(nsum,j,i)
        for j in range(len(ns)):
            curve = corrs1[j]
            fits = fits1[j]
            axes[3].semilogx(curve[:,0], curve[:,1]/curve[:normsize,1].mean(),
                             label=ns[j],color="C{}".format(ns[j]))
            axes[3].semilogx(curve[:,0], fits/curve[:normsize,1].mean(), color="k",linestyle='--')
    
            axes[4].semilogx(curve[:,0], curve[:,1],
                             label=ns[j],color="C{}".format(ns[j]))
            axes[4].semilogx(curve[:,0], fits, color="k",linestyle='--')
            
        axes[3].legend()
        
        ns, dm, ds = stack.get_param_coord(nsum,j,i)
        axes[5].cla()
        axes[5].errorbar(ns,dm,yerr=ds,capsize=5)
        """ns, dm, ds = stack.get_param_coord(nsum,i,j)
        line4[0].set_data(ns,dm)
        line4[1].set_data(ns,ds)"""
        
        fig.canvas.draw_idle()
        
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    im = axes[0].imshow(stack.downsample_image(nsum))
    axes[0].set_title("Intensity")
    im2 = axes[1].imshow(stack.parfit_dict[nsum][:,:,parn])
    axes[1].set_title("Diffusion coeff.")
    fig.colorbar(im2,ax=axes[1])
    
    trace = stack.traces_dict[nsum][0,0]
    
    dt = stack.dt
    xt = np.arange(trace.size)*dt
    line1, = axes[2].plot(xt,trace)
    
    
    
multiplot_stack(stack,4)
