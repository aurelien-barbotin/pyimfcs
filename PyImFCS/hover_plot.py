import matplotlib.pyplot as plt
import numpy as np
import math

plt.close('all')

plt.ion()
def multiplot_stack(stack,nsum, parn=1):
    
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
        
        
        curve = stack.correl_dicts[nsum][j,i]
        fits = stack.yh_dict[nsum][j,i]
        line2.set_data([curve[:,0], curve[:,1]])
        line3.set_data(curve[:,0], fits)
        axes[3].set_ylim(bottom=curve[:,1].min()*0.95, 
                         top=curve[:,1].max()*1.05)
        
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
    
    curve = stack.correl_dicts[nsum][0,0]
    fits = stack.yh_dict[nsum][0,0]
    line2, = axes[3].semilogx(curve[:,0], curve[:,1])
    line3, = axes[3].semilogx(curve[:,0], fits, color="k",linestyle='--')
    
    sums = stack.correl_dicts.keys()
    sums = sorted([w for w in sums if w<=nsum])
        
    
    plt.show()
    
multiplot_stack(stack,4)
