#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 14:21:21 2023

@author: aurelienb
"""
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
ax.imshow(np.random.rand(10,10))

def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))

cid = fig.canvas.mpl_connect('button_press_event', onclick)


"""def onclick(event):
    if event.inaxes not in axes[0:2]:
        return
    X_coordinate = event.xdata
    mutable_object['click'] = X_coordinate
    
    ii0,jj0=math.floor(event.xdata+0.5), math.floor(event.ydata+0.5)"""