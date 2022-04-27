#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 16:06:50 2021

@author: aurelien
"""
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    from SO: https://stackoverflow.com/questions/13685386/
    matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    
def sphere2ellipsoid(coords,b):
    a = 1/np.sqrt(b)
    mat = np.eye(3)*a
    mat[-1,-1] = b
    new_coords = np.dot(mat,coords)
    return new_coords

def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(xarray[:,:,frame_number], 
                              yarray[:,:,frame_number], 
                              zarray[:,:,frame_number], cmap="magma")
    
dt = 10**-3
nsteps = 300

motion_period = 50*dt
motion_freq = 2*np.pi/motion_period
motion_amplitude = 0.1
motion_z = 1-motion_amplitude + np.cos(motion_freq*np.arange(nsteps)*dt)*motion_amplitude

R=7
npts = 200
x_surface = np.concatenate([np.linspace(-R,R,npts)])
y_surface = np.concatenate([np.linspace(-R,R,npts)])

x_surface,y_surface = np.meshgrid(x_surface,y_surface)
"""r2 = x_surface**2+y_surface**2
x_surface = x_surface[r2<R**2]
y_surface = y_surface[r2<R**2]"""
z_surface = np.sqrt(R**2-(x_surface**2+y_surface**2))

z_surface[np.isnan(z_surface)]=0
#z_surface[]*=-1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
"""sf = ax.plot_surface(x_surface, y_surface, z_surface,
                   linewidth=0, antialiased=False)"""

# set_axes_equal(ax)

zarray = np.zeros([z_surface.shape[0], z_surface.shape[1], nsteps])
xarray = np.zeros_like(zarray)
yarray= np.zeros_like(zarray)

for i in range(nsteps):
    
    aa = sphere2ellipsoid(np.array([x_surface.reshape(-1),y_surface.reshape(-1),
                               z_surface.reshape(-1)]), motion_z[i])
    x1,y1,z1 = aa[0].reshape(npts,npts), aa[1].reshape(npts,npts), aa[2].reshape(npts,npts)
    zarray[:,:,i] = z1
    xarray[:,:,i] = x1
    yarray[:,:,i] = y1
    
zarray*=-1
zarray-=zarray.min(axis=(0,1))[np.newaxis,np.newaxis,:]

x,y = x_surface, y_surface
plot = [ax.plot_surface(x, y, zarray[:,:,0], color='0.75', rstride=1, cstride=1)]

ax.set_zlim(zarray.min(),zarray.max())
ax.set_xlim(xarray.min(),xarray.max())
ax.set_ylim(yarray.min(),yarray.max())

animate = animation.FuncAnimation(fig, update_plot, nsteps, fargs=(zarray, plot))
plt.show()
