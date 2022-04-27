# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 12:09:25 2021

@author: abarbotin
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import gaussian
import tifffile

plt.close('all')

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
    
def gauss(x,sig):
    return np.exp(-2*x**2/(sig**2))

def spherical2cart(R,phi,theta):
    # theta+=np.pi/3
    x = R*np.cos(phi)*np.sin(theta)
    y = R*np.sin(phi)*np.sin(theta)
    z = R*np.cos(theta)
    return np.array([x, y, z])

def cart2spherical(x,y,z):
    # theta+=np.pi/3
    theta = np.arccos(z)
    phi = np.arctan2(y,x)
    return np.array([phi, theta])

def spherical2cart_special(R,theta,phi):
    # theta+=np.pi/3
    x = R*np.cos(theta)
    if (1-np.cos(theta)**2-np.cos(phi)**2<0).any():
        print('negative value')
        vc = (1-np.cos(theta)**2-np.cos(phi)**2<0)
        print(np.count_nonzero(vc))
        sgns = np.sign(np.sin(theta)*np.sin(phi))
        print(np.count_nonzero(sgns<0),np.count_nonzero(sgns>0))
    y = R*np.sqrt(1-np.cos(theta)**2-np.cos(phi)**2)*np.sign(np.sin(theta)*np.sin(phi))
    z = R*np.cos(phi)+R
    return x, y, z

def sphere2ellipsoid(coords,b):
    a = 1/np.sqrt(b)
    mat = np.eye(3)*a
    mat[-1,-1] = b
    new_coords = np.dot(mat,coords)
    return new_coords

plot = True
save = True

psize = 0.16
sigma_psf = 0.2/psize
sigmaz = 4*sigma_psf
dz_tirf = 0.2 # um

dt = 1*10**-3 # s
D = 3 #um2/s

R = 8 #um
brightness = 18*10**3 #Hz/molecule

nsteps = 20000
nparts = 1000

# ------ GUV compression ----------------------

motion_period = 3500*dt
motion_freq = 2*np.pi/motion_period
motion_amplitude = 0.05
motion_z = 1-motion_amplitude + np.cos(motion_freq*np.arange(nsteps)*dt)*motion_amplitude

plt.figure()
plt.plot(motion_z)
plt.ylim(bottom=0)
pos0 = np.random.uniform(size = (nparts,2))*2*np.pi

# ---------- Calculation of positions-------------------
moves = np.random.normal(scale = np.sqrt(2*D*dt)/R,
                         size = (nsteps,nparts,2) )

moves[0] = 0

positions = np.cumsum(moves,axis=0)
positions = positions+pos0[np.newaxis,:,:]


p1 = pos0[:,0]
phis = moves[:,:,0]/np.sin(positions[:,:,1])
positions[:,:,0] = p1[np.newaxis,:] + np.cumsum(phis,axis=0)

positions = positions%(2*np.pi)
def cartesian_sampling():
    pos0_cart = spherical2cart(1,pos0[:,0], pos0[:,1])
    positions_cart = np.zeros((nsteps,3,nparts))
    positions_cart[0] = pos0_cart
    angles = [pos0.T]
    for j in range(nsteps-1):
        curpos_cart = positions_cart[j]
        curpos_angle = cart2spherical(curpos_cart[0], curpos_cart[1], curpos_cart[2])
        mv = moves[j].T
        mv[0]/=np.sin(curpos_angle[1])
        newpos_angle = curpos_angle+mv
        angles.append(newpos_angle)
        newpos_cart = spherical2cart(1, newpos_angle[0],  newpos_angle[1])
        positions_cart[j+1] = newpos_cart
        if j%500==0:
            print('Calculating position {}'.format(j))
    positions_cart = positions_cart.swapaxes(1,2)
    angles = np.asarray(angles)
    
    # x,y,z = spherical2cart(R,positions[:,:,0],positions[:,:,1])
    x,y,z = positions_cart[:,:,0], positions_cart[:,:,1], positions_cart[:,:,2]
    x*=R
    y*=R
    z*=R
    
    if plot:
        nn0 = 50
        plt.figure()
        plt.subplot(121)
        plt.plot(x[:nn0,0],"-o", label = "x")
        plt.plot(y[:nn0,0],"-o", label = "y")
        plt.plot(z[:nn0,0],"-o", label = "z")
        plt.legend()
        plt.subplot(122)
        plt.plot(angles[:nn0,0,50],"-o", label = "theta")
        plt.plot(angles[:nn0,1,50],"-o", label = "phi")
        print(angles[0,0,50]-angles[1,0,50], angles[0,1,50]-angles[1,1,50])
        plt.legend()
    return x,y,z

z,y,x = spherical2cart(R, positions[:,:,0], positions[:,:,1])
# x,z,y =  cartesian_sampling()
"""
coords = np.array([x,y,z])
aa = 1.2
coords_new=np.zeros_like(coords)
for j in range(coords_new.shape[2]):
    coords_new[:,:,j] = sphere2ellipsoid(coords[:,:,j], aa)
x,y,z = coords_new
z+=R/(aa**2)"""





if plot:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for j in range(nparts):
        ax.plot(x[::200,j],y[::200,j], zs=z[::200,j])
        
    set_axes_equal(ax)
    """
    from matplotlib import cm
    plt.figure()
    plt.plot(x[::200,j],y[::200,j])
    plt.scatter(x[::200,j],y[::200,j], c=cm.hot(np.arange(x[::200,j].size)), edgecolor='none')
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    j = np.random.choice(x.shape[1])
    ax.plot(x[::200,j],y[::200,j], zs=z[::200,j])"""
    npts = 200
    x_surface = np.concatenate([np.linspace(-R,R,npts)])
    y_surface = np.concatenate([np.linspace(-R,R,npts)])

    x_surface,y_surface = np.meshgrid(x_surface,y_surface)
    """r2 = x_surface**2+y_surface**2
    x_surface = x_surface[r2<R**2]
    y_surface = y_surface[r2<R**2]"""
    z_surface = np.sqrt(R**2-(x_surface**2+y_surface**2))
    #z_surface[]*=-1
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x_surface, y_surface, z_surface,
                       linewidth=0, antialiased=False)
    set_axes_equal(ax)

# ---------- making image ------------------

npix_img = 16
z_cutoff = 10*dz_tirf
stack = np.zeros((nsteps,npix_img*2+1, npix_img*2+1))

def coord2counts(x,y,z):
    zr = np.sqrt(2*np.log(2))*sigmaz
    omegaz = sigma_psf*np.sqrt(1+z/zr**2)
    frame = np.zeros((npix_img*2+1, npix_img*2+1))
    frame[x,y] = np.exp(-z/dz_tirf)
    
    frame = gaussian(frame,sigma = omegaz)*(sigma_psf/omegaz)**2
    return np.random.poisson(frame* brightness*dt )

for j in range(nsteps):
    if j%500==0:
        print("Processing frame {}".format(j))
        

    
    x1, y1, z1 = x[j], y[j], z[j]
    x1, y1, z1 = sphere2ellipsoid(np.array([x1,y1,z1]), motion_z[j])
    # print("z1",z1.min(),R*motion_z[j])
    z1+=R*motion_z[j]
    # round is necessary to ensure fair distribution of parts and not concentration in the centre
    positions_new = np.array([x1,y1]).T/psize + npix_img
    positions_new = np.round(positions_new).astype(int) 
    
    msk0 = np.logical_and(positions_new>=0,positions_new<npix_img*2+1).all(axis=1)
    msk3 = np.logical_and(msk0,z1<z_cutoff)
    positions_new = positions_new[msk3,:]
    znew = z1[msk3]
    for k in range(positions_new.shape[0]):
        frame = coord2counts(positions_new[k,0], positions_new[k,1],znew[k])
        stack[j]+=frame
    
    """for k in range(len(z1)):
        stack[j, positions_new[k, 0],positions_new[k, 1]]+=np.exp(-z1[k]/dz_tirf)
    stack[j] = gaussian(stack[j],sigma = sigma_psf)"""
    
"""stack = stack * brightness*dt 
stack = np.random.poisson(stack)"""

plt.figure()
plt.subplot(221)
plt.imshow(stack[0])
plt.subplot(222)
plt.imshow(stack[2000])
plt.subplot(223)
plt.imshow(stack[9999])
plt.subplot(224)
plt.imshow(stack.sum(axis=0))

if save:
    savepath = "/home/aurelien/Data/imFCS simulations/GUV/"
    tifffile.imsave(savepath+"simulationGUV_motion4_zpsf_1msexp_D{:.2f}.tif".format(D),stack)

def chord2arc(d,R):
    """d is chord size, R is radius"""
    return 2*R*np.arcsin(d/(2*R))

def get_chord(x,y,z,n0,n1):
    return np.sqrt( (x[n0]-x[n1])**2+(y[n0]-y[n1])**2+(z[n0]-z[n1])**2 )

def get_msds(ts,x,y,z):
    msds = list()
    for t in ts:
        chord = get_chord(x,y,z,0,t)
        msds.append( chord2arc(chord,R)**2 )
    return np.asarray(msds)

all_mms = list()
ts = np.arange(2,500,10)

plt.figure()
plt.subplot(121)
for j in range(nparts):
    
    mm = get_msds(ts,x[:,j],y[:,j],z[:,j])
    all_mms.append(mm)
    plt.plot(ts*dt,mm)
all_mms = np.asarray(all_mms)
mean_msd = all_mms.mean(axis=0)
plt.plot(ts*dt,mean_msd,"k--")

plt.subplot(122)
from scipy.stats import linregress

lr = linregress(ts*dt, mean_msd)

plt.plot(ts*dt,mean_msd,label='Mean msd')
plt.plot(ts*dt,lr.slope*ts*dt+lr.intercept,"k--", label="Fit")
plt.xlabel('Time (s)')
plt.ylabel('MSD (um2/s)')
plt.legend()

print('MSD: {:.2f} um/s'.format(lr.slope))

plt.figure()
plt.plot(stack.mean(axis=(1,2)))
plt.title("Stack intensntiy variation with time")
plt.xlabel('time')
plt.ylabel("Intensity variation")