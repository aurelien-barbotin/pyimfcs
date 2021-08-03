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
plot = True
save = False

psize = 0.16
sigma_psf = 0.2/psize
sigmaz = 4*sigma_psf
dz_tirf = 0.2 # um

dt = 1*10**-3 # s
D = 3 #um2/s

R = 8 #um

brightness = 18*10**4 #Hz/molecule

nsteps = 20000
nparts = 2000

pos0 = np.random.uniform(size = (nparts,2))*2*np.pi


moves = np.random.normal(scale = np.sqrt(2*D*dt)/R,
                         size = (nsteps,nparts,2) )

moves[0] = 0

positions = np.cumsum(moves,axis=0)
positions = positions+pos0[np.newaxis,:,:]


p1 = pos0[:,0]
phis = moves[:,:,0]/np.sin(positions[:,:,1])
positions[:,:,0] = p1[np.newaxis,:] + np.cumsum(phis,axis=0)

positions = positions%(2*np.pi)
"""mm = positions[:,:,0]>np.pi
positions[mm,0]-=np.pi
positions[mm,1]-=np.pi"""
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
z+=R

if plot:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for j in range(nparts):
        ax.plot(x[::200,j],y[::200,j], zs=z[::200,j])
        
    
    """
    from matplotlib import cm
    plt.figure()
    plt.plot(x[::200,j],y[::200,j])
    plt.scatter(x[::200,j],y[::200,j], c=cm.hot(np.arange(x[::200,j].size)), edgecolor='none')
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    j = np.random.choice(x.shape[1])
    ax.plot(x[::200,j],y[::200,j], zs=z[::200,j])"""



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
        
    msk1 = z[j]>=0
    
    x1, y1, z1 = x[j,msk1], y[j,msk1], z[j,msk1]
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
    tifffile.imsave(savepath+"simulationGUV_zpsf_1msexp_D{:.2f}.tif".format(D),stack)

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
