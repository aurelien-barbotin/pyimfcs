#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 12:53:03 2021

@author: aurelien
"""

import numpy as np
import matplotlib.pyplot as plt

from PyImFCS.fitting import make_Gim2D,make_Gim3D
plt.close('all')

gim1 = make_Gim3D(5,0.11,0.48,ginf=False)

x = np.logspace(-3,5)

ncurves=10
Ns = np.ones(ncurves)

dmin = 0.05
dmax = 50

Ds = np.linspace(dmin,dmax,ncurves)

curves = [gim1(x,Ns[j],Ds[j]) for j in range(ncurves)]

plt.figure()
for j in range(ncurves):
    plt.semilogx(x, curves[j])
plt.xlabel(r'$\tau\ (s)$')
plt.ylabel(r'$G(\tau)$')