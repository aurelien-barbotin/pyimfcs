#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 14:43:26 2022

@author: aurelienb
"""

import numpy as np
import matplotlib.pyplot as plt

from pyimfcs.class_imFCS import StackFCS
from pyimfcs.fitting import Fitter

path = "/home/aurelienb/Data/simulations/2022-11-08 14:28:58/stack.tif"
stack = StackFCS(path)
stack.load(light_version=True)


yscale = 0.16
sigmaxy = 0.2
parameters_dict = {"a":yscale, "sigma":sigmaxy}
ft = Fitter("2D",parameters_dict, ginf=True)

# ft.p0 = (1/6,3,0)
x,y = 8,8
nsum = 2
curves = stack.correl_dicts[nsum]
curve = curves[x,y]

popt,yh = ft.fit(curve)
plt.figure()
plt.semilogx(curve[:,0],curve[:,1])
plt.semilogx(curve[:,0],yh, color='k', linestyle='--')
