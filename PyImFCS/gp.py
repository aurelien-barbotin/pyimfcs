#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 17:27:53 2021

@author: aurelien
"""

import numpy as np
import matplotlib.pyplot as plt
import tifffile
from skimage.filters import threshold_otsu, gaussian

plt.close('all')
path = "/home/aurelien/Data/2021_06_08/Bacteria_NileRed/Image 34_2colorTIRF.tif"

stack = tifffile.imread(path).astype(float)

g = stack[0]
r = stack[1]

g = gaussian(g)
r = gaussian(r)
gp = (r-g)/(r+g)

tmap = g>threshold_otsu(g)

gp[~tmap]=np.nan


