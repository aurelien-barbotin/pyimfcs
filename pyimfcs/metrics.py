#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:52:01 2022

@author: aurelien
"""

import numpy as np
from scipy.ndimage import label

def new_chi_square(y,yh):
    if yh[0]==0:
        return np.inf
    diff = (y-yh)/yh[0]
    diffpos = diff>0
    diffneg = diff<0
    
    diff =  diff**2
    chi = 0
    lab, num_features = label(diffpos)
    for j in range(1,num_features+1):
        msk = lab==j
        chi+=diff[msk].sum()*(np.count_nonzero(msk)-1)
        
    lab, num_features = label(diffneg)
    for j in range(1,num_features+1):
        msk = lab==j
        chi+=diff[msk].sum()*(np.count_nonzero(msk)-1)
    
    return chi/yh.size

def intensity_threshold(th, intensities):
    """Claculates an intensity threshold in presence of baseline background"""
    ithr = th*(
        np.percentile(intensities,98)-np.percentile(intensities,2)
        )+np.percentile(intensities,2)
    return ithr