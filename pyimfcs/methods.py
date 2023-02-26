#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:46:51 2023

@author: aurelienb
"""
import numpy as np
from pyimfcs.metrics import intensity_threshold

def downsample_mask(im,nsum,hard_th=True):
    """Downsamples a boolean mask"""
    u,v = im.shape
    out = np.zeros((u//nsum,v//nsum))
    for j in range(u//nsum):
        for k in range(v//nsum):
            px = im[j* nsum:(j+1)*nsum, k * nsum:(k+1)*nsum ]
            if len(np.unique(px[px>0]))>1: # excludes pixels with 2 different masks in them
                px = 0
            elif hard_th:
                px = min(np.unique(px))
            else:
                px = max(np.unique(px))
            out[j,k] = px
    return out
            

def indices_intensity_threshold(mask,ith,intensities):
    """Takes a mask (label) image and applies the intensity thresholding function
    to each individual element of the labeled mask. Considers the maximum within
    the label instead of the global maximum to limit skewing caused by brighter
    elements in the FOV"""
    assert( (np.array(intensities.shape)==np.array(mask.shape)).all() )
    indices=np.zeros_like(mask)
    for maskval in np.unique(mask[mask>0]):
        # intensity threshold for each mask respectively
        ithr2 = ith*(
            np.percentile(intensities[mask==maskval],98)-np.percentile(intensities,2)
            )+np.percentile(intensities,2)
        tmp_mask2 = np.logical_and(mask==maskval,intensities>ithr2)
        indices[tmp_mask2]=maskval
    return indices