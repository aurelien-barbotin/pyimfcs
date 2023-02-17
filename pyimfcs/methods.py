#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:46:51 2023

@author: aurelienb
"""
import numpy as np

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
            
