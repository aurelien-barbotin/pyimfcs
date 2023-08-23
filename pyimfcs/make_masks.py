#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 16:16:41 2023

@author: aurelienb
"""

import numpy as np
import matplotlib.pyplot as plt
import tifffile
from skimage.filters import try_all_threshold, threshold_otsu

from skimage.feature import peak_local_max
from skimage.segmentation import watershed
import scipy.ndimage as ndi

def make_masks(path, plot=False, psize = None, celldiam = None, save = True):
    if path[-9:]=='_mask.tif':
        print(path,"Is a mask and not processed")
        return
    if celldiam is None or psize is None:
        raise ValueError('Please specify pixel size and cell diameter')
        
    stack = tifffile.imread(path)
    
    stackmean = stack.mean(axis=0)
    
    """plt.subplot(122)
    plt.hist(stackmean.reshape(-1),bins=50)"""
    
    threshold_value = threshold_otsu(stackmean)
    
    coordinates = peak_local_max(stackmean, min_distance=int(celldiam/psize/2),
                                 exclude_border=False,threshold_abs=threshold_value)
    

    markers = np.zeros(stackmean.shape, dtype=bool)
    markers[tuple(coordinates.T)] = True
    markers, _ = ndi.label(markers)
    
    masks = watershed(-stackmean,markers=markers,mask=stackmean>threshold_value)
    
    if plot:
        plt.figure()
        plt.subplot(211)
        plt.imshow(stackmean)
        for j in range(coordinates.shape[0]):
            x,y = coordinates[j]
            plt.plot(y,x,"rx")
        plt.subplot(212)
        plt.imshow(masks,cmap="tab20")
    if save:
        tifffile.imwrite(path.rstrip('.tif')+'_mask.tif',masks)
    return masks
path="/run/user/1000/gvfs/smb-share:server=data.micalis.com,share=proced/microscopy/ZEISS/Aurelien/2023_07_28/2_RNR-stat_1hinPBS/Image 95.tif"
make_masks(path,plot=True,psize=0.16,celldiam=1)

if __name__=="__main__":
    # Script to call in command line
    import argparse
    import glob
    import os
    
    parser = argparse.ArgumentParser(description=
         'Generate masks of individual cells in TIRF using a watershed algorithm')
    
    parser.add_argument("-p", "--pixel_size", help="Pixel size")
    parser.add_argument("-c", "--cell_diameter", help="Cell diameter")
    parser.add_argument("-n", "--name", help="Makes mask of a single file")
    args = parser.parse_args()
    
    if args.name is not None:
        files = [args.name]
    else:
        files = glob.glob('*.tif')
        
    for file in files:
        name = os.path.split(file)[-1]
        pixel_size=0.16
        if args.pixel_size is not None:
            pixel_size=float(args.pixel_size)
        cell_diameter = args.cell_diameter
        if cell_diameter is None:
            raise KeyError('Please specify a cell diameter')
        cell_diameter = float(cell_diameter)
        make_masks(file,psize=pixel_size,celldiam=cell_diameter)
