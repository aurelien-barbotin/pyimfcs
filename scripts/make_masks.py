#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:46:53 2023

@author: aurelienb
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import label

from tifffile import imread, imwrite
from skimage.filters import try_all_threshold, threshold_li

import glob

mainpath = "/home/aurelienb/Data/2023_05_12/1_BSLBs/"
files = glob.glob(mainpath+"*.tif")
files = list(filter(lambda x: x[-9:]!="_mask.tif",files))
for path in files:
    st = imread(path)
    img = st.mean(axis=0)
    
    th = img>threshold_li(img)
    label_img, nelts = label(th)
    
    plt.figure()
    plt.suptitle(path.split('/')[-1])
    plt.subplot(121)
    plt.imshow(img)
    plt.subplot(122)
    plt.imshow(label_img,cmap='tab20')
    imwrite(path[:-4]+"_mask.tif",label_img)

"""fig, ax = try_all_threshold(img, figsize=(10, 8), verbose=False)
plt.show()
"""