# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 12:35:33 2021

@author: abarbotin
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
path="C:/Users/abarbotin/Desktop/Charlene data/C1-Position6.tif kept stack kept stack_seg.npy"
uu=np.load(path,allow_pickle=True).item()

masks = uu["masks"]
img = uu["img"]
masks=masks.astype(float)
masks[masks==0] = np.nan

plt.figure()
plt.imshow(img,cmap='gray')
plt.imshow(masks%10,cmap="tab10",alpha=0.3)
plt.imshow(masks,cmap="tab10",alpha=0.003)