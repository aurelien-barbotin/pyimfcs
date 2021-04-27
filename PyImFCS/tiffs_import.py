# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 15:47:29 2021

@author: abarbotin
"""

import glob
import numpy as np
import tifffile
import os

def gnumber(fn):
    return int(os.path.split(fn)[-1].split("_t")[-1].rstrip(".TIF"))

files = glob.glob(r"C:\Users\abarbotin\Desktop\analysis_tmp\2020_04_07\test3_40percent/imFCS100percent_001_t[0-9]*")
numbers = [gnumber(w) for w in files]
files.sort(key = gnumber)
export = True

if export:
    frames = np.stack([tifffile.imread(w) for w in files])
    new_path = os.path.split(files[0])[0]+"_stack.tif"
    tifffile.imwrite(new_path,frames)

f = open(r"O:/microscopy/NIKON/Aurelien/2021_04_01_imFCS/100percentSucrose/imFCS20percent_002.nd","r")
out = ""
for line in f.readlines():
    out+=line
    
f.close()

from PIL import Image
from PIL.TiffTags import TAGS

with Image.open(files[0]) as img:
    meta_dict = {TAGS[key] : img.tag[key] for key in img.tag.keys()}
    
description = meta_dict['ImageDescription']
