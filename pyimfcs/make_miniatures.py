#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:09:25 2023

@author: aurelienb
"""
import tifffile
import glob

from pyimfcs.constants import datapath
"""
files = glob.glob(datapath+"2023_08_25/1_listeria/*/*/*.tif")

for file in files:
    stack = tifffile.imread(file)
    if len(stack.shape)<3:
        continue
    avg = stack.mean(axis=0)
    tifffile.imsave(file[:-4]+"_sum.tif",avg)"""

if __name__=="__main__":
    # Script to call in command line to check named metadata
    import argparse
    
    parser = argparse.ArgumentParser(description='Makes t-axis projections on 3D FCS stacks')
    parser.add_argument("-a", "--all", help="Projects a single tif file")
    parser.add_argument("-n", "--name", help="Projects all tif files in a folder")
    args = parser.parse_args()
    
    if args.all is not None:
        files = glob.glob(args.all+'/*.tif')
    elif args.name is not None:
        files = [args.name]
    else:
        raise ValueError('Please use the parameters -a or -n')
    for file in files:
        print('Processing file ',file)
        stack = tifffile.imread(file)
        if len(stack.shape)<3:
            continue
        avg = stack.mean(axis=0)
        tifffile.imsave(file[:-4]+"_sum.tif",avg)
