#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 15:22:12 2021

@author: aurelien
"""


import h5py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

# parameters: [N,D,offset]
def get_dict(h5f, dname):
    """Extract a dictionary stored in a h5 file under the keyword
    dname"""
    if dname not in h5f.keys():
        print("grougrou")
        raise KeyError('Parameter not in file')
    out_dic = {}
    ds = h5f[dname]
    for key in ds.keys():
        dd = ds[key][()]
        out_dic[int(key)] = dd
    return out_dic

# dicts_to_load  = ["correlations", "parameters_fits","yhat", "traces"]

files = glob.glob("/home/aurelien/Data/2021_11_02/BSLB/POPC/*.h5")

name = files[0]
h5f = h5py.File(name, "r")

path = h5f['parameters/path'][()]

parameters_fits = get_dict(h5f, "parameters_fits")
traces = get_dict(h5f, "traces")