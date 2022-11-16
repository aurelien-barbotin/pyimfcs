#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 14:17:00 2022

@author: aurelienb
"""
import glob
from pyimfcs.io import get_image_metadata, get_h5_metadata, extract_from_h5
import numpy as np
import os

def time_str_2_min(creation_date):
    creation_time = creation_date[-8:].split(':')
    crt_min = int(creation_time[0])*60+int(creation_time[1])+float(creation_time[2])/60
    return crt_min

def extract_timeseries_fromfiles(files, out_time, acqtime_min = 4.5, nsum=3,
                       time_threshold = 700, intensity_threshold = 0.8, chthr=0.03):
    """Same as the others but takes as argument a list of files and not directly a path"""
    out_time_min = time_str_2_min(out_time)
    
    all_diffs = []
    all_times = []
    for file in files:
        diffs = extract_from_h5(file,nsum=nsum, intensity_threshold=intensity_threshold, 
                                chi_threshold = chthr)
        if os.path.isfile(file[:-3]+".tif"):
            metadata = get_image_metadata(file[:-3]+".tif")
        else:
            metadata = get_h5_metadata(file)
        datekey = list(filter(lambda x: "Date" in x, metadata.keys()))
        assert len(datekey)==1
        creation_date = metadata[datekey[0]]
        crt = time_str_2_min(creation_date)
        all_times.append(crt-out_time_min+acqtime_min/2)
        all_diffs.append(diffs)
        
    all_times = np.array(all_times)
    difval = np.array([np.median(w) for w in all_diffs])
    diferr = np.array([np.std(w) for w in all_diffs])
    timerr = np.ones_like(difval)*acqtime_min/2
    
    msk = all_times<time_threshold
    difval = [x for _, x in sorted(zip(all_times[msk], difval[msk]))]
    diferr = [x for _, x in sorted(zip(all_times[msk], diferr[msk]))]
    all_times = sorted(all_times[msk])
    return all_times, difval, diferr, timerr[msk]


def extract_timeseries_split(path, out_time, acqtime_min = 4.5, nsum=3,
                       time_threshold = 700):
    """To extract time series in presence of reprocessed stacks, if for each 
    dataset the full stack was processed and saved as h5but also the first and 
    second half"""
    out_time_min = time_str_2_min(out_time)
    files = glob.glob(path+"*.h5")
    print(files)
    
    all_diffs = []
    all_times = []
    all_timerr= []
    for file in files:
        diffs = extract_from_h5(file,nsum=nsum, intensity_threshold=0.5, 
                                chi_threshold = 0.03)
        # if "_firsthalf" in file:
        tifname = "".join(file.split("_firsthalf"))
        tifname = "".join(tifname.split("_lasthalf"))
        metadata = get_image_metadata(tifname[:-3]+".tif")
        datekey = list(filter(lambda x: "Date" in x, metadata.keys()))
        assert len(datekey)==1
        creation_date = metadata[datekey[0]]
        crt = time_str_2_min(creation_date)
        all_diffs.append(diffs)
        t0 = crt-out_time_min + acqtime_min/2
        if "_firsthalf" in file:
            time = t0-acqtime_min/4
            timerr = acqtime_min/4
        elif "lasthalf" in file:
            time = t0 + acqtime_min/4
            timerr = acqtime_min/4
        else:
            time = t0
            timerr = acqtime_min/2
        all_times.append(time)
        all_timerr.append(timerr)
        
    all_times = np.array(all_times)
    difval = np.array([np.median(w) for w in all_diffs])
    diferr = np.array([np.std(w) for w in all_diffs])
    timerr = np.array(all_timerr)
    
    msk = all_times<time_threshold
    difval = [x for _, x in sorted(zip(all_times[msk], difval[msk]))]
    diferr = [x for _, x in sorted(zip(all_times[msk], diferr[msk]))]
    timerr = [x for _, x in sorted(zip(all_times[msk], timerr[msk]))]
    all_times = sorted(all_times[msk])
    return np.array(all_times), np.array(difval), np.array(diferr), np.array(timerr)
