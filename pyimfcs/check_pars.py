#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 17:00:52 2022

@author: aurelien
"""

import numpy as np

from pyimfcs.class_imFCS import get_image_metadata
from termcolor import colored

path = "/home/aurelien/Data/2022_01_06/4_NR12A/GP/Image 53_GP.tif"

gain_name = "Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|Detector|AmplifierGain"
exptime_name = "Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|CameraIntegrationTime"
tubelens_name = "Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|TubeLensPosition"
acqmode_name = "Information|Image|Channel|AcquisitionMode"
intensity_name = "Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|Attenuator|Transmission"
timespan_name = "Experiment|AcquisitionBlock|TimeSeriesSetup|Switch|SwitchAction|SetIntervalAction|Interval|TimeSpan|Value"
xresolution_name = 'Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX'
yresolution_name = 'Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY'
creationdate_name = 'Information|Document|CreationDate'
names_to_test = {"gain": gain_name,
                 "exposure time": exptime_name,
                 "tube lens": tubelens_name,
                 "acquisition mode":acqmode_name,
                 "laser intensity":intensity_name,
                 'xresolution': xresolution_name,
                 'yresolution':yresolution_name,
                 'creation date':creationdate_name
                 }
all_par_names = {"gain": gain_name,
                 "exposure time": exptime_name,
                 "tube lens": tubelens_name,
                 "acquisition mode":acqmode_name,
                 "laser intensity":intensity_name,
                 'time span': timespan_name,
                 'xresolution': xresolution_name,
                 'yresolution':yresolution_name,
                 'creation date':creationdate_name
                 }

def get_parameter(easy_name, metadata):
    """Given a human readable name, key to names_to_test, finds the corresponding
    parameter in a metadata dict"""
    name = names_to_test[easy_name] 
    if name in metadata:
        ref = metadata[name]
    else:
        candidates_list = list(filter(lambda x: name in x,metadata.keys()))
        if len(candidates_list)==1:
            ref = metadata[candidates_list[0]]
        else:
            ref = 'Not found'
    return ref

short = lambda x: x.split("/")[-1].split(".")[0]
def compare_channel_file(path, verbose = True):
    all_same = True

    metadata = get_image_metadata(path)
    nchannels = len(list(filter(lambda x: gain_name in x,metadata.keys())))
    outstring = "# channels: {}\n".format(nchannels)
    filename = path.split("/")[-1].split(".")[0]
    par_dict ={}
    for nn in names_to_test.keys():
        ref = metadata[names_to_test[nn]+" #1"]
        par_dict[nn] = ref
        for j in range(nchannels-1):
            to_test = metadata[names_to_test[nn]+" #{}".format(j+2)]
            if to_test != ref:
                all_same = False
                outstring+= colored("{}: have different values ".format(nn),"red")
            outstring+= "{} value: {}\n".format(nn, ref)
    if all_same:
        header = "{}: all channels have identical parameters".format(filename)
    else:
        header = "{}: channels have different parameters".format(filename)
    header = header + "\n" + "-"*len(header)+"\n"
    outstring = header+outstring
    if verbose:
        print(outstring)
    return all_same, par_dict

def get_singlechannel_file(path, verbose = True):
    
    metadata = get_image_metadata(path)
    par_dict ={}
    for nn in names_to_test.keys():
        name = names_to_test[nn]
        if name in metadata:
            ref = metadata[name]
        else:
            candidates_list = list(filter(lambda x: name in x,metadata.keys()))
            if len(candidates_list)==1:
                ref = metadata[candidates_list[0]]
            else:
                ref = 'Not found'
            
        par_dict[nn] = ref
    return par_dict

def compare_folder(files_list, verbose = False):
    # tested control positive and neg, works OK
    assert len(files_list)>1
    all_dicts =  []
    all_b1s = []
    for file in files_list:
        b1, dic = compare_channel_file(file, verbose = verbose)
        all_b1s.append(b1)
        all_dicts.append(dic)
        
    if not all(all_b1s):
        print("Some files have different parameters in different channels ")
    
    dic0 = all_dicts[0]
    keynames = dic0.keys()
    same_params_perfile = True
    
    refname=files_list[0]
    differents_string = ""
    
    for kn in keynames:
        for j,dic in enumerate(all_dicts):
            if dic0[kn] != dic[kn]:
                same_params_perfile = False
                differents_string+="{} and {} have different {} parameters\n".format(
                    short(refname),short(files_list[j]), kn)
    if same_params_perfile:
        outstring = "All files have the same acquisition parameters:\n"
        outstring+="-"*len(outstring)+"\n"
        for namep, parval in dic0.items():
            outstring+="{}: {}\n".format(namep, parval)
    else:
        outstring = differents_string
    print(outstring)
    return same_params_perfile

def compare_dicts(all_dicts):
    dic0 = all_dicts[0]
    keynames = dic0.keys()
    same_params_perfile = True
    
    differents_string = ""
    
    for kn in keynames:
        for j,dic in enumerate(all_dicts):
            if dic0[kn] != dic[kn]:
                same_params_perfile = False
    if same_params_perfile:
        outstring = "All files have the same acquisition parameters:\n"
        outstring+="-"*len(outstring)+"\n"
        for namep, parval in dic0.items():
            outstring+="{}: {}\n".format(namep, parval)
    else:
        outstring = differents_string
    print(outstring)
    return same_params_perfile

def compare_single_experiments(files_list):
    all_dicts = [get_image_metadata(w) for w in files_list]
    dic0 = all_dicts[0]
    keynames = names_to_test.values()
    same_params_perfile = True
    refname=files_list[0]
    
    differents_string = ""
    
    for kn in keynames:
        for j,dic in enumerate(all_dicts):
            if dic0[kn] != dic[kn]:
                same_params_perfile = False
                differents_string+="{} and {} have different {} parameters\n".format(
                    short(refname),short(files_list[j]), kn)
    if same_params_perfile:
        outstring = "All files have the same acquisition parameters:\n"
        outstring+="-"*len(outstring)+"\n"
        for namep, parval in dic0.items():
            outstring+="{}: {}\n".format(namep, parval)
    else:
        outstring = differents_string
    print(outstring)
    return same_params_perfile

def get_single_parameter(path, parname):
    """Returns the value of a parameter (simplified as in names_to_test) in a
    given tiff file
    
    Parameters:
        path (str): path to a tifffile from the Zeiss
        parname (str): name of the parameter of interest. Has to be in names_to_test
    Returns:
        ndarray: array of parameter values (one per channel)
    
    Tested on one parameter, one single_channel image and one 2-channel image"""
    metadata = get_image_metadata(path)
    
    key = all_par_names[parname]
    if key in metadata:
        return np.array([metadata[key]])
    
    elif key+' #1' in metadata:
        out = []
        nkey = key+' #1'
        jind = 1
        while nkey in metadata:
            out.append(metadata[nkey])
            jind+=1
            nkey = key+' #{}'.format(jind)
        return np.array(out)
    raise KeyError('Parameter not in file')

if __name__=="__main__":
    # Script to call in command line to check named metadata
    import json
    import argparse
    import glob
    import os
    
    parser = argparse.ArgumentParser(description='Find stats of a tif file generated with Zen. Calls the script pyimfcs.check_pars')
    parser.add_argument("-a", "--all", help="Finds stats on all tif files in folder")
    parser.add_argument("-n", "--name", help="Finds stats of a single file")
    parser.add_argument("-p", "--parameter", help="If you want to plot only one of the parameters")
    args = parser.parse_args()
    
    if args.all is not None:
        files = glob.glob(args.all+'/*.tif')
    elif args.name is not None:
        files = [args.name]
    else:
        raise ValueError('Please use the parameters -a or -n')
    for file in files:
        name = os.path.split(file)[-1]
        try:
            if args.parameter is None:
                print(name,'\n------')
                print(json.dumps(get_singlechannel_file(file),indent=2))
            else:
                print(name,end=": ")
                print(get_singlechannel_file(file)[args.parameter])
        except Exception as e:
            print("Metadata unreadable")
            print(e)
        print()
