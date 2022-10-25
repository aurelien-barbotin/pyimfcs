#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 16:06:15 2021

@author: aurelien
"""

import matplotlib.pyplot as plt

from pyimfcs.class_imFCS import StackFCS
from pyimfcs.blcorr import blexp_double_offset

from pyimfcs.fitting import Fitter
from pyimfcs.plotting import multiplot_stack
from pyimfcs.io import get_image_metadata

import os

def get_metadata_zeiss(file):
    metadata = get_image_metadata(file)
    dt = metadata['finterval']
    if 'Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1' in metadata:
        xscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'] * 10 ** 6
        yscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1'] * 10 ** 6
    else:
        xscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX'] * 10 ** 6
        yscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY'] * 10 ** 6
    
    return dt, xscale, yscale

def batch_bacteria_process(files,first_n = 3000, last_n = 0, nsums=[2,3], nreg=4000,
                           plot=False, default_dt= None, default_psize = None):
    """Processes a series of FCS exeriments saved as tiffs.
    Parameters:
        files (list): list of paths to tiff images
        first_n (int): first n frames to remove from every acquisition
        last_n (int): last n frames to remove from every acquisition
        nsums (list): of int, binning values to use for FCS
        nreg (int): averages frames n by n to determine the drift correction. 
            No drift correction if 0
        plot (bool): if True, plots intermediate results
        default_dt (float): if specified, default value to use for frame interval. 
            Units: s. Used if value could not be retrieved from metadata
        default_psize (float): if specified, default value to use for pixel size,
            in micrometers. Used only if could not be retrieved from metadata"""
    for path in files:
        print('Processing',path)
        stack = StackFCS(path, background_correction = True,                     
                             first_n = first_n, last_n = last_n, clipval = 0)
        if not stack.metadata_fully_loaded:
            print('Metadata was not loaded in file {}'.format(path.split(os.sep)[-1]))
            if default_dt is not None:
                stack.dt = default_dt
            if default_psize is not None:
                stack.xscale = default_psize
                stack.yscale = default_psize
        xscale = stack.xscale
        yscale = stack.yscale
        assert(xscale==yscale)
        stack.registration(nreg,plot=plot)
            
        stack.set_bleaching_function(blexp_double_offset)
        # stack.set_bleaching_function(bleaching_correct_sliding,wsize = 5000)
        
        curves_avg = stack.binned_average_curves(nsums,n_norm=2, plot=False)
        # !!! TO CHANGE
        sigmaxy = 0.2
        parameters_dict = {"a":yscale, "sigma":sigmaxy}
        ft = Fitter("2D",parameters_dict, ginf=True)
        
        stack.fit_curves(ft,xmax=None)
        
        chi_threshold = 0.03
        if plot:
            multiplot_stack(stack,nsums[-1], chi_threshold = chi_threshold)
            ttl = path.split(os.sep)[-1]
            plt.suptitle(ttl)
            
            plt.figure()
            plt.plot(stack.trace())
            plt.xlabel('Frame nr')
            plt.ylabel('Average Intensity')
            plt.suptitle(ttl)
        
        stack.save()
        print("saving stack")


