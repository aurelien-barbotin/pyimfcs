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
from pyimfcs.plotting import interactive_fcs_plot
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
                           plot=False, default_dt= None, default_psize = None, 
                           fitter = None, export_summaries = True, 
                           chi_threshold = 0.03, ith=0.8):
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
            in micrometers. Used only if could not be retrieved from metadata
        fitter (Fitter): object of the class Fitter
        export_summaries (bool): if True, saves summaries of the processing"""
    if export_summaries:
        export_path = os.path.split(files[0])[0]+"/summaries/"
        if not os.path.isdir(export_path):
            os.mkdir(export_path)
            
    for path in files:
        print('Processing',path)
        fname = "".join(path.split(os.sep)[-1].split('.')[:-1])
        export_folder = export_path+fname+"/"
        if export_summaries:
            if not os.path.isdir(export_folder):
                os.mkdir(export_folder)
        try:
            stack = StackFCS(path, background_correction = True,                     
                                 first_n = first_n, last_n = last_n, clipval = 0)
        except:
            continue
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
        if nreg>0:
            stack.registration(nreg,plot=plot or export_summaries)
            if export_summaries:
                plt.savefig(export_folder+"drift_correction.png")
                plt.close()
        stack.set_bleaching_function(blexp_double_offset)
        
        for nSum in nsums:
            stack.correlate_stack(nSum)
        if fitter is None:
            raise KeyError("Please specify a fitting method")
            """sigmaxy = 0.2
            parameters_dict = {"a":yscale, "sigma":sigmaxy}
            ft = Fitter("2D",parameters_dict, ginf=True)"""
        else:
            ft = fitter
        
        stack.fit_curves(ft,xmax=None)
        
        if plot:
            interactive_fcs_plot(stack,nsums[-1], chi_threshold = chi_threshold)
            ttl = path.split(os.sep)[-1]
            plt.suptitle(ttl)
            
            plt.figure()
            plt.plot(stack.trace())
            plt.xlabel('Frame nr')
            plt.ylabel('Average Intensity')
            plt.suptitle(fname)
        
        stack.save()
        
        if export_summaries:
            plt.figure()
            plt.subplot(121)
            for nsum in nsums:
                avgcorr = stack.average_curve(nSum=nsum,plot=False,
                                              chi_th=chi_threshold,ith=ith)
                plt.semilogx(avgcorr[:,0], avgcorr[:,1]/avgcorr[0,1], 
                     color="C{}".format(int(nsum)), label="nsum {}".format(int(nsum)))
            plt.axhline(0,color="k", linestyle='--')
            plt.xlabel(r'$\rm \tau$')
            plt.xlabel(r'$\rm G(\tau)$')
            plt.subplot(122)
            plt.plot(stack.trace())
            plt.xlabel('Frame nr')
            plt.ylabel('Average intensity')
            # add diffusion coefficients
            plt.tight_layout()
            plt.savefig(export_folder+"averages.png")
            plt.close()
        print("saving stack")


