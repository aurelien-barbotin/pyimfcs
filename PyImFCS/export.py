#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 18:17:27 2022

@author: aurelien

Downstream methods for easy results export
"""
import h5py

import numpy as np
import pandas as pd

from PyImFCS.io import get_dict, save_as_excel
from PyImFCS.class_imFCS import new_chi_square

def get_fit_error(files,nsums = None, intensity_threshold = None, chi_threshold = None):
    """Retrieves and calculates diffusion coefficients and non-linear fit errors
    in a series of *.h5 FCS files.
    
    Parameters:
        files (list): filenames (str) of h5 files
        nsums (list): Optional. if specified, list of binning values to export. If not, finds
            the available binning values in the dataset, returns an error if not
            all datasets were processed similarly
        intensity_threshold (float): optional. If specified, intensity threshold
            in fraction of max intensity. Has to be between 0 and 1
        chi_threshold (float): optional. If specified, removes all curves below 
            a certain quality thresold
    Returns:
        list: [diffusion coefficients, error_metric]. Each item is a list of dictionaries."""
    all_curves = []
    all_fits = []
    all_chis = []
    all_chis_new = []
    all_diffs = []
    all_labels = []
    
    diffs_out = list()
    chis_out = list()
    
    # if nsums is not specified; take the nsums in every file and check that 
    # it is consistent across the dataset
    check_nsums = False
    nsums_backup = None
    if nsums is None:
        check_nsums = True
    
    for k,file in enumerate(files):
        h5f = h5py.File(file,mode='r')
        dict_diffcoeff = get_dict(h5f, "parameters_fits")
        dict_curves = get_dict(h5f, "correlations")
        dict_curves_fits = get_dict(h5f, "yhat")
        dict_traces = get_dict(h5f,'traces')
        h5f.close()
        
        if check_nsums:
            nsums = list(dict_curves.keys())
            if nsums_backup is None:
                nsums_backup = nsums
            if sorted(nsums)!=sorted(nsums_backup):
                raise KeyError('Not all files were processed with the same parameters')
            
        diffs_out.append( dict(zip(nsums,[[] for w in nsums])))
        chis_out.append(dict(zip(nsums,[[] for w in nsums])))
        
        for jj, nsum in enumerate(nsums):
            diffcoeffs = dict_diffcoeff[nsum][:,:,1]
            curves = dict_curves[nsum]
            curves_fits = dict_curves_fits[nsum]
            traces = dict_traces[nsum]
            
            intensities = traces.mean(axis=2).reshape(-1)
            # xcurves = curves[:,:,:,0]
            
            ycurves = curves[:,:,:,1]
            
            chis = np.sqrt(((ycurves-curves_fits)**2).sum(axis=2) )
            
            chis_new = np.zeros_like(chis)
            
            for i in range(chis_new.shape[0]):
                for j in range(chis_new.shape[1]):
                    
                    chis_new[i,j] = new_chi_square(ycurves[i,j], curves_fits[i,j])
        
            curves_reshaped = curves.reshape((curves.shape[0]*curves.shape[1],curves.shape[2],2))
            fits_reshaped = curves_fits.reshape((curves.shape[0]*curves.shape[1],curves.shape[2]))
            
            msk = np.ones_like(intensities, dtype = bool)
            msk[diffcoeffs.reshape(-1)<0] = 0
            
            if intensity_threshold is not None:
                msk = np.logical_and(msk,
                                     intensities>(intensities.max()*intensity_threshold))
            if chi_threshold is not None:
                msk = np.logical_and(msk, chis_new.reshape(-1)<chi_threshold)
            
            chis_new = chis_new.reshape(-1)[msk]
            curves_reshaped = curves_reshaped[msk]
            fits_reshaped = fits_reshaped[msk]
            diffs = diffcoeffs.reshape(-1)[msk]
            all_curves.append(curves_reshaped)
            all_fits.append(fits_reshaped)
            all_chis.append(chis.reshape(-1)[msk])
            all_chis_new.append(chis_new)
            all_diffs.append(diffs)
            # all_intensities.append()
            all_labels.append(file.split('/')[-1]+"_{}".format(nsum))
            
            diffs_out[k][nsum] = diffs
            chis_out[k][nsum]= chis_new
    return diffs_out, chis_out

def merge_fcs_results(files, out_name, intensity_threshold = None, chi_threshold = None):
    """Wrapper function to merge all experiment results in a single excel file"""
    
    if len(files)==0:
        raise ValueError('No files selected')
        
    all_diffs,all_chis = get_fit_error(files, nsums = None, intensity_threshold=intensity_threshold, 
                                       chi_threshold=chi_threshold)
    assert(len(all_diffs)>0)
    nsums = sorted(list(all_diffs[0].keys()))
    for w in all_diffs:
        assert(sorted(list(w.keys())) ==nsums for w in all_diffs)
        
    parameters_dict = {"chi_threshold": chi_threshold,
                       "intensity_threshold":intensity_threshold}
    
    save_as_excel(out_name,files,nsums,all_diffs,all_chis, 
                  parameters_dict=parameters_dict)
    
def merge_excels(fileslist_list, out_name, keep_single_indices = False, conditions = None):
    """Merge ther esults of FCS experiments in a single file.
    Parameters:
        fileslist_list (list): list of lists. Each item of the main list contains 
            the excel files corresponding to one condition
        out_name (str): save name
        keep_single_indices (bool): if true, gives an indes to single acquisition. 
            If False, single index to an excel file.
        conditions (list): if specified, gives specific condition names 
            (e.g control or experiment) to acquisitions"""
    
    if len(fileslist_list)>1 and conditions is None:
        return ValueError('Please specify condition names')
    assert( conditions is None or len(conditions)==len(fileslist_list))
    
    for inumber, files in enumerate(fileslist_list):
        excels = [pd.ExcelFile(w) for w in files]
        
        # xl0 = excels[0]
        
        all_names = [w.sheet_names for w in excels]
        names0 = all_names[0]
        kept_names = list()
        for name in names0:
              if np.all([name in sublist for sublist in all_names]):
                  kept_names.append(name)
        
        all_dfs = {}
        for name in kept_names:
            dfs = []
            maxindex = 0
            for j, xl in enumerate(excels):
                df = xl.parse(sheet_name=name)
                if name=="parameters":
                    fname = files[j]
                    df["file"] = fname
                    
                else:
                    if keep_single_indices:
                        df['repeat']+=maxindex
                        maxindex = df['repeat'].values.max()+1
                        # print("maxindex", maxindex)
                    else:
                        df['repeat'] = maxindex
                        maxindex += 1
                if conditions is not None:
                    df["condition"] = conditions[inumber]
                dfs.append(df)
                    
            dfs = pd.concat(dfs)
            all_dfs[name] = dfs
            
    with pd.ExcelWriter(out_name) as writer:  
            for name in kept_names:
                df_pars = all_dfs[name]
                df_pars.to_excel(writer, sheet_name = name)