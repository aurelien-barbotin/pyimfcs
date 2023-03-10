#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 18:17:27 2022

@author: aurelien

"""
import numpy as np
import pandas as pd

from pyimfcs.class_imFCS import StackFCS

def summarise_df(df):
    """Merges results of a Dataframe containing multiple experiments"""
    means = df.groupby('repeat')['D [µm²/s]'].mean()
    std = df.groupby('repeat')['D [µm²/s]'].std()
    medians = df.groupby('repeat')['D [µm²/s]'].median()
    count = df.groupby('repeat')['D [µm²/s]'].count()
    repeat = df.groupby('repeat')['repeat'].last()
    filename = df.groupby('repeat')['filename'].last()
    binning = df.groupby('repeat')['binning'].last()

    out_df = pd.DataFrame.from_dict({
                        "filename": filename,
                        "binning": binning,
                        "repeat":repeat,
                        "Mean":means,
                        "Stdev": std,
                        "Median":medians,
                        "Count":count})
    return out_df
    
def merge_fcs_results(out_name, files, ith = None, 
                      chi_threshold = None,use_mask=False):
    """Wrapper function to merge all experiment results in a single excel file"""
    
    # Step1: browse through every file and extract what we want
    if len(files)==0:
        raise ValueError('No files selected')
    nsums = None
    all_dfs={}
    for nfile,file in enumerate(files):
        print(file)
        stack = StackFCS(file,load_stack=False)
        stack.load()
        stack_res = stack.extract_results(ith=ith,
                                        chi_threshold=chi_threshold,use_mask=use_mask)
        diffs = stack_res["diffusion_coefficients"]
        chis = stack_res["non_linear_chis"]
        nmolecules = stack_res["number_molecules"]
        indices = stack_res['indices']
        # check that nsums are identical in all files
        nsums_tmp = sorted(diffs.keys())
        if nsums is None:
            nsums = nsums_tmp
            for ns in nsums:
                all_dfs[ns]=[]
        else:
            for j in range(max(len(nsums),len(nsums_tmp))):
                if nsums[j]!=nsums_tmp[j]:
                    raise KeyError('All files were not processed identically')
                    
        # populates dataframes
        for nsum in nsums:
            fname = file
            diff = diffs[nsum]
            chi = chis[nsum]
            nmol = nmolecules[nsum]
            indice = indices[nsum].astype(int)
            repeats_arr = np.full(diff.size, nfile)
            name_arr = np.full(diff.size, fname)
            nsum_arr = np.full(diff.size, nsum)
            
            out_arr = np.array([name_arr,repeats_arr, diff, nsum_arr, chi, nmol,indice]).T
            
            df = pd.DataFrame(out_arr, columns = 
                              ["filename", "repeat","D [µm²/s]","binning",
                               "fit error", "N","label"])
            df = df.astype({'filename':"str",
                           "repeat":"int",
                           "D [µm²/s]":"float",
                           "fit error":"float",
                           "N": "float",
                           "label":"int"})
            all_dfs[nsum].append(df)
        
    parameters_dict = {"chi_threshold": chi_threshold,
                       "intensity_threshold":ith}
    # Extracts global parameters
    for nsum in nsums:
        all_dfs[nsum] = pd.concat(all_dfs[nsum])
        knm = "nsum {}".format(nsum)
        parameters_dict[knm+"_median"] = np.median(all_dfs[nsum]["D [µm²/s]"].values)
        parameters_dict[knm+"_mean"] = np.mean(all_dfs[nsum]["D [µm²/s]"].values)
        parameters_dict[knm+"_std"] = np.std(all_dfs[nsum]["D [µm²/s]"].values)
    
    # Writes in excel file
    with pd.ExcelWriter(out_name+".xlsx") as writer:  
        dfs_total = []
        for nsum in nsums:
            dfpooled = all_dfs[nsum]
            dfs_total.append(summarise_df(all_dfs[nsum]))
            dfpooled.to_excel(writer, sheet_name = "nsum {}".format(nsum))
        dfs_total = pd.concat(dfs_total)
        dfs_total.to_excel(writer, sheet_name = "summaries all")
        df_pars = pd.DataFrame(parameters_dict, index=[0]).T
        df_pars.to_excel(writer, sheet_name = "parameters")