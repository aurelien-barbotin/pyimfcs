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
    valid_fractions = df.groupby("repeat")["valid fraction"].mean()
    out_df = pd.DataFrame.from_dict({
                        "filename": filename,
                        "binning": binning,
                        "repeat":repeat,
                        "Mean":means,
                        "Stdev": std,
                        "Median":medians,
                        "Count":count,
                        "valid fraction": valid_fractions})
    return out_df
    
def merge_fcs_results(out_name, files, ith = None, 
                      chi_threshold = None,use_mask=False):
    """Wrapper function to merge all experiment results in a single excel file"""
    
    # Step1: browse through every file and extract what we want
    if len(files)==0:
        raise ValueError('No files selected')
    nsums = None
    all_dfs={}
    descriptions = []
    for nfile,file in enumerate(files):
        print(file)
        stack = StackFCS(file,load_stack=False)
        stack.load()
        
        description = stack.describe()
        descriptions.append(description)
        
        stack_res = stack.extract_results(ith=ith,
                                        chi_threshold=chi_threshold,use_mask=use_mask)
        diffs = stack_res["diffusion_coefficients"]
        chis = stack_res["non_linear_chis"]
        nmolecules = stack_res["number_molecules"]
        indices = stack_res['indices']
        square_errors = stack_res['square_errors']
        valid_fraction=stack_res["valid_fraction"]
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
            fname = description['filename']
            diff = diffs[nsum]
            chi = chis[nsum]
            nmol = nmolecules[nsum]
            square_error = square_errors[nsum]
            indice = indices[nsum].astype(int)
            repeats_arr = np.full(diff.size, nfile)
            name_arr = np.full(diff.size, fname)
            nsum_arr = np.full(diff.size, nsum)
            valid_fraction_arr = np.ones_like(diff)*valid_fraction[nsum]
            out_arr = np.array([name_arr,repeats_arr, diff, nsum_arr, chi, 
                                nmol,indice, square_error,valid_fraction_arr]).T
            
            df = pd.DataFrame(out_arr, columns = 
                              ["filename", "repeat","D [µm²/s]","binning",
                               "fit error", "N","label","square error","valid fraction"])
            df = df.astype({'filename':"str",
                           "repeat":"int",
                           "D [µm²/s]":"float",
                           "fit error":"float",
                           "N": "float",
                           "label":"int",
                           "square error": "float",
                           "valid fraction":"float"})
            all_dfs[nsum].append(df)
            
    parameters_dict = {"chi_threshold": chi_threshold,
                       "intensity_threshold":ith}
    # Extracts global parameters
    global_summaries = []
    global_names = []
    for nsum in nsums:
        all_dfs[nsum] = pd.concat(all_dfs[nsum])
        knm = "nsum {}".format(nsum)
        global_names.append(knm)
        sum_dict={"Median":np.median(all_dfs[nsum]["D [µm²/s]"].values),
                  "Mean": np.mean(all_dfs[nsum]["D [µm²/s]"].values),
                  "stdev":np.std(all_dfs[nsum]["D [µm²/s]"].values),
                  "valid fractions": np.mean(all_dfs[nsum]["valid fraction"].values)
                  }
        global_summaries.append(sum_dict)

    # Writes in excel file
    if not out_name.endswith(".xlsx"):
        out_name = out_name + ".xlsx"
    with pd.ExcelWriter(out_name) as writer:  
        dfs_total = []
        for nsum in nsums:
            dfpooled = all_dfs[nsum]
            dfs_total.append(summarise_df(all_dfs[nsum]))
            dfpooled.to_excel(writer, sheet_name = "nsum {}".format(nsum))
        dfs_total = pd.concat(dfs_total)
        dfs_total.to_excel(writer, sheet_name = "summaries file by file")
        dfs_description=pd.DataFrame.from_records(descriptions)
        dfs_description.to_excel(writer, sheet_name = "stack parameters")
        dfs_global = pd.DataFrame.from_records(global_summaries,index=global_names)
        dfs_global.to_excel(writer, sheet_name = "summary")
        df_pars = pd.DataFrame(parameters_dict, index=[0]).T
        df_pars.to_excel(writer, sheet_name = "export parameters")