#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 18:17:27 2022

@author: aurelien

"""
import numpy as np
import pandas as pd

from pyimfcs.class_imFCS import StackFCS
from pyimfcs.fitting import possible_fit_parameters

def summarise_df(df, val_keys=['D [µm²/s]']):
    """Merges results of a Dataframe containing multiple experiments"""
    # General stuff independent of fitting model
    repeat = df.groupby('repeat')['repeat'].last()
    filename = df.groupby('repeat')['filename'].last()
    binning = df.groupby('repeat')['binning'].last()
    valid_fractions = df.groupby("repeat")["valid_fraction"].mean()
    
    out_dict={
            "filename": filename,
            "binning": binning,
            "repeat":repeat,
            "valid_fraction": valid_fractions,
            }
    # summarises every parameter
    for val_key in val_keys:
        out_dict[val_key+" mean"] = df.groupby('repeat')[val_key].mean()
        out_dict[val_key+" std"] = df.groupby('repeat')[val_key].std()
        out_dict[val_key+" median"] = df.groupby('repeat')[val_key].median()
        out_dict[val_key+"_end"]="" # stupid hack for better vis
    out_dict["count"] = df.groupby('repeat')[val_key].count()
    out_df = pd.DataFrame.from_dict(out_dict)
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
        res_keys=list(stack_res.keys())
        # check that nsums are identical in all files
        nsums_tmp = sorted(stack_res[res_keys[0]].keys())
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
            out_dict = dict(zip(stack_res.keys(),[stack_res[w][nsum] for w in stack_res.keys()]))
            out_dict['filename'] = description['filename']
            out_dict['binning'] = nsum
            out_dict['repeat'] = nfile
            df = pd.DataFrame.from_dict(out_dict)
            
            # specifies types in the Dataframe
            default_type="float"
            specific_types={'filename':"str",
                           "repeat":"int",
                           "label":"int",}
            types_dict={}
            for k in out_dict.keys():
                if k in specific_types.keys():
                    types_dict[k] = specific_types[k]
                else:
                    types_dict[k] = default_type
            all_dfs[nsum].append(df)
            
    parameters_dict = {"chi_threshold": chi_threshold,
                       "intensity_threshold":ith}
    # Extracts global parameters
    global_summaries = []
    global_names = []
    parameters_to_summarise = list(filter(
        lambda x: x in possible_fit_parameters,all_dfs[nsums[0]][0].keys()))
    for nsum in nsums:
        all_dfs[nsum] = pd.concat(all_dfs[nsum])
        knm = "nsum {}".format(nsum)
        global_names.append(knm)
        sum_dict={}
        for pp in parameters_to_summarise:
            sum_dict[pp+" Median"] = np.median(all_dfs[nsum][pp].values)
            sum_dict[pp+" mean"] = np.mean(all_dfs[nsum][pp].values)
            sum_dict[pp+" std"] = np.std(all_dfs[nsum][pp].values)
            sum_dict[pp+" end"] = ""
        sum_dict[pp+" valid fractions"] = np.mean(all_dfs[nsum]["valid_fraction"].values)
        global_summaries.append(sum_dict)

    # Writes in excel file
    if not out_name.endswith(".xlsx"):
        out_name = out_name + ".xlsx"
    with pd.ExcelWriter(out_name) as writer:
        dfs_total = []
        for nsum in nsums:
            dfpooled = all_dfs[nsum]
            dfs_total.append(summarise_df(all_dfs[nsum],val_keys=parameters_to_summarise))
            dfpooled.to_excel(writer, sheet_name = "nsum {}".format(nsum))
        dfs_total = pd.concat(dfs_total)
        dfs_total.to_excel(writer, sheet_name = "summaries file by file")
        dfs_description=pd.DataFrame.from_records(descriptions)
        dfs_description.to_excel(writer, sheet_name = "stack parameters")
        dfs_global = pd.DataFrame.from_records(global_summaries,index=global_names)
        dfs_global.to_excel(writer, sheet_name = "summary")
        df_pars = pd.DataFrame(parameters_dict, index=[0]).T
        df_pars.to_excel(writer, sheet_name = "export parameters")