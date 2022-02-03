#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 09:15:05 2022

@author: aurelien
"""
import pandas as pd
import glob
import numpy as np

path="/home/aurelien/Documents/Analysis/SLB_analysis/POPC_SLBS/"
out_name = "test1.xlsx"

keep_single_indices = True #if true, gives an indes to single acquisition. If False, single index to an excel file.
files = glob.glob(path+"*.xlsx")

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
                
fileslist = [files, ['/home/aurelien/Documents/Analysis/SLB_analysis/DOPC_SLBs/2021_11_03.xlsx']]
merge_excels(fileslist,"test2.xlsx",conditions=['POPC', 'DOPC'])