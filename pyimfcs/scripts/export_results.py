#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 14:58:18 2023

@author: aurelienb
"""

import glob
import os

from pyimfcs.export import merge_fcs_results
"""
path="/home/aurelienb/Documents/Projects/Species/Staph/stat/2023_09_28"
print()
for x in os.walk(path):
    print(x)
    
subfolders = [ f.path for f in os.scandir(path) if f.is_dir() ]
"""
if __name__=='__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description=
         'Exports results stored in h5 files')
    
    parser.add_argument("-p", "--path", help="Target directory", default='.')
    
    parser.add_argument('-i','--intensity_threshold', 
                        help='intensity threshold',default=0.8, type=float)
    parser.add_argument('-c','--chi_threshold', 
                        help='quality metric threshold',default=0.015,type=float)
    parser.add_argument('-m','--use_mask', 
                        help='boolean, if True uses mask for export if found',
                        default=True, type=bool)
    
    parser.add_argument('-e','--export_folder', 
                        help='Optional argument. if specified, export all the\
excels in given export folder (renames them accordingly)',
                        default=None)
    parser.add_argument('-n','--savename', 
                        help='name of excel file to write. Overwritten if export_folder is provided',
                        default='results')
    
    parser.add_argument('-r','--recursive', 
                        help='Boolean True by default. If True, exports also results in all subfolders.',
                        default=True, type=bool)
    parser.add_argument('-d','--distances', 
                        help='''If provided, calculates in every stack the distance to a reference label (first argument),\
expecting n classes (2nd argument)''',
                        default=None, type=int,nargs='+')
    args = parser.parse_args()
    
    folder=args.path
    #print(folder)
    subfolders = glob.glob(folder+'/**/',recursive=args.recursive)
    #print(subfolders)
    for subfolder in subfolders:
        files = glob.glob(subfolder+'/*.h5')
        if len(files)>0:
            print('exporting subfolder ',subfolder[len(folder):])
            subfolder_name = subfolder[len(folder):].replace(os.sep,'_').lstrip('_').rstrip('_')
            if args.export_folder is not None:
                out_name = args.export_folder+"/"+subfolder_name+".xlsx"
            else:
                out_name = subfolder+"/"+args.savename+".xlsx"
            merge_fcs_results(out_name, files, 
                              ith=float(args.intensity_threshold),
                              chi_threshold=float(args.chi_threshold),
                              use_mask=args.use_mask, distance_labels=args.distances)
    