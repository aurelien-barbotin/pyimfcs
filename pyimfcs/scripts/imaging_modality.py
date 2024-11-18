#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:26:34 2024

@author: aurelienb
"""

from pyimfcs.io import get_metadata_bioformats
import czifile
import glob
import re

def remove_tags(text):
    return re.sub('<[^<]+?>', '', text).lstrip(' ')

def get_acqmode(file):
    img = czifile.CziFile(file)
    metadata = img.metadata()
    description = metadata.split('\n')
    angle = float(remove_tags(list(filter(lambda x: 'Angle' in x, description))[0]))
    ffilter = remove_tags(list(filter(lambda x: '<Filter>' in x, description))[0])
    if "Transmission Light" in ffilter:
        return "BF"
    elif angle==0:
        return 'epi'
    else:
        return 'TIRF'
    
if __name__=='__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description=
         'Exports results stored in h5 files')
    
    parser.add_argument("-p", "--path", help="Target directory", default='.')
    parser.add_argument("-f", "--filter", 
                        help="Filters results: shows only selected imaging types. If not specified, shows all imaging types", 
                        default=None,
                        choices=['BF','epi','TIRF'])
    parser.add_argument('-r','--recursive', 
                        help='Boolean False by default. If True, analyses also files in all subfolders.',
                        default=False, type=bool)
    args = parser.parse_args()
    
    path=args.path
    filterval = args.filter
    files = sorted(glob.glob(path+"/*.czi"))
    if args.recursive:
        files = sorted(glob.glob(path+"/**/*.czi",recursive=args.recursive))
    print('{} files found'.format(len(files)))
    for file in files:
        fname = file[len(path):]
        if filterval is None:
            print("file",fname,':', get_acqmode(file))
        elif filterval==get_acqmode(file):
            print("file",fname,':', get_acqmode(file))
       
