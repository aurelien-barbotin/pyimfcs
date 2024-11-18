#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 11:28:24 2024

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
    meta_dict = {}
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
    
path="/run/user/1001/gvfs/smb-share:server=data.micalis.com,share=proced/microscopy/ZEISS/Aurelien/2024_08_01/1_NRNY/"

files = sorted(glob.glob(path+"/*.czi"))
print('{} files found'.format(len(files)))
for file in files:
    # dic = get_metadata_bioformats(file)
    
    img = czifile.CziFile(file)
    meta_dict = {}
    metadata = img.metadata()
    description = metadata.split('\n')
    # dirty hack
    acq_modes = list(filter(lambda x: '<AcquisitionMode>' in x, description))[:1]
    if len(acq_modes)>1:
        acq_mode = 'undefined'
    else:
        acq_mode = acq_modes[0]
    # acq_mode = 
    acq_mode = remove_tags(acq_mode)
    angle = remove_tags(list(filter(lambda x: 'Angle' in x, description))[0])
    ff = remove_tags(list(filter(lambda x: '<Filter>' in x, description))[0])
    
    if get_acqmode(file)=="epi":
        print('------')
        print(angle)
        print(ff)
        print("file",file.split('/')[-1],': ', get_acqmode(file))