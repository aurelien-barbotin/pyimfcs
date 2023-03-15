#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 16:53:05 2023

@author: aurelienb
"""

from pyimfcs.class_imFCS import StackFCS

path='/run/user/1000/gvfs/smb-share:server=data.micalis.com,share=proced/microscopy/ZEISS/Aurelien//2023_03_14/2_44_37/2p5X/summaries/Image 49.h5'

stack = StackFCS(path,load_stack=False)
stack.fcs_curves_dict={"test":0.5}
stack.load()