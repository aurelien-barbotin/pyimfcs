#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 14:56:04 2023

@author: aurelienb
"""
from PyQt5.QtWidgets import QApplication, QFileDialog,QMessageBox
from pyimfcs.gui import FCS_Visualisator, ParametersDialog, LoadingWindow

import glob
import os
from pyimfcs.make_masks import make_masks
from pyimfcs.process import batch_bacteria_process

global app
app = QApplication([])
#win = FCS_Visualisator()
#win.process_measurements()

fd = QFileDialog()
fd.setAcceptDrops(True)

folder = str(fd.getExistingDirectory(None, "Select Directory",os.getcwd()))
if folder=="":
    pass

files = glob.glob(folder+"/*.tif")
if len(files)>0:
    pdial = ParametersDialog(files)
    if pdial.exec():
        msg1 = LoadingWindow()
        app.processEvents()
        msg1.show()
        app.processEvents()
        parameters_dict = pdial.model_parameter_dict
        
        if pdial.make_masks_bool:
            for file in files:
                print('Generating mask for file: ',file)
                make_masks(file,psize=pdial.psize,
                           celldiam=pdial.masks_diameter)
        parameters_dict["a"] = pdial.psize
        # fitter = Fitter(parameters_dict)
        
        batch_bacteria_process(files, first_n = pdial.first_n,
                               last_n = pdial.last_n,nsums = pdial.nsums,
                               nreg = pdial.nreg, default_dt = pdial.dt, 
                               default_psize = pdial.psize, 
                               default_fitparams=parameters_dict)
        
        msg1.close()
        msg = QMessageBox()
        msg.setText('Processing Finished')
        msg.exec_()
        

app.exec_()
