#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 10:58:38 2023

@author: aurelienb
"""

from PyQt5.QtWidgets import (QDialog,QDialogButtonBox,QGridLayout,
                     QApplication, QPushButton, QComboBox)

from pyimfcs.fitting import fit_parameters_dict, fit_parameter_types

fit_par_types = sorted(list(fit_parameters_dict.keys()))

class ModelDialog(QDialog):
    def __init__(self):
        super().__init__()
        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        
        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
    
        self.loadPushButton = QPushButton("load")
        self.savePushButton = QPushButton("save")
        self.modelTypesComboBox = QComboBox()
        
        self.layout = QGridLayout()
        self.setLayout(self.layout)
        self.make_base_layout()
                
    def make_base_layout(self):
        self.layout.addWidget(self.buttonBox,0,0,1,2)
        self.layout.addWidget(self.savePushButton,1,0,1,1)
        self.layout.addWidget(self.loadPushButton,1,1,1,1)
        self.show()
        
    def accept(self):
        
        self.layout.addWidget(QPushButton("Push"), self.j, 0)
        self.show()
        super().accept()

    def load(self,fname):
        pass
    
    def save(self, fname):
        pass
    
    def reject(self):
        super().reject()

    def clear_layout(self):
        for i in reversed(range(self.layout.count())): 
            self.layout.itemAt(i).widget().setParent(None)
        self.make_base_layout()
        
app = QApplication([])
win = ModelDialog()
win.show()
app.exec_()
