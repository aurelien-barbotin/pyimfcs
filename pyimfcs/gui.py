#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 10:41:12 2022

@author: aurelienb
"""
import numpy as np
import glob
import os
import tifffile
import h5py
import sys
import json

from PyQt5 import QtCore
from PyQt5.QtWidgets import (QDialog, QWidget, QApplication,QListWidgetItem,
                             QPushButton, QLineEdit, QLabel, QSpinBox, QMessageBox)

from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QGridLayout,QGroupBox,QCheckBox,QHBoxLayout
from PyQt5.QtWidgets import QListWidget,QFileDialog, QComboBox, QDialogButtonBox

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.image import AxesImage

from pyimfcs.plotting import interactive_fcs_plot
from pyimfcs.class_imFCS import StackFCS
from pyimfcs.export import merge_fcs_results
from pyimfcs.process import batch_bacteria_process, get_metadata_zeiss
from pyimfcs.fitting import Fitter
from pyimfcs.splash_screen import LoadingWindow
from pyimfcs.make_masks import make_masks

BUNDLE_DIR = getattr(sys, '_MEIPASS', os.path.abspath(os.path.dirname(__file__)))

class ExperimentListWidget(QListWidget):
   """Class designed to contain the different correction rounds. Each correction
   element is stored in itemList. A correction item is a list of 
   [str modes,int number, str filename]"""
   to_analyse = QtCore.Signal(str)
   
   def __init__(self,*args,**kwargs):
       super().__init__(*args, **kwargs)
       #self.itemClicked.connect(self.Clicked)
       def itemChangedSlot(item):
           try:
               self.to_analyse.emit(item.data(QtCore.Qt.UserRole))
           except:
               self.to_analyse.emit(None)
       self.currentItemChanged.connect(itemChangedSlot)
      
   def addToList(self,file):
       item = QListWidgetItem( self)
       item.setText(os.path.split(file)[-1])
       item.setData(QtCore.Qt.UserRole, file)

   def fill(self,folder):
       index = self.currentRow()
       self.clear()
       files = glob.glob(folder+"/*.h5")
       files.sort()
       for file in files:
            self.addToList(file)
       try:
            self.setCurrentRow(index)
       except Exception as e:
            print(e)
            
   def get_filenames(self):
        items = []
        for x in range(self.count()):
            items.append(self.item(x).data(QtCore.Qt.UserRole))
        return items
    
class MatplotlibWindow(QDialog):
    plot_clicked = QtCore.Signal(int)
    back_to_plot = QtCore.Signal()
    image_clicked = QtCore.Signal(np.ndarray)
    onclickf = None
    onclick_function = lambda x: print('no thing to display')
    def __init__(self, parent=None):
        super().__init__(parent)

        # a figure instance to plot on
        self.figure = Figure(figsize = (20,15))

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

        self.canvas.mpl_connect('button_press_event', self.onclick)
        self.canvas.mpl_connect('pick_event', self.onpick)
        
    def make_axes(self,n=2):
        self.axes = [self.figure.add_subplot(2,1,1),
                      self.figure.add_subplot(2,1,2)]
        
    def onclick(self,event):
        if self.onclick_function is None:
            return
        self.onclick_function(event)
       
    def onpick(self,event):
        artist = event.artist
        if isinstance(artist, AxesImage):
            im = artist
            A = im.get_array()
            self.image_clicked.emit(A)
            
    def plot(self):
        self.canvas.draw()
        
class FCS_Visualisator(QWidget):
    onclick_function = None
    current_stack = None
    
    def __init__(self,*args,**kwargs):
        """This method is necessary innit? Mate"""
        super().__init__(*args, **kwargs)
        self.setAcceptDrops(True)
        self.folder = "."
        
        self.newAnalystButton = QPushButton("Load")
        self.refreshButton = QPushButton("Refresh")
        self.trashButton = QPushButton("Trash")
        self.exportButton = QPushButton("Export")
        self.processingPushButton = QPushButton("Process Files")
        self.maskPushButton = QPushButton('Set masks')
        
        self.newAnalystButton.clicked.connect(lambda :
            self.loadFiles(str(QFileDialog.getExistingDirectory(self, "Select Directory"))))
        self.refreshButton.clicked.connect(self.refreshFileList)
        self.trashButton.clicked.connect(self.trash_measurement)
        self.exportButton.clicked.connect(self.export_measurements)
        self.processingPushButton.clicked.connect(self.process_measurements)
        self.maskPushButton.clicked.connect(self.set_masks)
        
        self.current_mode = None
        self.expListWidget = ExperimentListWidget()
        self.plotBox = MatplotlibWindow()
        
        self.make_metrics_tab()
        # self.make_processing_tab()
        
        self.imageComparisonWidget = QWidget()
        self.imageComparisonWidgetOn=False
        self.connects()
        
        self.expListWidget.fill(".")
            
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(self.newAnalystButton,0,0,1,1)
        self.grid.addWidget(self.refreshButton,0,1,1,1)
        self.grid.addWidget(self.trashButton,0,2,1,1)
        self.grid.addWidget(self.exportButton,0,3,1,1)
        self.grid.addWidget(self.processingPushButton,0,4,1,1)
        self.grid.addWidget(self.maskPushButton,0,5,1,1)
        
        self.grid.addWidget(self.expListWidget,1,0,9,1)
        
        self.grid.addWidget(self.plotBox,1,1,10,10)
        self.grid.addWidget(self.metrics_tab,10,0,1,1)
        
        
    def dragEnterEvent(self, e):
        e.accept()
    
    def dropEvent(self,e):
        for url in e.mimeData().urls():
            url = str(url.toLocalFile())
            if url[-3:]==".h5" or url[-3:]=="npy" or url[-3:]=="SIN":
                url = "/".join(url.split("/")[:-1])
            self.loadFiles(url)
            
    def connects(self):
        #Connecting various signals
        self.expListWidget.to_analyse.connect(self.update_plot)
        self.plotBox.plot_clicked.connect(self.update_interactive_plot)
        self.plotBox.back_to_plot.connect(self.update_plot)
        
    def disconnects(self):
        self.expListWidget.to_analyse.disconnect(self.update_plot)
        self.plotBox.plot_clicked.disconnect(self.update_interactive_plot)
        self.plotBox.back_to_plot.disconnect(self.update_plot)

    def export_measurements(self):
        filename = str(QFileDialog.getSaveFileName(self, "Select File name", 
                                    os.path.join(os.getcwd(),"results.xlsx"),
                                    filter="*.xlsx")[0])
        if filename == "":
            return
        
        files = self.expListWidget.get_filenames()
        thr = None
        tht  = self.thresholdLineEdit.text()
        if tht.replace('.','',1).isdigit():
            thr = float(tht)
        
        ith = None
        intensity_threshold_tmp = self.intensityLineEdit.text()
        if intensity_threshold_tmp.replace('.','',1).isdigit():
            ith = float(intensity_threshold_tmp)
        use_mask=self.useMaskCheckBox.isChecked()
        print('Start exporting measurements ...')
        self.exportButton.setEnabled(False)
        msg1 = LoadingWindow()
        msg1.show()
        app.processEvents()
        merge_fcs_results(filename, files,
              ith = ith, chi_threshold = thr, 
              use_mask=use_mask)
        msg1.close()
        self.exportButton.setEnabled(True)
        print('Done exporting measurments')
    
    def process_measurements(self):
        fd = QFileDialog()
        fd.setAcceptDrops(True)
        
        folder = str(fd.getExistingDirectory(self, "Select Directory",os.getcwd()))
        if folder=="":
            return
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
                
                self.loadFiles(folder=folder)
                msg1.close()
                msg = QMessageBox()
                msg.setText('Processing Finished')
                msg.exec_()
                
    def set_masks(self):
        print('setting masks')
        folder = self.folder
        files = glob.glob(folder+"/*.h5")
        for file in files:
            maskname = file[:-3]+"_mask.tif"
            print(maskname)
            if os.path.isfile(maskname):
                mask=tifffile.imread(maskname)
                print("found mask for",file)
                # !!! Might struggle with filenames
                hf = h5py.File(file, 'a')
                if 'parameters/mask' in hf:
                    del hf['parameters/mask']
                hf['parameters/mask'] = mask
                hf.close()
                
    def trash_measurement(self):
        try:
            file = self.expListWidget.currentItem().data(QtCore.Qt.UserRole)
        except:
            return
        path,name = os.path.split(file)
        trashDir = path+"/trash"
        if not os.path.isdir(trashDir):
            os.mkdir(trashDir)
        os.rename(file,trashDir+"/"+name)
        self.refreshFileList()
        
    def update_plot(self,file=None,extract=False, load_stack = True):
        self.current_mode = None
        if file is None:
            try:
                file = self.expListWidget.currentItem().data(QtCore.Qt.UserRole)
            except:
                return 
        fig = None
        if not extract:
            fig = self.plotBox.figure
        self.plotBox.figure.clf()
        
        light_version= self.lightDisplayCheckBox.isChecked()
            
        if load_stack:
            print('Load stack')
            current_index = int(self.binningComboBox.currentIndex())
            self.current_stack = StackFCS(file, load_stack = False)
            self.current_stack.load(light_version =light_version)
            nsums = self.current_stack.fit_results_dict.keys()
            _ = self.update_binnings(list(nsums))
            new_index = 0
            if current_index<len(nsums):
                new_index = current_index
            self.binningComboBox.setCurrentIndex(new_index)
        nsum = int(self.binningComboBox.currentText())
        
        vmax = None
        vt  = self.vmaxLineEdit.text()
        if vt.replace('.','',1).isdigit():
            vmax = float(vt)
        
        vmin = None
        vt  = self.vminLineEdit.text()
        if vt.replace('.','',1).isdigit():
            vmin = float(vt)
            
        chi_thr = None
        tht  = self.thresholdLineEdit.text()
        if tht.replace('.','',1).isdigit():
            chi_thr = float(tht)
        
        ith = None
        intensity_threshold_tmp = self.intensityLineEdit.text()
        if intensity_threshold_tmp.replace('.','',1).isdigit():
            ith = float(intensity_threshold_tmp)
        use_mask = self.useMaskCheckBox.isChecked()
        self.plotBox.figure.clf()
        self.onclick_function = interactive_fcs_plot(self.current_stack, fig=fig, 
                        nsum = nsum, vmax=vmax, vmin=vmin ,chi_threshold=chi_thr, 
                        light_version=light_version, ith=ith,use_mask=use_mask)
        self.plotBox.onclick_function = self.onclick_function
        
          
    def update_interactive_plot(self,mode):
        self.plotBox.figure.clf()
        self.plotBox.plot()
    
    def update_binnings(self,nsums):

        self.binningComboBox.disconnect()
        self.binningComboBox.clear()
        self.binningComboBox.addItems([str(w) for w in nsums])
        self.binningComboBox.currentIndexChanged.connect(lambda x: self.update_plot(load_stack=False))
        return 0
    
    def loadFiles(self,folder=None):
        if folder is None:
            folder = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        if folder:
            self.disconnects()
            self.expListWidget.fill(folder)
            self.connects()
            self.folder = folder
            os.chdir(folder)
            
    def refreshFileList(self):
        self.loadFiles(self.folder)
        self.update_plot()
        
    def make_metrics_tab(self):
        top = QGroupBox('Display options')
        toplay = QGridLayout()
        top.setLayout(toplay)
        
        self.binningComboBox = QComboBox(self)
        self.binningComboBox.addItem("2")
        self.binningComboBox.addItem("4")
        self.binningComboBox.currentIndexChanged.connect(lambda x: self.update_plot(load_stack=False))
        
        self.vmaxLineEdit = QLineEdit("None")
        self.vmaxLineEdit.editingFinished.connect(lambda : self.update_plot(load_stack=False))
        
        self.vminLineEdit = QLineEdit("None")
        self.vminLineEdit.editingFinished.connect(lambda : self.update_plot(load_stack=False))
        
        self.thresholdLineEdit = QLineEdit("0.015")
        self.thresholdLineEdit.editingFinished.connect(lambda : self.update_plot(load_stack=False))
        
        self.intensityLineEdit = QLineEdit("0.8")
        self.intensityLineEdit.editingFinished.connect(lambda : self.update_plot(load_stack=False))
        
        self.lightDisplayCheckBox = QCheckBox("Light version")
        self.lightDisplayCheckBox.setChecked(True)
        self.lightDisplayCheckBox.toggled.connect(lambda : self.update_plot(load_stack= not self.lightDisplayCheckBox.isChecked()))
        
        self.useMaskCheckBox = QCheckBox("Use masks")
        self.useMaskCheckBox.setChecked(True)
        self.useMaskCheckBox.toggled.connect(lambda : self.update_plot(load_stack=False))
        
        toplay.addWidget(QLabel("Binning"),0,0)
        toplay.addWidget(self.binningComboBox,0,1)
        toplay.addWidget(self.vmaxLineEdit,1,1)
        toplay.addWidget(QLabel("Max diff. shown"),1,0)
        toplay.addWidget(self.vminLineEdit,2,1)
        toplay.addWidget(QLabel("Min diff. shown"),2,0)
        
        toplay.addWidget(self.thresholdLineEdit,3,1)
        toplay.addWidget(QLabel("Chi threshold"),3,0)
        toplay.addWidget(self.intensityLineEdit,4,1)
        toplay.addWidget(QLabel("Intensity threshold (0-1)"),4,0)
        toplay.addWidget(self.lightDisplayCheckBox,5,0,1,2)
        toplay.addWidget(self.useMaskCheckBox,6,0,1,2)
        self.metrics_tab = top
        
    def make_processing_tab(self):
        """Unused yet. Maybe later"""
        top = QGroupBox('Fitting options')
        toplay = QHBoxLayout()
        top.setLayout(toplay)
        name = "test"
        self.changeFitterPushButton = QPushButton("Change fit function")
        toplay.addWidget(QLabel('Fitter name: '))
        toplay.addWidget(self.changeFitterPushButton)
        self.processing_tab = top
        
class ParametersDialog(QDialog):

    def __init__(self, files):
        super().__init__()
        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        
        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.find_models()
        self.load_params()
        
        if len(files)>0:
            file = files[0]
            try:
                self.dt, self.xscale, self.yscale = get_metadata_zeiss(file)
            except:
                print('Data not loaded')
                pass
            
        self.layout = QGridLayout()
        self.first_nSpinBox = QSpinBox()
        self.first_nSpinBox.setMaximum(10**5)
        self.first_nSpinBox.setValue(self.first_n)
        self.last_nSpinBox = QSpinBox()
        self.last_nSpinBox.setMaximum(10**5)
        self.last_nSpinBox.setValue(self.last_n)
        self.nsumsLineEdit = QLineEdit(str(self.nsums)[1:-1])
        self.dtLineEdit = QLineEdit("{:.2f}".format(self.dt*10**3))
        self.psizeLineEdit = QLineEdit("{:.4f}".format(self.xscale))
        
        self.registrationSpinBox = QSpinBox()
        self.registrationSpinBox.setMaximum(10**5)
        self.registrationSpinBox.setMinimum(0)
        self.registrationSpinBox.setValue(self.nreg)
        
        self.makeMasksCheckBox = QCheckBox("Make masks")
        self.makeMasksCheckBox .setChecked(self.make_masks_bool)
        
        self.masksDiameterLineEdit = QLineEdit(str(self.masks_diameter))
        self.masksDiameterLineEdit.setText(str(self.masks_diameter))
        
        shift_row= 1 # to move everything up or down if needed
        self.layout.addWidget(QLabel('Fitting model'), shift_row-1, 0, 1, 1)
        self.layout.addWidget(self.modelsComboBox, shift_row-1, 1, 1, 1)
        self.layout.addWidget(QLabel('Remove first n frames'),0+shift_row,0)
        self.layout.addWidget(self.first_nSpinBox,0+shift_row,1)
        self.layout.addWidget(QLabel('Remove Last n frames'),1+shift_row,0)
        self.layout.addWidget(self.last_nSpinBox,1+shift_row,1)
        self.layout.addWidget(QLabel('Binning values'),2+shift_row,0)
        self.layout.addWidget(self.nsumsLineEdit,2+shift_row,1)
        self.layout.addWidget(QLabel('Frame interval (ms)'),3+shift_row,0)
        self.layout.addWidget(self.dtLineEdit,3+shift_row,1)
        self.layout.addWidget(QLabel('Pixel size (µm)'),4+shift_row,0)
        self.layout.addWidget(self.psizeLineEdit,4+shift_row,1)
        self.layout.addWidget(QLabel('Registration pooling value'),5+shift_row,0)
        self.layout.addWidget(self.registrationSpinBox,5+shift_row,1)
        
        self.layout.addWidget(self.makeMasksCheckBox,6+shift_row,0,1,2)
        self.layout.addWidget(QLabel('Expected mask diameter (µm):'),7+shift_row,0)
        self.layout.addWidget(self.masksDiameterLineEdit)
        
        self.layout.addWidget(self.buttonBox,10,0,1,2)
        self.setLayout(self.layout)
    
    def update_params(self):
        self.load_model()
        self.default_parameter_dialog['fit_model_name'] = self.modelsComboBox.currentText()
        
        self.first_n = self.first_nSpinBox.value()
        self.default_parameter_dialog['first_n'] = self.first_n
        
        self.last_n = self.last_nSpinBox.value()
        self.default_parameter_dialog['last_n'] = self.last_n
        
        nsums= self.nsumsLineEdit.text()
        self.nsums = [int(w) for w in nsums.split(",")]
        self.default_parameter_dialog['nsums']=self.nsums
        
        self.nreg= self.registrationSpinBox.value()
        self.default_parameter_dialog['nreg'] = self.nreg
        
        self.dt = float(self.dtLineEdit.text())*10**-3 # conversion in s
        self. default_parameter_dialog['dt'] = self.dt #reconversion in ms for saving
        
        self.psize = float(self.psizeLineEdit.text())
        self.default_parameter_dialog['psize'] = self.psize
        
        self.xscale=self.psize
        self.default_parameter_dialog['xscale'] = self.xscale
        
        self.yscale=self.psize
        self.default_parameter_dialog['yscale'] = self.yscale
        
        self.make_masks_bool = self.makeMasksCheckBox.isChecked()
        self.default_parameter_dialog['make_masks_bool'] = self.make_masks_bool
        self.default_parameter_dialog['masks_diameter'] = self.masks_diameter
        self.save_params()
      
        
    def find_models(self):
        self.models = glob.glob(os.path.join(BUNDLE_DIR,"models/*.json"))
        self.modelnames = [os.path.split(w)[-1][:-5] for w in self.models]
            
        self.modelsComboBox = QComboBox()
        for mname in self.modelnames:
            self.modelsComboBox.addItem(mname)
            
    def load_model(self):
        """Loads the model currently selected stored in a json file"""
        current_index = int(self.modelsComboBox.currentIndex())
        model = self.models[current_index]
        
        with open(model, "r") as f:
            self.model_parameter_dict = json.load(f)
            
    def set_model_combobox_fromname(self,name):
        try:
            index = self.modelnames.index(name)
        except:
            print('Saved model was not found')
            index = 0
        self.modelsComboBox.setCurrentIndex(index)
        
    def load_params(self):
        print('loading params')
        try:
            with open(os.path.join(BUNDLE_DIR,"data/default_parameters_fordialog.json")
                      ,"r") as f:
                self.default_parameter_dialog = json.load(f)
        except:
            self.default_parameter_dialog = {"first_n":15000,
                                         "last_n":0,
                                         "nsums":[2,3],
                                         "nreg":4000,
                                         "dt":1,
                                         "xscale":1,
                                         "yscale":1,
                                         "make_masks_bool":False,
                                         "masks_diameter":1}
            
        self.first_n = self.default_parameter_dialog['first_n']
        self.last_n = self.default_parameter_dialog['last_n']
        self.nsums = self.default_parameter_dialog['nsums']
        self.nreg = self.default_parameter_dialog['nreg']
        self.dt = self.default_parameter_dialog['dt']
        self.xscale=self.default_parameter_dialog['xscale']
        self.yscale=self.default_parameter_dialog['yscale']
        self.set_model_combobox_fromname(self.default_parameter_dialog['fit_model_name'])
        
        # legacy of previous versions when these keys were not part of the program
        try:
            self.make_masks_bool = self.default_parameter_dialog['make_masks_bool']
            self.masks_diameter = self.default_parameter_dialog['masks_diameter']
        except:
            self.make_masks_bool = False
            self.masks_diameter = 1
            self.default_parameter_dialog['make_masks_bool'] = self.make_masks_bool
            self.default_parameter_dialog['masks_diameter'] = self.masks_diameter
            
    def save_params(self):
        with open(os.path.join(BUNDLE_DIR,"data/default_parameters_fordialog.json")
                  ,"w") as f:
            json.dump(self.default_parameter_dialog,f,indent=2)
            
    def accept(self):
        self.update_params()
        super().accept()
        
app = QApplication([])
#files = glob.glob("/home/aurelienb/Data/2022_09_22/*.tif")
#win = ParametersDialog(files)

win = FCS_Visualisator()
win.show()
app.exec_()
