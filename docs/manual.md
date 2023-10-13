# pyimfcs User's manual

This document explains how to use the package pyimfcs once the package is installed and the setup is calibrated, that is that a fitting model accounting for PSF size was derived. For calibration and setting up the fitting model, see [calibration](../blob/main/docs/calibration.md).

Assuming that the package was properly installed, let's first start a session of the GUI of pyimfcs:

	python -m pyimfcs.gui

If you installed pyimfcs in a virtual environment, please don't forget to activate it first!

## Processing data

In this step, pyimfcs will:
- Take as an input a series of tiff image stacks (ordered as follows : t,x,y)
- Perform correlations within desired binning value(s)
- fit the corresponding curves
- (optional) create segmentation masks within the image and save them in the folder where the dataset is located
- Fetch masks in the folder, if there are any, and assign them to the corresponding image stack
- Save each processed stack in a h5 file. If the raw data stack was named 'stack1.tif', the processed data will be saved as 'stack1.h5'

Here is how to process data in details:

### Open the GUI

As described in [README](../blob/main/README.md). Once open the software looks like this:

![A screenshot of the graphical interface of pyimfcs](https://github.com/aurelien-barbotin/imFCS/blob/main/images/manual/1_opening.png)

### Select a folder with 3D stacks to process

Click on 'Process Files' at the top of the GUI. This will open a file dialog. Browse towards the folder you want to process. It can directly contain 3D tif stacks to process, or contain multiple subfolders themselves containing tif stacks. 

![A screenshot of an image stack for FCS, opened with ImageJ](https://github.com/aurelien-barbotin/imFCS/blob/main/images/manual/2_stack_example.png)

### Select processing settings

Once you browsed to the right folder, click on 'Open'. This will prompt the opening of a dialog in which you can select the processing options

![A screenshot of a the processing dialog](https://github.com/aurelien-barbotin/imFCS/blob/main/images/manual/3_processing_dialog.png)

- Fitting model: Choose among pre-defined (see [calibration](../blob/main/docs/calibration.md)) fitting models for FCS curves
- Remove first n frames: remove n frames at the beginning of each stack. Set to 0 to keep all frames
- Remove last n frames: remove n frames at the end of each stack. Set to 0 to keep all frames
- Binning values: groups pixels by binning values before processing. It is possible to use several binning values, separated by a coma.
- Frame interval (ms): the time between frames that pyimfcs will use by **default**. This value will be replaced by the value found in the stack metadata, if any (for now, only metadata from Zeiss image files can be read).
- Pixel size (µm): similarly, the value pyimfcs will use by default.
- Registration pooling value: if there is drift during acquisition, pyimfcs an use an image registration algorithm to compensate. It will find shifts between frames averaged in time (to have sufficient signal). If this value **n** is set to >0, frames will be averaged **n*** by **n**, the shift between these averaged frames will be calculated and the whole stack will be inversely shifted in order to correct for the drift. If set to 0, no registration is performed.
- Process files in subfolders: if checked, pyimfcs will find all subfolders in the current folder and process the image stacks it finds. If not, processes only in the selected folder
- Mkae masks: If checked, will generate masks from the averaged image stack using a Watershed algorithm. The expected mask diameter parameter is then used to separate adjacent objects based on their size (e.g if Expected mask diameter is set to 1 µm, the masks created will find two objects if they are adjacent only if they are separated by more than 1 µm).

## Data visualisation and export



## Names of output files

## Scripts (advanced users)

I wrote some scripts that were frequently needed so you don't have to.

