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

- Fitting model: Choose among pre-defined (see [calibration](https://github.com/aurelien-barbotin/pyimfcs/blob/main/docs/calibration.md)) fitting models for FCS curves
- Remove first n frames: remove n frames at the beginning of each stack. Set to 0 to keep all frames
- Remove last n frames: remove n frames at the end of each stack. Set to 0 to keep all frames
- Binning values: groups pixels by binning values before processing. It is possible to use several binning values, separated by a coma.
- Frame interval (ms): the time between frames that pyimfcs will use by **default**. This value will be replaced by the value found in the stack metadata, if any (for now, only metadata from Zeiss image files can be read).
- Pixel size (µm): similarly, the value pyimfcs will use by default.
- Registration pooling value: if there is drift during acquisition, pyimfcs an use an image registration algorithm to compensate. It will find shifts between frames averaged in time (to have sufficient signal). If this value **n** is set to >0, frames will be averaged **n*** by **n**, the shift between these averaged frames will be calculated and the whole stack will be inversely shifted in order to correct for the drift. If set to 0, no registration is performed.
- Process files in subfolders: if checked, pyimfcs will find all subfolders in the current folder and process the image stacks it finds. If not, processes only in the selected folder
- Make masks: If checked, will generate masks from the averaged image stack using a Watershed algorithm. The expected mask diameter parameter is then used to separate adjacent objects based on their size (e.g if Expected mask diameter is set to 1 µm, the masks created will find two objects if they are adjacent only if they are separated by more than 1 µm).

## Data visualisation and export

### Visualisation
Once the processing is done, it is possible (and recommended !) to inspect the processed dataset. Processed stacks are stored on a file with a h5 extension, bearing the same name as the original stack. For instance, the processing results of 'stack1.tif' are stored in a file 'stack1.h5'. You can inspect any processed stack by clicking on its name in the list on the top left of the GUI.

![A screenshot of the pyimfcs GUI with a dataset open](https://github.com/aurelien-barbotin/imFCS/blob/main/images/manual/4_data_processed.png)

Different display options can be found on the bottom left of the GUI. 
- Binning: Select the binning value for which you want to display the results.
- Max diff. shown: upper limit on diffusion coefficient you want to display. All values above will be discarded (shown in white in the top right panel). If you do not want to discard any outlier, set this to None
- Min diff. shown: lower limit on diffusion coefficient you want to display. Similarly, set to None to ignore
- Chi threshold: Maximum allowed value for the fit error metric. Any curves with higher fit error will be discarded. Set to a high value (e.g. 1) if you do not want to filter curves based on fit quality (not recommended). Recommended thresholds are between 0.015-0.03
- Intensity threshold: excludes any pixels from the analysis that are below (Intensity threshold)*(max value). Max value is defined as 98th percentile of the image intensity to limit the effect of outliers. Set to None to not use an intensity threshold. The selected pixels will be displayed within red lines in the 'Intensity' panel
- Light version: untick if you want to see more information on screen (mainly intensity timetrace at each point). Useful to understand if something went wrong but makes display much slower
- Use masks: tick if you want to account for masks (if any found in the h5 files).

If you want to use your own masks, but are unhappy with the masks generated automatically, you can create them manually. They have to be 2d images of the same size as the original stack (e.g if the image stack is 50000 frames x 128 x 10 pixels, the mask must have a size of 128x10 pixels). Masks contain only integer values, with pixels having the same value belonging to the same object. To create a mask for the file 'stack1.tif', place it in the same folder under the name 'stack1_mask.tif' then click on 'Set masks' in the GUI.

If a mask is present in the processed dataset and the 'Use mask' option is selected, the outlines of different masks will be shown in different colors in the 'Intensity' panel.

### Export

Results are exported folder by folder. To export all the results within a folder, load these results in the GUI of pyimfcs. You can then click on 'export' to generate an excel file containing all the results.

## Naming convention for output files

pyimfcs can associate files of raw data (tif files), processed data (h5 files) and masks (tif files ending with _mask.tif) based on filenames. A stack named **stack1.tif** is associated with the following files:
- stack1_mask.tif : the mask
- stack1.h5 : the processed data

## Scripts (advanced users)

I wrote some scripts that were frequently needed so you don't have to. They are located in pyimfcs/scripts. The script named 'export_results.py' can ba called from the command line and is used to generate an excel file summarising the results within a folder or a series of folders.

