# PyimFCS

A package to analyse imFCS data.It features, in a graphical user interface:
- calculation of correlation curves in an image or a series of images, for each pixel and for different resampling values
- Fitting (simple 2D model only for now) and extraction of diffusion coefficients
- Filtering of curves based on intensity, fitting quality, and position (possibility to use a mask)
- Saving/Loading processed datasets in hdf5 format
- Extraction of data in a human-readable excel file.

![A screenshot of the graphical interface of pyimfcs](https://github.com/aurelien-barbotin/imFCS/blob/main/images/screenshot_imFCS.png)

## Installation:

pyimfcs was last tested on python 3.10 but should work with many other python versions. To install, simply run:

	pip install -e path_to_pyimfcs

where `path_to_pyimfcs` is the path to the folder containing `setup.py`

## Usage:

It is recommended to use pyimfcs with its gui, for this simply run

	python -m pyimfcs.gui

A user manual can be found in [docs/manual.md](../blob/main/docs/manual.md)

Using pyimfcs on one's custom setup requires setting up fitting models. The procedure to do this can be found in [docs/calibration.md](../blob/main/docs/calibration.md)

## Credit:


If you used this package for your research, please cite: 

Barbotin, A.; Billaudeau, C.; Sezgin, E.; Carballido Lopez, R. Quantification of Membrane Fluidity in Bacteria Using TIR-FCS; preprint; Biophysics, 2023. https://doi.org/10.1101/2023.10.13.562271.

