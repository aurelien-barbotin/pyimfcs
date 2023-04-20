PyimFCS
=======

A package to analyse imFCS data.It features, in a graphical user interface:
- calculation of correlation curves in an image or a series of images, for each pixel and for different resampling values
- Fitting (simple 2D model only for now) and extraction of diffusion coefficients
- Filtering of curves based on intensity, fitting quality, and position (possibility to use a mask)
- Saving/Loading processed datasets in hdf5 format
- Extraction of data in a human-readable excel file.

Installation:
-------------
pyimfcs was last tested on python 3.10 but should work with many other python versions. To install, simply run:
`pip install -e path_to_pyimfcs`
where `path_to_pyimfcs` is the path to the folder containing `setup.py`

Usage:
------
It is recommended to use pyimfcs with its gui, for this simply run

	python -m pyimfcs.gui
