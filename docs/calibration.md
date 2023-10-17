# Calibration of a fitting model

In conventional point FCS, the only assumption needed to fit data is that the observation area is Gaussian (at least in 2D). In imaging FCS, this is different: the detection area is a convolution of the pixel area by the point spread function (PSF). To account for this, the size of the detection PSF is a parameter of the standard imFCS 2D diffusion fitting model (see ref [1]).

## Defining a new fitting model
In pyimfcs, fitting models are stored in a '.json' file in the folder 'pyimfcs/models/'. To create a fitting model from an existing template using the parameters of your experiments, you need to create a json file containing the necessary information. It is recommended to use one of the files already existing as a template. For example, if you create a file 'default_model.json' in 'pyimfcs/models/' containing this:

	{
	  "mtype": "2D"
	  "sigma": 0.19,
	  "ginf": true,
	}

pyimfcs will add a fitting model of type '2D', with a PSF is standard deviation 0.19 µm, and that accepts convergence towards non-zero values (parameter 'ginf'). The possible model types can be found as keys to the dictionary 'fit_functions' in 'pyimfcs/fitting.py'. They currently include: "2D", "3D","2D_2c" (2 component diffusion in 2D), "2D_anisotropic" (2D diffusion in e.g a rod, see ref. [2]) "2D_spherical" (isotropic diffusion in a sphere of small dimensions as compared to the PSF size, see also ref. [2]).

## Implementation of a new fitting class

This requires minimum Python programming skills. If you are unhappy with current classes (e.g if you want to implement a 3-components fitting model), you can create a new fitting class by editing the file 'pyimfcs/fitting.py':
- Implement the fitting function maker (use for instance gim2D as a template)
- Add this function maker to the fit_functions dictionary. Example if you want to add 3-component diffusion model, add to the dictionary the following key:value pair: '"2D_3C":gim2D_3components'
- Add the parameters required by the fit function maker to the dictionary fit_parameters_dict
- Specify the names of the fitting parameters by updating the dictionary fit_result_names_dict, with the key:value pair "parameter name":position in output


## References
[1]  Ries, J., E.P. Petrov, and P. Schwille. 2008. Total Internal Reflection Fluorescence Correlation Spectroscopy: Effects of Lateral Diffusion and Surface-Generated Fluorescence. Biophysical Journal. 95:390–399.
[2] Barbotin, A.; Billaudeau, C.; Sezgin, E.; Carballido Lopez, R. Quantification of Membrane Fluidity in Bacteria Using TIR-FCS; preprint; Biophysics, 2023. https://doi.org/10.1101/2023.10.13.562271.

