# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 14:57:50 2021

@author: abarbotin
"""
import tifffile
import multipletau
import numpy as np
import matplotlib.pyplot as plt

import h5py
import os

from skimage.filters import threshold_otsu

from pyimfcs.shift_correction import stackreg
from pyimfcs.io import get_image_metadata
from pyimfcs.metrics import new_chi_square, intensity_threshold
from pyimfcs.methods import downsample_mask, indices_intensity_threshold
from pyimfcs.fitting import Fitter

# dictionary of old names used in StackFCS. Keys are old names, values new ones
legacy_dic_names = {"fcs_curves_dict":"correlations", 
                    "traces_dict":"traces" ,
                    "fit_results_dict":"parameters_fits", 
                    "yh_dict":"yhat", 
                    "thumbnails_dict":"thumbnails", 
                    "metadata_dict":"metadata",
                    "fitting_parameters_dict":"fitting_parameters_dict"}

class StackFCS(object):
    # Dictionaries to save and load
    dic_names = ["fcs_curves_dict", "traces_dict", "fit_results_dict", 
                 "yh_dict", "thumbnails_dict", "metadata_dict","fitting_parameters_dict"]
    # parameters to save
    parameters_names = ["dt", "xscale", "yscale", "path", "nreg", "shifts",
                        "first_n", "last_n", "clipval", "bl_kernel_size","mask"]
    
    default_psize=1
    default_dt=1
    def __init__(self, path, background_correction=True,
                 blcorrf=None, first_n=0, last_n=0, fitter=None, dt=None,
                 remove_zeroes=False, clipval=0, load_stack=True):
        
        self.load_stack = load_stack
        self.first_n = first_n
        self.last_n = last_n
        self.path = path
   
        self.mask = None
        if load_stack:
            self.stack = tifffile.imread(path)
            if len(self.stack.shape) == 2:
                raise ValueError('imFCS data must be 3-dimensional')
            self.stack = self.stack[self.first_n:self.stack.shape[0] - self.last_n]
            maskpath = path.rstrip('.tif')+"_mask.tif"
            if os.path.isfile(maskpath):
                print('Mask found at ',maskpath)
                self.mask = tifffile.imread(maskpath)
        else:
            self.stack = np.zeros((10, 50, 50))

        self.fitter = fitter

        self.threshold_map = None
        self.intensity_threshold_f = intensity_threshold # the intensity threshold function
        
        # removes clipval points before and after the intensity timetrace
        # before correlation. To remove artefacts from bleaching correction
        # or image registration
        self.clipval = clipval

        if background_correction:
            self.stack = self.stack - self.stack.min()

        self.blcorrf = blcorrf
        self.bl_kernel_size = 'None'
        # shift correction
        self.nreg = 0
        self.shifts = np.zeros(1)
        # resuts dictionaries
        self.fcs_curves_dict = {}
        self.traces_dict = {}
        self.fit_results_dict = {}
        self.yh_dict = {}
        self.chisquares_dict = {}
        self.thumbnails_dict = {}
        self.square_err_dict = {}
        self.fitting_parameters_dict = {}
        self.distance_to_masks_dict = {}
        
        self.metadata_dict = {}
        self.metadata_fully_loaded = False
        if dt is None and load_stack:
            try:
                metadata = get_image_metadata(path)
                self.metadata_dict = metadata
                dt = metadata['finterval']
                if 'Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1' in metadata:
                    xscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'] * 10 ** 6
                    yscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1'] * 10 ** 6
                else:
                    xscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX'] * 10 ** 6
                    yscale = metadata['Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY'] * 10 ** 6
                self.xscale = xscale
                self.yscale = yscale
                if xscale != yscale:
                    raise ValueError('Not square pixels')
                self.metadata_fully_loaded = True
            except:
                print('error loading metadata')
                dt = self.default_dt
                self.xscale = self.default_psize
                self.yscale = self.default_psize
        else:
            self.xscale = self.default_psize
            self.yscale = self.default_psize
            dt = self.default_dt

        self.dt = dt

        if remove_zeroes:
            print("Achtung! Removing zeroes in the intensity timetrace may lead to artefacts!")
            trace = self.stack.sum(axis=(1, 2))
            if np.any(trace == 0):
                print('frames removed')
                self.stack = self.stack[trace != 0, :, :]

    def save(self, name=None,exclude_list=[]):
        """Saves StackFCS as an h5 file.
        Parameters:
            name: name of target file. If None, named after the original tif file.
            exclude_list (list): name of dictionaries (from dic_names) that should not be saved"""
        if name is None:
            name = os.path.splitext(self.path)[0] + ".h5"
        if not name.endswith(".h5"):
            name += ".h5"
        if os.path.isfile(name):
            print("Removing existing file with same name")
            os.remove(name)
            
        h5f = h5py.File(name, "w")

        for dname in self.dic_names:
            if dname in exclude_list:
                continue
            dic = getattr(self, dname)
            for key, item in dic.items():
                if type(item)==bytes:
                    try:
                        item = item.decode('utf-8')
                    except:
                        pass
                if type(item)!=list and type(item)!=bytes:
                    h5f[dname + "/" + str(key)] = item
                else:
                    print("not encoded")
        for pn in self.parameters_names:
            par =  getattr(self, pn)
            if par is not None:
                h5f["parameters/" + pn] = par
        if self.blcorrf is not None:
            h5f["blcorrf"] = self.blcorrf.__name__
        else:
            h5f["blcorrf"] = 'None'

    def load(self, name=None, light_version=False):
        if name is None:
            name = os.path.splitext(self.path)[0] + ".h5"

        h5f = h5py.File(name, "r")
        
        for j in range(len(self.dic_names)):
            dname = self.dic_names[j]
            out_dic = getattr(self,dname)
            if dname not in h5f.keys():
                if legacy_dic_names[dname] in h5f.keys():
                    dname=legacy_dic_names[dname]
                else:
                    print("Warning: key {} not in loaded file".format(dname))
                    continue
            if light_version and dname=="traces_dict":
                continue
            ds = h5f[dname]
            for key in ds.keys():
                dd = ds[key][()]
                try:
                    out_dic[int(key)] = dd
                except:
                    out_dic[key] = dd
        for par in h5f["parameters"].keys():
            setattr(self, par, h5f["parameters"][par][()])

    def registration(self, nreg, plot=False):
        self.stack, shifts = stackreg(self.stack, nreg, plot=plot)
        self.nreg = nreg
        self.shifts = shifts

    def set_threshold_map(self, th_map):
        self.threshold_map = th_map

    def set_bleaching_function(self, blcorrf, wsize=None):
        if wsize is not None:
            self.blcorrf = lambda x: blcorrf(x, wsize=wsize)
        else:
            self.blcorrf = blcorrf
        self.bl_kernel_size = str(wsize)

    def correlate_stack(self, nSum):
        """Only method that correlates """
        if nSum > self.stack.shape[1] or nSum > self.stack.shape[2]:
            raise ValueError("Trying to sum more pixels than present in the stack")
        if nSum not in self.fcs_curves_dict.keys():
            print("correlating stack, binning {}".format(nSum))
            correls = []
            traces = []
            u, v, w = self.stack.shape
            thumbnail = np.zeros((v//nSum, w//nSum))
            for i in range(v//nSum):
                ctmp = []
                trtmp = []
                for j in range(w//nSum):
                    trace = self.stack[:, i * nSum:i * nSum + nSum,
                            j * nSum:j * nSum + nSum].mean(axis=(1, 2))
                    thumbnail[i,j]=trace.mean()
                    if self.blcorrf is not None:
                        trace = self.blcorrf(trace)

                    if self.clipval > 0:
                        trace = trace[self.clipval:-self.clipval]

                    corr = multipletau.autocorrelate(trace, normalize=True, 
                                                     deltat=self.dt,m=16)[1:]
                    ctmp.append(corr)
                    trtmp.append(trace)
                correls.append(ctmp)
                traces.append(trtmp)
            correls = np.asarray(correls)
            
            self.fcs_curves_dict[nSum] = correls
            self.traces_dict[nSum] = np.asarray(traces)
            self.thumbnails_dict[nSum] = thumbnail

    def get_curve(self, i0=0, j0=0, nSum=1):
        self.correlate_stack(nSum)

        correls = self.fcs_curves_dict[nSum]
        correl = correls[i0, j0]
        return correl

    def get_all_curves(self, nSum=1, spacing=0, npts=None, plot=True):
        self.correlate_stack(nSum)
        correls = self.fcs_curves_dict[nSum]

        if self.threshold_map is not None:
            pass

        u, v = correls.shape[0], correls.shape[1]
        spc = np.percentile(correls[:, :, :, 1], 95) * spacing
        if plot:
            plt.figure()
            for j in range(u):
                for k in range(v):
                    corr = correls[j, k]
                    plt.semilogx(corr[:, 0], corr[:, 1] + (j * u + k) * spc)
        return correls

    def get_correlation_dict(self):
        return self.fcs_curves_dict

    def average_curve(self, nSum=1, plot=False, chi_th = None, 
                      ith = None, normalise_with_N = False):
        self.correlate_stack(nSum)
        # calculates mask
        th_map = None
        th_map = np.ones_like(self.thumbnails_dict[nSum]).astype(bool)
        if chi_th is not None:
            if nSum not in self.chisquares_dict.keys():
                self.calculate_chisquares()
            chimap = self.chisquares_dict[nSum]
            th_map = np.logical_and(th_map, chimap<chi_th)
        if ith is not None:
            thr = intensity_threshold(ith,self.thumbnails_dict[nSum])
            th_map = np.logical_and(th_map,self.thumbnails_dict[nSum]>thr)
        self.set_threshold_map(th_map)
        
        correls = self.fcs_curves_dict[nSum].copy()
        if normalise_with_N:
            correls[:, :, :, 1] = correls[:, :, :, 1]*self.fit_results_dict[nSum][:,:,0][:,:,np.newaxis]
        if self.threshold_map is not None:
            cs = correls[self.threshold_map]
            avg = cs[:, :, 1].mean(axis=(0))
        else:
            avg = correls[:, :, :, 1].mean(axis=(0, 1))

        if plot:
            plt.figure()
            plt.semilogx(correls[0, 0, :, 0], avg)
            plt.axhline(0, linestyle="--", color="k")
            plt.title("Average of binning {}".format(nSum))
        return np.array([correls[0, 0, :, 0], avg]).T

    def trace(self, plot=False):
        tr = self.stack.sum(axis=(1, 2))
        if plot:
            xtr = np.arange(tr.size) * self.dt
            plt.figure()
            plt.plot(xtr, tr)
            plt.xlabel("time")
            plt.ylabel("Counts")
        return tr

    def plot_sum_img(self):
        plt.figure()
        plt.imshow(self.stack.mean(axis=0))

    def binned_average_curves(self, sum_list, plot=True, n_norm=8):
        if plot:
            fig, axes = plt.subplots(1, 2)
        all_corrs = []
        for sl in sum_list:
            corr = self.average_curve(nSum=sl)
            all_corrs.append(corr)
            if plot:
                axes[0].semilogx(corr[:, 0], corr[:, 1],
                                 label="Binning {}".format(sl))
                axes[1].semilogx(corr[:, 0], corr[:, 1] / corr[:n_norm, 1].mean(),
                                 label="Binning {}".format(sl))
        if plot:
            axes[1].axvline(self.dt * self.stack.shape[0] / 100, color="k",
                            label="Max unbiased transit time")
            axes[1].axvline(self.dt * 10, color="k",
                            label="Min unbiased transit time")
            axes[1].legend()
            axes[0].set_title("Raw curves")
            axes[1].set_title("Normalised curves")
            axes[0].set_xlabel(r"$\rm \tau (A.U)$")
            axes[0].set_ylabel(r"$\rm G(\tau)$")
            axes[1].set_xlabel(r"$\rm \tau (A.U)$")
            axes[1].set_ylabel(r"$\rm G(\tau)$")
            axes[0].axhline(0, color="k")
            axes[1].axhline(0, color="k")
        return all_corrs

    def fit_curves(self, fitter, xmax=None, xmin=None):
        self.fitter = fitter
        nsums = self.fcs_curves_dict.keys()
        for nsum in nsums:
            self.fitter.set_sum(nsum)
            # updates corels size so that everyone has the same size
            if xmax is not None:
                c0 = self.fcs_curves_dict[nsum][0,0,:,0]
                index0 = np.max(np.where((c0-xmax)<0))
                self.fcs_curves_dict[nsum] = self.fcs_curves_dict[nsum][:,:,:index0,:]
            if xmin is not None:
                c0 = self.fcs_curves_dict[nsum][0,0,:,0]
                index0 = np.min(np.where((c0-xmin)>0))
                self.fcs_curves_dict[nsum] = self.fcs_curves_dict[nsum][:,:,index0:,:]
            correls = self.fcs_curves_dict[nsum]
            popts = []
            yhs = []
            chisquares = []
            
            for j in range(correls.shape[0]):
                popt_tmp = []
                yh_tmp = []
                chisquares_tmp = []
                for k in range(correls.shape[1]):
                    corr = correls[j, k]
                    popt, yh = fitter.fit(corr)
                    chi = new_chi_square(corr[:, 1], yh)
                    popt_tmp.append(popt)
                    yh_tmp.append(yh)
                    chisquares_tmp.append(chi)

                popts.append(popt_tmp)
                yhs.append(yh_tmp)
                chisquares.append(chisquares_tmp)

            self.fit_results_dict[nsum] = np.array(popts)
            self.yh_dict[nsum] = np.array(yhs)
            self.chisquares_dict[nsum] = np.array(chisquares)
        self.fitting_parameters_dict = fitter.parameters_dict
        
    def calculate_chisquares(self):
        nsums = self.fcs_curves_dict.keys()
        for nsum in nsums:
            correls = self.fcs_curves_dict[nsum]
            fits = self.yh_dict[nsum]
            chisquares = []
            for j in range(correls.shape[0]):
                chisquares_tmp = []
                for k in range(correls.shape[1]):
                    corr = correls[j, k]
                    yh = fits[j, k]
                    chi = new_chi_square(corr[:, 1], yh)
                    chisquares_tmp.append(chi)
                chisquares.append(chisquares_tmp)
            self.chisquares_dict[nsum] = np.array(chisquares)

    def calculate_square_error(self):
        """Calculates the square error to the fit"""
        nsums = self.fcs_curves_dict.keys()
        for nsum in nsums:
            correls = self.fcs_curves_dict[nsum]
            fits = self.yh_dict[nsum]
            square_err = []
            for j in range(correls.shape[0]):
                square_err_tmp = []
                for k in range(correls.shape[1]):
                    corr = correls[j, k]
                    yh = fits[j, k]
                    
                    if yh[0]==0:
                        err=np.inf
                    else:
                        err = np.mean(( (corr[:, 1]- yh)/yh[0])**2)
                    square_err_tmp.append(err)
                square_err.append(square_err_tmp)
            self.square_err_dict[nsum] = np.array(square_err)
            
    def parameter_map(self, nsum=None, parn=1):
        print('Caution! You are using a resampled parameter map')
        if nsum is None:
            nsum = min(self.fit_results_dict.keys())
        out = self.fit_results_dict[nsum][:, :, parn]
        out = np.repeat(out, nsum, axis=0)
        out = np.repeat(out, nsum, axis=1)
        return out

    def get_param_threshold(self, nsum, thf=None, plot=False, parn=1):

        img = self.stack.mean(axis=0).astype(float)
        if thf is None:
            thf = threshold_otsu
        thresholded = img > thf(img)

        uu = thresholded.shape[0] - thresholded.shape[0] % nsum
        vv = thresholded.shape[1] - thresholded.shape[1] % nsum

        to_keep = self.get_threshold_map(nsum, thf=thf)
        if plot:
            plt.figure()
            plt.subplot(121)
            plt.imshow(img)

            plt.subplot(122)
            plt.imshow(to_keep)

        return self.fit_results_dict[nsum][:uu, :vv, parn][to_keep]

    def get_param_coord(self, nsum, i0, j0, parn=1, exclude_neg=True):
        """Get the value of given parameters for all binning values below nsum"""
        sums = self.fcs_curves_dict.keys()
        sums = sorted([w for w in sums if w <= nsum])
        ds_means = list()
        ds_std = list()
        for ns in sums:
            i00 = int(np.ceil(i0 * nsum / ns))
            i01 = int(np.floor((i0 + 1) * nsum / ns))

            j00 = int(np.ceil(j0 * nsum / ns))
            j01 = int(np.floor((j0 + 1) * nsum / ns))

            ds = self.fit_results_dict[ns][i00:i01, j00:j01, parn]
            if exclude_neg:
                ds = ds[ds >= 0]
            ds_means.append(np.mean(ds))
            ds_std.append(np.std(ds))
        sums = np.asarray(sums)
        ds_means = np.asarray(ds_means)
        ds_std = np.asarray(ds_std)
        mask = ~np.isnan(ds_means)
        return sums[mask], ds_means[mask], ds_std[mask]

    def describe(self):
        """Returns a dictionary of parameters describing the stack. Includes 
        physical parameters like pixel size and fitting parameters"""
        out_dict=self.fitting_parameters_dict.copy()
        out_dict['psize [nm]'] = self.xscale*10**3
        out_dict['dt [ms]'] = self.dt*10**3
        path,fname = os.path.split(self.path)
        def sdec(f): 
            if type(f)==bytes: return f.decode('utf-8') 
            else: return f
        out_dict["path"] = sdec(path)
        out_dict["filename"] = sdec(fname)
        
        return out_dict
    
    def get_acf_coord(self, nsum, i0, j0, parn=1, average=True):

        sums = self.fcs_curves_dict.keys()
        sums = sorted([w for w in sums if w <= nsum])
        all_corrs = list()
        all_yhs = list()
        all_ns = list()
        for ns in sums:
            i00 = int(np.ceil(i0 * nsum / ns))
            i01 = int(np.floor((i0 + 1) * nsum / ns))

            j00 = int(np.ceil(j0 * nsum / ns))
            j01 = int(np.floor((j0 + 1) * nsum / ns))
            corrs = self.fcs_curves_dict[ns][i00:i01, j00:j01]
            yhs = self.yh_dict[ns][i00:i01, j00:j01]
            if average:
                corrs = corrs.mean(axis=(0, 1))
                yhs = yhs.mean(axis=(0, 1))
            if not np.isnan(corrs).all():
                all_corrs.append(corrs)
                all_yhs.append(yhs)
                all_ns.append(ns)
        sums = np.asarray(sums)
        all_corrs = all_corrs
        all_yhs = all_yhs
        return all_ns, all_corrs, all_yhs

    def downsample_image(self, nsum):
        u, v, w = self.stack.shape

        img = self.stack.mean(axis=0)
        out = np.zeros((v // nsum, w // nsum))
        for i in range(v // nsum):
            for j in range(w // nsum):
                px = img[i * nsum:i * nsum + nsum, j * nsum:j * nsum + nsum].mean()
                out[i, j] = px
        return out

    def downsample_time(self, ndown):
        """Downsamples an image stack in time"""
        nframes = (self.stack.shape[0]//ndown)*ndown
        self.stack = self.stack[:nframes]
        self.stack = self.stack.reshape((ndown,nframes//ndown,self.stack.shape[1],
                            self.stack.shape[2]),order="F").mean(axis=0)
        self.dt = self.dt*ndown
    
    def distance2mask(self, index, nindices = None):
        """Calculates distances of every point in each mask to the centroid of a 
        given mask"""
        nsums = self.fcs_curves_dict.keys()
        if nindices is None:
            raise KeyError("Please indicate what is the number of classes in your masks")
            
        self.distance_to_masks_dict = dict(zip(nsums,[[] for w in nsums]))
        self.mask_index = index
        
        for nsum in nsums:
            maxind = self.mask.max()
            ncells = int(np.ceil(maxind/nindices))
            # indices = sorted(list(filter(lambda x: x!=0, np.unique(self.mask))))
            out_map = np.zeros_like(self.thumbnails_dict[nsum])
            for ncell in range(ncells):
                indices = np.arange(1,nindices+1) + ncell*nindices
                msk = downsample_mask(self.mask, nsum)
                index_ref = index + ncell*nindices
                centroid = np.array([np.mean(x) for x in np.where(msk==index_ref)])
                for ind in indices:
                    coords = np.array(np.where(ind==msk))
                    for j in range(coords.shape[1]):
                        vec = coords[:,j]
                        dist = np.linalg.norm(vec-centroid)
                        out_map[*vec] = dist
            self.distance_to_masks_dict[nsum] = out_map

    def extract_results(self, ith = None, 
                      chi_threshold = None, use_mask=True):
        """Extracts results like diffusion coefficient, chisquares etc. Meant to replace
        io.get_fit_error. Returns a single dictionary, each key being a quantity
        measured"""
        nsums = self.fcs_curves_dict.keys()
        
        mk_outdic= lambda:dict(zip(nsums,[[] for w in nsums]))
        # results: a dictionary containing results dictionaries
        
        results = {"non_linear_chis":mk_outdic(),
                   "indices":mk_outdic(),
                   "square_errors":mk_outdic(),
                   "intensities":mk_outdic(),
                   "valid_fraction":mk_outdic(), # fraction of curves removed from chi
                   }
        
        if len(self.distance_to_masks_dict )>0:
            results["distance_to_mask_{}".format(self.mask_index)] = mk_outdic()
            
        # populates results dict with every interesting fitting parameter
        if self.fitter is None:
            self.fitter = Fitter(self.fitting_parameters_dict)
        parameter_names_dict = self.fitter.get_parameter_names_dict()
        
        for k in parameter_names_dict:
            results[k] = mk_outdic()
            
        self.calculate_chisquares()
        self.calculate_square_error()
        
        if self.mask is None:
            use_mask=False
            
        for jj, nsum in enumerate(nsums):
            curves = self.fcs_curves_dict[nsum]
            curves_fits = self.yh_dict[nsum]
            intensities = self.thumbnails_dict[nsum].reshape(-1)
            square_errors = self.square_err_dict[nsum].reshape(-1)
            
            ycurves = curves[:,:,:,1]
            
            chis = np.sqrt(((ycurves-curves_fits)**2).sum(axis=2) )
            
            chis_new = np.zeros_like(chis)
            
            for i in range(chis_new.shape[0]):
                for j in range(chis_new.shape[1]):
                    
                    chis_new[i,j] = new_chi_square(ycurves[i,j], curves_fits[i,j])
        
            curves_reshaped = curves.reshape((curves.shape[0]*curves.shape[1],curves.shape[2],2))
            fits_reshaped = curves_fits.reshape((curves.shape[0]*curves.shape[1],curves.shape[2]))
            
            msk = np.ones_like(intensities, dtype = bool)
            
            msk[(self.fit_results_dict[nsum]==-1).all(axis=-1).reshape(-1)] = False
            indices=np.zeros_like(intensities)
            # set mask for measurements. msk is boolean
                
            if chi_threshold is not None:
                msk = np.logical_and(msk, chis_new.reshape(-1)<chi_threshold)
            # mask of values we want to keep, regardless of validity
            msk_for_values =np.ones_like(intensities,dtype=bool)
            if ith is not None and not use_mask:
                ithr = self.intensity_threshold_f(ith,intensities)
                msk = np.logical_and(msk,
                                     intensities>ithr)
                indices[msk]=1
                msk_for_values = intensities>ithr
            if use_mask and self.mask is not None:
                # !!! add selection of hard mask
                indices = indices_intensity_threshold(
                    downsample_mask(self.mask, nsum).reshape(-1),ith,intensities)
                msk = np.logical_and(msk,indices>0)
                msk_for_values = indices>0
            if np.count_nonzero(msk_for_values)==0:
                valid_measures=0
            else:
                valid_measures = np.count_nonzero(msk)/np.count_nonzero(msk_for_values)
            chis_new = chis_new.reshape(-1)[msk]
            curves_reshaped = curves_reshaped[msk]
            fits_reshaped = fits_reshaped[msk]
            
            indices = indices.reshape(-1)[msk]
            square_errors = square_errors.reshape(-1)[msk]
            intensities = intensities.reshape(-1)[msk]
            # For loop populates
            for k in parameter_names_dict.keys():
                ind = parameter_names_dict[k]
                val = self.fit_results_dict[nsum][:,:,ind]
                results[k][nsum] = val.reshape(-1)[msk]
            results["non_linear_chis"][nsum] = chis_new
            results["indices"][nsum] = indices
            results["square_errors"][nsum] = square_errors
            results["intensities"][nsum] = intensities
            results["valid_fraction"][nsum] = valid_measures
            if len(self.distance_to_masks_dict )>0:
                distances_to_masks =  self.distance_to_masks_dict[nsum].reshape(-1)
                results["distance_to_mask_{}".format(self.mask_index)][nsum] = distances_to_masks[msk]
        return results
    
    def plot_parameter_maps(self, nsums, parn=1, cmap="jet", vmin=None,
                            vmax=None, maxval=None):
        assert len(nsums) >= 1
        assert len(nsums) <= 5

        nr = 2
        nc = len(nsums)

        fig, axes = plt.subplots(nr, nc, sharex="col", sharey="col")

        parmaps = list()
        for j in range(len(nsums)):
            nsum = nsums[j]
            parmap = self.fit_results_dict[nsum][:, :, parn].copy()
            if maxval is not None:
                parmap[parmap > maxval] = np.nan
            im = self.downsample_image(nsum)
            parmaps.append(parmap)
            ax0 = axes[0, j]
            ax1 = axes[1, j]
            extent = [0, im.shape[1] * nsum * self.xscale, 0,
                      im.shape[0] * nsum * self.xscale]
            im0 = ax0.imshow(im, cmap="gray", extent=extent)
            ax0.set_title('Binning {}'.format(nsum))
            im1 = ax1.imshow(parmap, cmap=cmap, vmin=vmin, vmax=vmax,
                             extent=extent)
            ax1.set_xlabel("x [µm]")
            ax1.set_ylabel("y [µm]")

            fig.colorbar(im1, ax=ax1)
            fig.colorbar(im0, ax=ax0)

    def plot_amplitudes(self, sum_list):
        averages = self.binned_average_curves(sum_list, plot=False)
        nvals = []
        for j in range(len(sum_list)):
            nn = averages[j][:8, 1].mean()
            nvals.append(nn)
        plt.figure()
        plt.plot(sum_list, np.sqrt(nvals), label="measured")
        plt.plot(sum_list, np.sqrt(nvals)[0] * sum_list[0] / np.array(sum_list),
                 color="k", linestyle="--", label="theory")
        plt.xlabel("Binning")
        plt.ylabel("Sqrt Curve amplitude")
        plt.xscale('log')
        plt.legend()

    def plot_curve(self, i0=0, j0=0, nSum=1):
        correl = self.get_curve(i0=i0, j0=j0, nSum=nSum)

        plt.figure()
        plt.semilogx(correl[:, 0], correl[:, 1])

    def plot_D(self, show_acceptable=True):
        nsums = sorted(self.fit_results_dict.keys())
        nsums = np.asarray(nsums)
        ds_means = list()
        ds_std = list()

        if show_acceptable:
            psize = self.fitter.parameters_dict["a"]
            if "sigmaxy" in self.fitter.parameters_dict:
                sigmaxy = self.fitter.parameters_dict["sigmaxy"]
            else:
                sigmaxy = self.fitter.parameters_dict["sigma"]
            sigmaxy = sigmaxy * np.sqrt(8 * np.log(2))  # fwhm
            observation_sizes = np.sqrt((nsums * psize) ** 2 + sigmaxy ** 2)
            if self.fitter.name == "2D":
                factor = 4
            elif self.fitter.name == "3D":
                factor = 6
            dmax = observation_sizes ** 2 / (factor * self.dt * 10)
            dmins = observation_sizes ** 2 / (factor * self.dt * self.stack.shape[0] / 100)
        for ns in nsums:
            ds = self.fit_results_dict[ns][:, :, 1]
            # !!! check threshold maps here
            ds_means.append(np.median(ds))
            ds_std.append((np.percentile(ds, 75) - np.percentile(ds, 25)) / 2)
        ds_means = np.asarray(ds_means)
        ds_std = np.asarray(ds_std)

        plt.figure()
        if show_acceptable:
            plt.plot(nsums, dmins, color="gray")
            plt.plot(nsums, dmax, color="gray")
        plt.errorbar(nsums, ds_means, yerr=ds_std, capsize=5)
        plt.xlabel("Binning size")
        plt.ylabel("D (um2/s)")
        ymin = np.min(ds_means - ds_std)
        ymax = np.max(ds_means + ds_std)
        plt.ylim(bottom=ymin, top=ymax)
        return nsums, ds_means, ds_std

    def plot_intensity_correlation(self, nsum, parn=1):
        parmap = self.fit_results_dict[nsum][:, :, parn]
        parameters = list()
        intensities = list()
        stack_avg = self.stack.mean(axis=0)
        for j in range(parmap.shape[0]):
            for k in range(parmap.shape[0]):
                parameters.append(parmap[j, k])
                intensities.append(stack_avg[j * nsum:(j + 1) * nsum,
                                   k * nsum:(k + 1) * nsum].mean())
        plt.figure()
        plt.scatter(parameters, intensities)
        plt.xlabel("Parameter")
        plt.ylabel("Intensity")
        return np.asarray(parameters), np.asarray(intensities)

    def plot_fits(self, nSum, maxcurves=None, dz=0.2):
        curves = self.fcs_curves_dict[nSum]
        chisquares = self.chisquares_dict[nSum]
        sp1 = np.asarray(curves.shape)
        fits = self.yh_dict[nSum]
        fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
        jj = 0

        indices = np.arange(sp1[0] * sp1[1])
        indices1 = indices // sp1[1]
        indices2 = indices - (indices1) * sp1[1]
        if maxcurves is not None:
            indices = np.random.choice(np.arange(sp1[0] * sp1[1]), maxcurves)

        for i in indices:
            j = indices1[i]
            k = indices2[i]
            corr = curves[j, k]
            chi = chisquares[j, k]
            yh = fits[j, k]
            a = corr[:3, 1].mean()
            axes[0].semilogx(corr[:, 0], corr[:, 1] / a + dz * jj)
            axes[0].semilogx(corr[:yh.size, 0], yh / a + dz * jj, color="k", linestyle="--")
            axes[1].semilogx(corr[:yh.size, 0], yh / a - corr[:yh.size, 1] / a + dz * jj,
                             label="curve ({},{}), chisquare {}".format(j, k, chi))
            jj += 1
        axes[1].legend()

    def plot_fits_ordered(self, nSum, maxcurves=None, dz=-0.2, order_dict=None):
        """Plots FCS curves ordered following values contained in a specific 
        dictionary"""
        assert order_dict is not None

        curves = self.fcs_curves_dict[nSum].reshape(-1, *self.fcs_curves_dict[nSum].shape[2:])
        parameters = order_dict[nSum].reshape(-1)
        fits = self.yh_dict[nSum].reshape(-1, *self.yh_dict[nSum].shape[2:])
        fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
        jj = 0

        indices = np.arange(curves.shape[0])
        if maxcurves is not None and maxcurves < curves.shape[0]:
            indices = np.random.choice(curves.shape[0], maxcurves)

        indices = [x for _, x in sorted(zip(parameters[indices], indices))]

        for i in indices:
            corr = curves[i]
            param = parameters[i]
            yh = fits[i]
            a = yh[0]
            axes[0].semilogx(corr[:, 0], corr[:, 1] / a + dz * jj)
            axes[0].semilogx(corr[:yh.size, 0], yh / a + dz * jj, color="k", linestyle="--")
            axes[1].semilogx(corr[:yh.size, 0], yh / a - corr[:yh.size, 1] / a + dz * jj,
                             label="corder parameter {}".format(param))
            jj += 1
        axes[1].legend()

    def plot_random_intensity(self, nSum=None):
        if nSum is None:
            nSum = min(self.traces_dict.keys())
        traces_arr = self.traces_dict[nSum]
        trace_raw = self.stack.mean(axis=(1, 2))
        u, v = traces_arr.shape[:2]
        u1 = np.random.choice(u)
        v1 = np.random.choice(v)
        trace = traces_arr[u1, v1]
        trace_raw = self.stack[:, u1 * nSum:(u1 + 1) * nSum, v1 * nSum:(v1 + 1) * nSum].mean(axis=(1, 2))

        plt.figure()
        plt.plot(trace, label="Corrected")
        plt.axhline(trace.mean(), color="k", linestyle="--")
        plt.plot(trace_raw, label="Raw average intensity")
        plt.xlabel("Time (frames)")
        plt.ylabel("Intensity")
        plt.legend()

    def plot_taus(self, show_acceptable=True):
        nsums = sorted(self.fit_results_dict.keys())
        nsums = np.asarray(nsums)
        ds_means = list()
        ds_std = list()

        if show_acceptable:
            psize = self.fitter.parameters_dict["a"]
            if "sigmaxy" in self.fitter.parameters_dict:
                sigmaxy = self.fitter.parameters_dict["sigmaxy"]
            else:
                sigmaxy = self.fitter.parameters_dict["sigma"]
            sigmaxy = sigmaxy * np.sqrt(8 * np.log(2))  # fwhm
            observation_sizes = np.sqrt((nsums * psize) ** 2 + sigmaxy ** 2)

        nsums = list(nsums)
        for j, ns in enumerate(nsums):
            ds = self.fit_results_dict[ns][:, :, 1]
            #TODO here
            """if self.threshold_map is not None:
                th = self.get_threshold_map(ns)
                ds = ds[th]"""
            if len(ds) == 0:
                nsums.pop(j)
                continue
            taus = observation_sizes[j] ** 2 / ds
            ds_means.append(np.median(taus))
            ds_std.append((np.percentile(taus, 75) - np.percentile(taus, 25)))
        ds_means = np.asarray(ds_means)

        plt.figure()
        if show_acceptable:
            dmax = np.ones_like(nsums) * (self.dt * 10)
            dmins = np.ones_like(nsums) * self.dt * self.stack.shape[0] / 100
            plt.plot(nsums, dmins, color="gray")
            plt.plot(nsums, dmax, color="gray")
        plt.errorbar(nsums, ds_means, yerr=ds_std, capsize=5)
        plt.xlabel("Binning size")
        plt.ylabel("tau (s)")
        return nsums, ds_means, ds_std
