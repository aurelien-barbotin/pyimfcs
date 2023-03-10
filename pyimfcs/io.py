import h5py

import numpy as np
import pandas as pd
import tifffile

from pyimfcs.metrics import new_chi_square, intensity_threshold
from pyimfcs.methods import downsample_mask


# parameters: [N,D,offset]
def get_dict(h5f, dname):
    """Extract a dictionary stored in a h5 file under the keyword
    dname"""
    if dname not in h5f.keys():
        print("grougrou")
        raise KeyError('Parameter not in file')
    out_dic = {}
    ds = h5f[dname]
    for key in ds.keys():
        dd = ds[key][()]
        try:
            out_dic[int(key)] = dd
        except:
            out_dic[key] = dd
    return out_dic

def get_diffusion_df(name, repeat, condition, nsum):
    """Extracts diffusion coefficients from a h5 file and formats it into a pandas
    Dataframe for future extraction"""
    h5f = h5py.File(name, "r")
    
    parameters_fits = get_dict(h5f, "parameters_fits")
    
    diff = parameters_fits[nsum][:,:,1].reshape(-1)
    
    condition_arr = np.full(diff.size, condition)
    repeats_arr = np.full(diff.size, repeat)
    name_arr = np.full(diff.size, name)
    
    out_arr = np.array([diff, condition_arr, repeats_arr, name_arr]).T
    
    df = pd.DataFrame(out_arr, columns = ["D [µm2/s]", "condition", "repeat", "name"])
    h5f.close()
    
    return df

def summarise_df(df):
    """Merges results of a Dataframe containing multiple experiments"""
    means = df.groupby('repeat')['D [µm²/s]'].mean()
    std = df.groupby('repeat')['D [µm²/s]'].std()
    medians = df.groupby('repeat')['D [µm²/s]'].median()
    count = df.groupby('repeat')['D [µm²/s]'].count()
    repeat = df.groupby('repeat')['repeat'].last()
    filename = df.groupby('repeat')['filename'].last()
    binning = df.groupby('repeat')['binning'].last()

    out_df = pd.DataFrame.from_dict({
                        "filename": filename,
                        "binning": binning,
                        "repeat":repeat,
                        "Mean":means,
                        "Stdev": std,
                        "Median":medians,
                        "Count":count})
    return out_df

def save_as_excel(out_name,files,nsums,all_diffs,all_chis, all_ns, parameters_dict= {}):
    """Saves the results of several FCS experiments in a single excel file.
    Parameters:
        out_name (string): output file name. Extension is added later
        files (list): list of all filenames in the dataset
        nsums (list): integers, different pixel binnings
        all_diffs (list): list of diffusion coefficients dictionaries, same 
            size as files. Keys are the binning values (typically 2,3,4 or 2,4,8)
        all_chis (list): same format as all_diffs, containing fit error metric
        all_ns (list): same format as all_diffs, containing numbers of molecules
        parameters_dict (dict): dictionary of misc processing parameters like 
            chi_threshold or intensity_threshold
        
        TODO: refactor all list of dicts into a single parameter"""

    if out_name[-5:]==".xlsx":
        out_name=out_name[:-5]
    
    all_dfs = {}
    for nsum in nsums:
        df0 = []
        
        for j in range(len(files)):
            fname = files[j]
            diff = all_diffs[j][nsum]
            chis = all_chis[j][nsum]
            nmols = all_ns[j][nsum]
            
            repeats_arr = np.full(diff.size, j)
            name_arr = np.full(diff.size, fname)
            nsum_arr = np.full(diff.size, nsum)
            out_arr = np.array([name_arr,repeats_arr, diff, nsum_arr, chis, nmols]).T
            df = pd.DataFrame(out_arr, columns = 
                              ["filename", "repeat","D [µm²/s]","binning",
                               "fit error", "N"])
            df = df.astype({'filename':"str",
                           "repeat":"int",
                           "D [µm²/s]":"float",
                           "fit error":"float",
                           "N": "float"})
            df0.append(df)
        
        all_dfs[nsum] = pd.concat(df0)
        knm = "nsum {}".format(nsum)
        parameters_dict[knm+"_median"] = np.median(all_dfs[nsum]["D [µm²/s]"].values)
        
    with pd.ExcelWriter(out_name+".xlsx") as writer:  
        dfs_total = []
        for nsum in nsums:
            dfpooled = all_dfs[nsum]
            dfs_total.append(summarise_df(all_dfs[nsum]))
            dfpooled.to_excel(writer, sheet_name = "nsum {}".format(nsum))
        dfs_total = pd.concat(dfs_total)
        dfs_total.to_excel(writer, sheet_name = "summaries all")
        
        # parameters_dict = {"chi_threshold":new_chis_threshold,
        #                    "intensity_threshold": intensity_threshold}
        df_pars = pd.DataFrame(parameters_dict, index=[0]).T
        df_pars.to_excel(writer, sheet_name = "parameters")

def get_image_metadata(path):
    img = tifffile.TiffFile(path)
    meta_dict = img.imagej_metadata
    description = meta_dict.pop('Info')
    description = description.split('\n')
    for d in description:
        if len(d)>1 and '=' in d:
            oo = d.split('=')
            if len(oo)==2:
                k, val = oo
            elif len(oo)>2:
                k = oo[0]
                val = "=".join(oo[1:])
            k = k.strip(' ')
            val = val.strip(' ')
            try:
                meta_dict[k] = float(val)
            except:
                meta_dict[k] = val
    return meta_dict

def get_h5_metadata(file_h5):
    h5f = h5py.File(file_h5,'r')
    metadata = get_dict(h5f, "metadata")
    for k in metadata.keys():
        if type(metadata[k])==bytes:
            metadata[k] = metadata[k].decode('utf-8')
    return metadata

def get_metadata_bioformats(path):
    import czifile
    if path[-4:]!=".czi":
        raise ValueError('Incorrect file format')
    img = czifile.CziFile(path)
    metadata = img.metadata()
    
def save_tiff_withmetadata(outname,st, meta_dict):
    
    writer = tifffile.TiffWriter(outname,imagej=True)
    writer.write(st,metadata=meta_dict)
    
def merge_excels(fileslist_list, out_name, keep_single_indices = False, 
                 conditions = None, chi_threshold = None):
    """Merge ther esults of FCS experiments in a single file.
    Parameters:
        fileslist_list (list): list of lists. Each item of the main list contains 
            the excel files corresponding to one condition
        out_name (str): save name
        keep_single_indices (bool): if true, gives an indes to single acquisition. 
            If False, single index to an excel file.
        conditions (list): if specified, gives specific condition names 
            (e.g control or experiment) to acquisitions
        cho_threshold (float): if specified, eliminates all results with chi 
            above this value"""
    
    if len(fileslist_list)>1 and conditions is None:
        return ValueError('Please specify condition names')
    assert( conditions is None or len(conditions)==len(fileslist_list))
    
    if not out_name.endswith('.xlsx'):
        out_name+=".xlsx"
    
    # unravel
    files = [w for flist in fileslist_list for w in flist]
    excels = [pd.ExcelFile(w) for w in files]
    if conditions is not None:
        conditions_list = [conditions[j] for j, flist in enumerate(fileslist_list) for w in flist]

    
    all_names = [w.sheet_names for w in excels]
    names0 = all_names[0]
    kept_names = list()
    for name in names0:
          if np.all([name in sublist for sublist in all_names]):
              kept_names.append(name)
    
    all_dfs = {}
    for name in kept_names:
        dfs = []
        maxindex = 0
        for j, xl in enumerate(excels):
            df = xl.parse(sheet_name=name, index_col = 0)
            if name=="parameters":
                fname = files[j]
                df["file"] = fname
                
                if chi_threshold is not None:
                    cur_thr = df.loc["chi_threshold"]
                    if cur_thr.ndim>1:
                        cur_thr = cur_thr[0].values.max()
                    else:
                        cur_thr = cur_thr[0]
                    
                    if cur_thr<chi_threshold:
                        raise ValueError(
                    'Manually specified threshold is above already existing ones')
                    df['chi_threshold'] = chi_threshold
            else:
                if keep_single_indices:
                    df['repeat']+=maxindex
                    maxindex = df['repeat'].values.max()+1
                else:
                    df['repeat'] = maxindex
                    maxindex += 1
            if conditions is not None:
                df["condition"] = conditions_list[j]
            dfs.append(df)
                
        dfs = pd.concat(dfs)
        if chi_threshold is not None and name!="parameters":
            dfs = dfs[dfs['fit error']<chi_threshold]
        all_dfs[name] = dfs
    all_dfs["parameters"]
    with pd.ExcelWriter(out_name) as writer:  
            for name in kept_names:
                df_pars = all_dfs[name]
                df_pars.to_excel(writer, sheet_name = name)
                