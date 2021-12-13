import h5py

import numpy as np
import pandas as pd

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
        out_dic[int(key)] = dd
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


def save_as_excel(out_name,files,nsums,all_diffs,all_chis):
    """Saves the results of several FCS experiments in a single excel file.
    Parameters:
        out_name (string): output file name. Extension is added later
        file (list): list of all filenames in the dataset
        nsums (list): integers, different pixel binnings
        all_diffs (list): list of diffusion coefficients dictionaries, same 
            size as files. Keys are the binning values (typically 2,3,4 or 2,4,8)"""
    extension = ".xlsx"
    if not out_name.endswith(extension):
        out_name+=extension
    
    all_dfs = {}
    for nsum in nsums:
        df0 = []
        
        for j in range(len(files)):
            fname = files[j]
            diff = all_diffs[j][nsum]
            chis = all_chis[j][nsum]
            repeats_arr = np.full(diff.size, j)
            name_arr = np.full(diff.size, fname)
            nsum_arr = np.full(diff.size, nsum)
            out_arr = np.array([name_arr,repeats_arr, diff, nsum_arr, chis]).T
            df = pd.DataFrame(out_arr, columns = ["name", "repeat","D [µm²/s]","binning","fit error"])
            df0.append(df)
        all_dfs[nsum] = pd.concat(df0)
        
    with pd.ExcelWriter(out_name) as writer:  
        dfs_total = []
        for nsum in nsums:
            dfpooled = all_dfs[nsum]
            dfs_total.append(dfpooled)
            dfpooled.to_excel(writer, sheet_name = "nsum {}".format(nsum))
        dfs_total = pd.concat(dfs_total)
        dfs_total.to_excel(writer, sheet_name = "all")
        
        parameters_dict = {"chi_threshold":new_chis_threshold,
                           "intensity_threshold": intensity_threshold}
        df_pars = pd.DataFrame(parameters_dict, index=[0]).T
        df_pars.to_excel(writer, sheet_name = "parameters")
