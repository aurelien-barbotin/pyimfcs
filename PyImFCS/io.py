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
