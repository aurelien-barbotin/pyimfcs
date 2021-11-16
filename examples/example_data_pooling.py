import h5py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

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

# dicts_to_load  = ["correlations", "parameters_fits","yhat", "traces"]

files = glob.glob("/home/aurelien/Data/2021_11_03/SLB_POPC/*.h5")
files.sort()

condition = "4 ms"
repeat = 1
name = files[0]

nsum = 2


df = get_diffusion_df(name, repeat, condition, nsum)

files_4ms = [w for w in files if '1ms' not in w]
files_1ms = [w for w in files if '1ms' in w]

counters = np.array([0,0])
all_dfs = []

for j in range(len(files)):
    name = files[j]
    print('Processing file {}'.format(name))
    if "1ms" in name:
        condition = '1 ms'
        repeat = counters[0]
        counters[0]+=1
        df = get_diffusion_df(name, repeat, condition, nsum)
    else:
        condition = '4 ms'
        repeat = counters[1]
        counters[1]+=1
        df = get_diffusion_df(name, repeat, condition, nsum)
    all_dfs.append(df)

final_df = pd.concat(all_dfs)
final_df.to_csv("diffusion_results_pooling_{}.csv".format(nsum))
