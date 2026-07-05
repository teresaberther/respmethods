#!/usr/bin/env python3

import numpy as np
import os
import copy
import time

from scipy.stats import zscore
from scipy.signal import decimate
from multiprocessing import Pool

# custom
from respymethods.GenerateSurrogates import iaaft
from respymethods.PhaseExtraction import two_point_interp
from respymethods.Helpers import plv, prune_nan


def load_data(fname):
    """
    convenience function to be changed according to the data format 
    it should return a 1 dimensional array 
    """
    
    x = np.load(fname)
    return x


def parallel_process(dict_in):
    """
    dict_in = {
            'x_norm',
            'phase_vect',
            'n_iter_process',
            }
    """
    x_norm = dict_in['x_norm']
    phase_vect = dict_in['phase_vect']
    n_iter_process = dict_in['n_iter_process']
    crit_plv = dict_in['crit_plv']
    
    surr_data = np.zeros((len(x_norm), n_iter_process))
    
    i=0;
    while i<n_iter_process:
        
        surr_ = iaaft(x_norm)
        theta_ = two_point_interp(surr_)
        tmp1, tmp2 = prune_nan(phase_vect, theta_)
        this_plv = plv(tmp1, tmp2)
        
        if this_plv<crit_plv:
            surr_data[:, i] = surr_
            i+=1         
    
    return surr_data


def optimize_cpu_dist(n_iter, maxusage=.75):
    """
    optimize the number of cpus based on the number of permutations to be 
    performed.
    """
    n_cpus_tot = os.cpu_count()
    n_cpus_avail = np.floor(n_cpus_tot*maxusage)

    mod = n_iter%n_cpus_avail    
    while mod:
        n_cpus_avail -= 1
        mod = n_iter%n_cpus_avail

    return int(n_cpus_avail)
        
    

if __name__ == '__main__':
    """
    The script should take as input:
        path-to-input-file
        path-to-output-file
        n iterations
        n cpus (optional)
    """
    
    onset_time = time.perf_counter()

    in_file = "data/sampledata.npy"
    out_file = "data/outdata.npy"
    n_iter = 1000
        
    # other configs that might be put as input in some future versions
    # for the time being brutally defaulted here
    downsmpl_rate = 10 # if >13, not recomendations from sp.signal.decimate doc
    fs = 1000
    win_len = int(.4*fs)
    crit_z = 1 # !!! still unused !!!. Does it make sense to have such limit only for peaks?
    crit_plv = .1
    
    # load our signal
    x = load_data(in_file)
    
    # overhead operations 
    
    ################### 1. Preprocessing
    
    # 1.1 Smoothing
    x_smooth = np.convolve(x, np.ones(win_len), 'same')/win_len
    # 1.2 Downsampling
    x_down = decimate(x_smooth, downsmpl_rate)
    # 1.3 Z scoring
    x_norm = zscore(x_down)

    ################### 2. Phase computation
    
    phase_vect = two_point_interp(x_norm)


    ################### 3. Prepare for multiprocessing
    pool_size = optimize_cpu_dist(n_iter)
    
    input_processes = []
    for i in range(pool_size):
        
        tmp_dict = {
            'x_norm' : copy.deepcopy(x_norm),
            'phase_vect' : copy.deepcopy(phase_vect),
            'n_iter_process' : int(n_iter / pool_size),
            'crit_plv' : crit_plv
            }
        
        input_processes.append(tmp_dict)
    
    with Pool(pool_size) as p:
        result = p.map(parallel_process, input_processes)
        
    out_array = np.concat(result, axis=-1)
    np.save(out_file, out_array)

    offset_time = time.perf_counter()
    elapsed_time = np.round(offset_time-onset_time, 3)
    
    print("\nDone\n")
    print(f"{n_iter} iterations: {elapsed_time} s")
    
