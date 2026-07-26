from itertools import product

import numpy as np
from numpy.ma.core import diff
from scipy.stats import zscore
import os
import copy
import time
import matplotlib.pyplot as plt
from scipy.signal import decimate
from scipy.io import loadmat

# respymethods imports
import respymethods.RespStats as rmstats
import respymethods.Simulation as rmsim
from respymethods.GenerateSurrogates import iaaft
from respymethods.PhaseExtraction import two_point_interp
from respymethods.Helpers import plv, prune_nan
from respymethods.Simulation import create_phase_bins


def test_readtoml():
    assert len(rmsim.create_dict_from_toml("tests/sim_pars.toml")) == 11


def help_call_datagen():
    tmp_dict = rmsim.create_dict_from_toml("tests/sim_pars.toml")
    # generation of all the parameters combination
    combinations = product(*tmp_dict.values())
    cfg_list = [dict(zip(tmp_dict.keys(), c)) for c in combinations]
    # use only the first item of the list for testing
    dict_pars = cfg_list[0]
    n_reps_sims = dict_pars.pop("n_reps_sims")
    stat_type = dict_pars.pop("stat_type")
    mc_method = dict_pars.pop("mc_method")
    smpl, perms, theta = rmsim.create_sample(**dict_pars)
    return dict_pars, smpl, perms, theta


def test_samplegeneration():
    dict_pars, smpl, perms, theta = help_call_datagen()
    assert smpl.shape[0] == dict_pars["n_phase_bins"]
    assert smpl.shape[1] == dict_pars["n_subjects"]
    assert perms.shape[0] == dict_pars["n_phase_bins"]
    assert perms.shape[1] == dict_pars["n_subjects"]
    assert perms.shape[2] == dict_pars["n_perms"]


def test_circperm():
    dict_pars, smpl, perms, theta = help_call_datagen()
    aggr_dat = np.concatenate((perms, smpl[:, :, np.newaxis]), axis=2)
    z_along_permutations = zscore(aggr_dat, axis=2)
    # compute t statistics
    permemp = (
        np.sqrt(dict_pars["n_subjects"])
        * z_along_permutations.mean(axis=1)
        / z_along_permutations.std(axis=1)
    )
    # select empirical and permuted data, and run circular permutation on it
    emp = permemp[:, -1]
    perm = permemp[:, :-1]
    [bounds, summary] = rmstats.circ_perm(emp, perm, theta, alternative="two_sided")
    print(summary)
    assert len(summary["idxs"]) == 2

def test_preprocessing():

    # Parameters
    # other configs that might be put as input in some future versions
    # for the time being brutally defaulted here
    orgfs = 600; # original sampling frequency of raw data
    fs = 100; # final sampling frequency after downsampling
    downsmpl_rate = int(orgfs/fs) # if >13, not recomendations from sp.signal.decimate doc
    win_len = int(.4*fs)
    nbin = 60

    # Load and display data
    resp_trace = loadmat("../_exampledata/resp_raw_001.mat")["data"][0]
    # load behavioral events
    events = loadmat("../_exampledata/events.mat")["events"][0][0][0][0]
    # load template data
    template_HR = np.load("./tests/template_HR.npy")

    # 1. Preprocessing
    # 1.1
    x_norm = zscore(resp_trace)
    # 1.2 Downsampling
    x_down = decimate(x_norm, downsmpl_rate)
    # 1.3
    x_smooth = np.convolve(x_down, np.ones(win_len), 'same')/win_len
    # 1.4 Resample behavioral events
    dict_events = {"idx_sample" : np.asarray(events[0])//downsmpl_rate,
                   "HvsM" : events[1],
                  }
    dict_events.update({"resp" : x_smooth[dict_events["idx_sample"]]})
    # 2. Phase computation

    phase_vect = two_point_interp(x_smooth)
    dict_events.update({"resp_phase_angle" : phase_vect[dict_events["idx_sample"]]})

    # 3. Bin data
    phw = (np.pi/10)                                 # The angular half-width of the bin (the distance center-to-edge).
    pb = np.linspace(-np.pi,np.pi-2*np.pi/nbin,nbin) # vector containing centre of each phase bin (= phase bin vector)
    bin_idxs = create_phase_bins(dict_events["resp_phase_angle"], pb, phw)
    computed_HR = []
    for ibin in bin_idxs:
        computed_HR.append(np.mean(dict_events["HvsM"][ibin]))

    computed_HR = np.asarray(computed_HR)
    diff_tmp = np.abs(template_HR-computed_HR).sum()

    assert diff_tmp < 1e6
