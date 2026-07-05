from itertools import product

import numpy as np
from scipy.stats import zscore

# respymethods imports
import respymethods.RespStats as rmstats
import respymethods.Simulation as rmsim


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
