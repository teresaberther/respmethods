import numpy as np
from scipy.stats import zscore, iqr
import pandas as pd
from copy import copy

# pack configs combinations
from itertools import product

# respmethods
import Simulation as sim
from RespStats import circ_perm, omni_perm

from multiprocessing import Pool
import random


def run_batch_sims(batch_sims):
    
    updated_dicts = []
    for idict in batch_sims:
        
        n_reps_sims = idict.pop("n_reps_sims")
        stat_type = idict.pop("stat_type")
        mc_method = idict.pop("mc_method")
        
        for irep in range(n_reps_sims):
        
            tmp_dict = copy(idict)    
            
            try:
                smpl, perms, theta = sim.create_sample(**tmp_dict)                
                # note: dimensions here are phase, subj, perms (differently from the Tutorial)
                
                aggr_dat = np.concatenate((perms, smpl[:, :, np.newaxis]), 
                                          axis=2)
                tmp_z = zscore(aggr_dat, axis=2)
            
                match stat_type.lower():             
                    case "tval":
                        permemp = (np.sqrt(tmp_dict["n_subjects"])*
                                   tmp_z.mean(axis=1)/tmp_z.std(axis=1)) # along axis 0 in Tutorial...
                    case "robust":
                        permemp = np.median(tmp_z, axis=1)/iqr(tmp_z, axis=1)

            
                emp = permemp[:, -1]
                perm = permemp[:, 0:-1]
            
                match mc_method.lower():                
                    case "circ":                                    
                        [bounds, summary] = circ_perm(emp, perm, theta)
                    case "omni":
                        [bounds, summary] = omni_perm(emp, perm, theta)
                        
                out = sim.sdt(summary, tmp_dict["effect_mag"], 
                              tmp_dict["n_phase_bins"])    
            except:
                out = {"H" : np.nan,
                       "FA" : np.nan}            
        
            for key, value in out.items():
                tmp_dict[key] = value
                
            tmp_dict["rep_idx"] = irep
            tmp_dict["stat_type"] = stat_type
            tmp_dict["mc_method"] = mc_method
            
        
            updated_dicts.append(tmp_dict) 
            
    DF = pd.DataFrame(updated_dicts)
    
    return DF


def distribute_sims(cfg_list, n_procs):

    if len(cfg_list)%(n_procs):
        n_procs -= 1
    
    n_simXproc = len(cfg_list)//(n_procs)

    batched_sims = []
    for iproc in range(n_procs):        
        single_batch = []
        for isim in range(n_simXproc):
            this_sim = cfg_list.pop()
            single_batch.append(this_sim)
        batched_sims.append(single_batch)
    
    if cfg_list:
        batched_sims.append(cfg_list)

    return batched_sims


if __name__ == "__main__":

    n_procs = 248
    compose_name = f"palma-{n_procs}"
    # n_procs = 10
    # compose_name = "sim-cfg"


    fc = sim.create_dict_from_toml(f"{compose_name}.toml")
    
    # generation of all the parameters combination
    combinations = product(*fc.values())
    cfg_list = [dict(zip(fc.keys(), c)) for c in combinations]
    random.shuffle(cfg_list)
    
    batched_sims = distribute_sims(cfg_list, n_procs)

    with Pool(len(batched_sims)) as p:
        DFs_bag = p.map(run_batch_sims, batched_sims)

    DF = pd.concat(DFs_bag)
    DF.to_csv(f"{compose_name}.csv")
