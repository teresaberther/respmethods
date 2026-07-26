#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 15:31:50 2026

@author: ebalestr
"""


import pandas as pd
from pathlib import Path
import random
import pickle

# HQ
from hyperqueue import Client, Job
from hyperqueue.ffi.protocol import ResourceRequest

# pack configs combinations
from itertools import product

# respmethods
import Simulation as sim


def split_list(data: list[dict], n: int) -> list[list[dict]]:
    """
    Splits a list of dictionaries into N sublists of equal (or near-equal) size.

    Args:
        data: The list of dictionaries to split.
        n:    The number of sublists to create.

    Returns:
        A list of N sublists.

    Raises:
        ValueError: If n is less than 1 or greater than the length of data.
    """
    if n < 1:
        raise ValueError("n must be at least 1.")
    if n > len(data):
        raise ValueError("n cannot be greater than the length of the list.")

    # Base size of each chunk and how many chunks get an extra element
    chunk_size, remainder = divmod(len(data), n)

    result = []
    start = 0
    for i in range(n):
        # The first 'remainder' chunks get one extra element
        end = start + chunk_size + (1 if i < remainder else 0)
        result.append(data[start:end])
        start = end

    return result


#%%

sim_name = "par-space-1M"

fc = sim.create_dict_from_toml(f"{sim_name}.toml")

# generation of all the parameters combination
combinations = product(*fc.values())
cfg_list = [dict(zip(fc.keys(), c)) for c in combinations]
random.shuffle(cfg_list)

tot_combs = len(cfg_list)
n_procs, job_reps = 16, 45
njobs = int(tot_combs / (job_reps*n_procs))

jobs_list = split_list(cfg_list, njobs)
cfg_path_list = []
for idx, cfg in enumerate(jobs_list):
    
    job_fold = f"./{sim_name}/job-{idx}"
    Path(job_fold).mkdir(parents=True,
                         exist_ok=True)
    
    fname = f"{job_fold}/cfg.pkl"
    with open(fname, "wb") as f:
        pickle.dump(cfg, f)
    
    cfg_path_list.append(job_fold)
        
#%% HQ submission

client = Client()

################################################################
cpu_res = ResourceRequest(cpus=n_procs)

n_batches = 7
mini_batches = split_list(cfg_path_list, n_batches) 

for ibatch in mini_batches: 

    task_list = []
    for job_path in ibatch:

        cmd = ["uv", "run", "exec-job.py", job_path, str(n_procs)]
        job = Job(max_fails=5)

        job.program(cmd,
                    resources=cpu_res,
                    )    
        task_list.append(client.submit(job))

    client.wait_for_jobs(task_list)

#%% concatenate all DF in one
# and resubmit jobs that did not produce output

dfs, repeat_task_list = [], []
exceptions = 0
for job_path in cfg_path_list:
    file = f"{job_path}/DF.csv"
    try:
        df = pd.read_csv(file)
        dfs.append(df)
    except:
        print(f"  ❌ {file} ")
        print("Resubmitting...")
        cmd = ["uv", "run", "exec-job.py", job_path, str(n_procs)]
        job = Job()
        job.program(cmd,
                    resources=cpu_res,
                    )    
        repeat_task_list.append(client.submit(job))
        exceptions += 1

final_df = pd.concat(dfs, ignore_index=True)
final_df.to_csv(f"{sim_name}-failed-{exceptions}.csv", index=False)

client.wait_for_jobs(repeat_task_list)

#%% last round

dfs =  []
exceptions = 0
for job_path in cfg_path_list:
    file = f"{job_path}/DF.csv"
    try:
        df = pd.read_csv(file)
        dfs.append(df)
    except:
        print(f"  ❌ {file} ")
        exceptions += 1

final_df = pd.concat(dfs, ignore_index=True)
final_df.to_csv(f"{sim_name}-run2-failed-{exceptions}.csv", index=False)
