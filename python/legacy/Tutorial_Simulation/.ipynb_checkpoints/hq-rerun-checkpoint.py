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
import os

# HQ
from hyperqueue import Client, Job
from hyperqueue.ffi.protocol import ResourceRequest

# pack configs combinations
from itertools import product

# respmethods
import Simulation as sim



sim_name = "par-space-1M"
n_procs, job_reps = 16, 45
nques = 96

SIMFOLD = f"./{sim_name}"
missing, all_paths = [], []
for subfolder in os.listdir(SIMFOLD):
    path = os.path.join(SIMFOLD, subfolder)
    if os.path.isdir(path):
        all_paths.append(path)
    if os.path.isdir(path) and not os.path.exists(os.path.join(path, "DF.csv")):
        missing.append(path)

print(f"{len(missing)}/{len(all_paths)} will be run")

#%% HQ submission

client = Client()
cpu_res = ResourceRequest(cpus=n_procs)

task_list = []
job_acc, batch_acc = 0, 0
for job_path in missing:

    job_acc += 1
    
    cmd = ["uv", "run", "exec-job.py", job_path, str(n_procs)]
    job = Job()

    job.program(cmd,
                resources=cpu_res,
                )    
    task_list.append(client.submit(job))

    if (job_acc%nques == 0) or (job_acc==len(missing)):
        print(f"waiting execution of batch {batch_acc}")
        client.wait_for_jobs(task_list)
        batch_acc += 1
        task_list = []


#%% last round

dfs =  []
exceptions = 0
for job_path in all_paths:
    file = f"{job_path}/DF.csv"
    try:
        df = pd.read_csv(file)
        dfs.append(df)
    except:
        print(f"  ❌ {file} ")
        exceptions += 1

final_df = pd.concat(dfs, ignore_index=True)
final_df.to_csv(f"{sim_name}-failed-{exceptions}.csv", index=False)
