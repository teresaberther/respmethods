#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 16:47:28 2025

@author: ebalestr
"""

from scipy.io import loadmat
from scipy.stats import zscore
import matplotlib.pyplot as plt
from respymethods.RespStats import circ_perm
import numpy as np


tmp_surr = loadmat("data/surrbinhr.mat")
sbinhr = tmp_surr["sbinhr"]

nsubj, ntheta, nperm = sbinhr.shape

tmp_dat = loadmat("data/binhr.mat")
binhr = tmp_dat["binhr"]

theta = np.linspace(-np.pi, np.pi, ntheta)

aggr_dat = np.concatenate((sbinhr, binhr[:, :, np.newaxis]), axis=2)

tmp_z = zscore(aggr_dat, axis=2)

permemp_t = np.sqrt(nsubj)*tmp_z.mean(axis=0)/tmp_z.std(axis=0)

emp_t = permemp_t[:, -1]
perm_t = permemp_t[:, 0:-1]

[summary, bounds] = circ_perm(emp_t, perm_t, theta)

