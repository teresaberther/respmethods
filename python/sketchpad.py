#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  8 11:35:06 2025

@author: ebalestr
"""

import numpy as np
from scipy.io import loadmat
from scipy.fft import fft, ifft
from scipy.signal import detrend
import time, copy
import matplotlib.pyplot as plt

# custom
from GenerateSurrogates import IAAFT

binhr = loadmat("../_exampledata/binhr.mat")["binhr"].T
pb = loadmat("../_exampledata/phasebinvect.mat")["pb"]
prev_shuffled = loadmat("../_exampledata/surrbinhr.mat")["sbinhr"]

# quicktest: scipy fft or numpy fft?
# result: scipy fft slightly better in both accuracy & speed. 

QUICKTEST=False

if QUICKTEST:

    NPERM = int(1e6)
    
    s = time.process_time()
    ERR_A = []
    strt_hr = copy.deepcopy(binhr)
    for n in range(NPERM):
        A = np.fft.fft(strt_hr)
        strt_hr = np.real(np.fft.ifft(A))
        ERR_A.append(np.abs(binhr - strt_hr).sum())
    print(time.process_time()-s)
    print(np.mean(ERR_A))
    print(ERR_A[-1])
    
    s = time.process_time()
    ERR_B = []
    strt_hr = copy.deepcopy(binhr)
    for n in range(NPERM):
        B = fft(strt_hr)
        strt_hr = np.real(ifft(B))
        ERR_B.append(np.abs(binhr - strt_hr).sum())
    
    print(time.process_time()-s)
    print(np.mean(ERR_B))
    print(ERR_B[-1])

# create surrogates

NPERM = 100
indat = detrend(binhr)

perm_mat = IAAFT(indat, NPERM)
