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

NPERM = 10

# detrend data (I preferred not to bloat the IAAFT with stuff that can be more 
# easily and intuitively implemented by the user externally)
indat = detrend(binhr, axis=0)

perm_mat = IAAFT(indat, NPERM)


# prelim stats
nsubjs = indat.shape[1]
tvals = np.sqrt(nsubjs)*indat.mean(axis=1)/indat.std(axis=1)

# sketch localclust here
masklength = np.abs(tvals) > 2
anglevect = pb.copy()

masklength = 1*masklength[:, np.newaxis]
cmplx_repr = masklength*np.exp(anglevect * 1j)

X, Y = np.meshgrid(cmplx_repr, cmplx_repr)
dists = np.round(np.abs(X-Y), 5)
shortlenarc = dists == np.unique(dists)[1]

