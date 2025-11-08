#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate surrogate data
"""

import numpy as np
from scipy.fft import fft, ifft

def IAAFT(indat, nperm, maxiter=1e4, converr=0):
    
    npt, nent = indat.shape
    seqpt = np.arange(npt)

    # initialize the array of surrogate data with 0s 
    perm_mat = np.repeat(0*indat[:, :, np.newaxis], nperm, axis=-1) 

    # original amplitude values
    orig_amp = np.abs(fft(indat, axis=0))

    # original idxs values
    orig_srt = np.sort(indat, axis=0)

    for iperm in range(nperm):
        
        shffld_idxs, shffld_array = np.zeros(indat.shape, dtype=int), np.zeros(indat.shape)

        for ient in range(nent):
            
            mask_pos_ent = np.random.permutation(seqpt)
            shffld_idxs[:, ient] = mask_pos_ent
            shffld_array[:, ient] = indat[shffld_idxs[:, ient], ient] 
            
        itercount, mismatch = 0, 1    
        while (itercount<int(maxiter)) and (mismatch!=converr):
                    
            fft_shffld = fft(shffld_array, axis=0)
            phi_shffld = np.angle(fft_shffld)
            unit_phvect = np.exp(phi_shffld * 1j)
            z_n = np.real(ifft(orig_amp * unit_phvect, axis=0)) # approximate to machine precision
            
            z_idxs = np.argsort(z_n, axis=0)
            
            for ient in range(nent):
                z_n[z_idxs[:, ient], ient] = orig_srt[:, ient]
                
            p_diff = ((shffld_idxs != z_idxs).sum(axis=0))/npt    
            mismatch = p_diff.sum()
            
            # #debug
            # print(p_diff)
            
            shffld_idxs = z_idxs.copy()
            
            # # previous stop rule (Trento people)
            # diff_ = ((shffld_array-z_n)**2).sum(axis=0)
            # diff_norm = diff_/(z_n**2).sum(axis=0)
            # print(diff_norm)
            
            shffld_array = z_n.copy()
            
            itercount+=1
            
        perm_mat[:, :, iperm] = shffld_array

    return perm_mat