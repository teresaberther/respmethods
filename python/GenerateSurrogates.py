#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate surrogate data
"""

import numpy as np
from scipy.fft import fft, ifft
from tqdm import tqdm

def IAAFT(indat, nsurr, maxiter=int(1e5)):
    """
    The function generates, via Iterative Amplitude Adjusted Fourier Transform 
    (IAFFT) n surrogates of time series with the same spectral characteristics 
    of the original data.
    
    Parameters
    ----------
    indat : array of M timepoints x N observations (e.g. subjects)
    nsurr : number of surrogates to be generated
    maxiter : max iterations for the algorithm to converge (default=1E5)

    Returns
    ----------
    surrdat : array of M timepoints x N observation x S surrogates 
    
    
    References
    
    - Schreiber, T., & Schmitz, A. (1996). Improved surrogate data for 
      nonlinearity tests. Physical review letters, 77(4), 635.
    - Lancaster, G., Iatsenko, D., Pidde, A., Ticcinelli, V., & Stefanovska, 
      A. (2018). Surrogate data for hypothesis testing of physical systems. 
      Physics Reports, 748, 1-60.
      
    """
    
    
    npt, nent = indat.shape
    seqpt = np.arange(npt)

    # initialize the array of surrogate data with 0s 
    surrdat = np.repeat(0*indat[:, :, np.newaxis], nsurr, axis=-1) 

    # original amplitude values
    orig_amp = np.abs(fft(indat, axis=0))

    # original idxs values
    orig_srt = np.sort(indat, axis=0)

    for isurr in tqdm(range(nsurr)):
        
        shffld_idxs, shffld_array = np.zeros(indat.shape, dtype=int), np.zeros(indat.shape)

        for ient in range(nent):
            
            mask_pos_ent = np.random.permutation(seqpt)
            shffld_idxs[:, ient] = mask_pos_ent
            shffld_array[:, ient] = indat[shffld_idxs[:, ient], ient] 
            
        itercount, mismatch = 0, 1    
        while (itercount<int(maxiter)) and (mismatch!=0):
                    
            fft_shffld = fft(shffld_array, axis=0)
            phi_shffld = np.angle(fft_shffld)
            unit_phvect = np.exp(phi_shffld * 1j)
            z_n = np.real(ifft(orig_amp * unit_phvect, axis=0)) # approximate to machine precision
            
            z_idxs = np.argsort(z_n, axis=0)
            
            for ient in range(nent):
                z_n[z_idxs[:, ient], ient] = orig_srt[:, ient]
                
            p_diff = ((shffld_idxs != z_idxs).sum(axis=0))/npt    
            mismatch = p_diff.sum()                        
            shffld_idxs = z_idxs.copy()
            
            # # previous stop rule (Trento people)
            # diff_ = ((shffld_array-z_n)**2).sum(axis=0)
            # diff_norm = diff_/(z_n**2).sum(axis=0)
            # print(diff_norm)
            
            shffld_array = z_n.copy()
            
            itercount+=1
            
        surrdat[:, :, isurr] = shffld_array

    return surrdat