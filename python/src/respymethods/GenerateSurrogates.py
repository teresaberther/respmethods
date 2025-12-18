#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate surrogate data
"""

import numpy as np
from pyfftw.interfaces.scipy_fft import rfft, irfft

def iaaft(indat, maxiter=int(1e5), convthreshold=1e-6):
    """
    The function generates, via Iterative Amplitude Adjusted Fourier Transform 
    (IAFFT) 1 surrogatewith the same spectral characteristics of the original 
    data.
    
    Parameters
    ----------
    indat : array of M timepoints 
    maxiter : max iterations for the algorithm to converge (default=1E5)

    Returns
    ----------
    surrdat : array of M timepoints
    
    
    References
    
    - Schreiber, T., & Schmitz, A. (1996). Improved surrogate data for 
      nonlinearity tests. Physical review letters, 77(4), 635.
    - Lancaster, G., Iatsenko, D., Pidde, A., Ticcinelli, V., & Stefanovska, 
      A. (2018). Surrogate data for hypothesis testing of physical systems. 
      Physics Reports, 748, 1-60.
      
    """
    
    
    npt = len(indat)
    seqpt = np.arange(npt)

    # original amplitude values
    orig_amp = np.abs(rfft(indat, axis=0))

    # original idxs values
    orig_srt = np.sort(indat, axis=0)
        
    shffld_idxs = np.random.permutation(seqpt)
    shffld_array = indat[shffld_idxs] 
        
    itercount, mismatch = 0, 1    
    while (itercount<int(maxiter)) and (mismatch>convthreshold):
                
        fft_shffld = rfft(shffld_array, axis=0)
        phi_shffld = np.arctan2(fft_shffld.imag, fft_shffld.real)
        unit_phvect = np.exp(phi_shffld * 1j)
        z_n = np.real(irfft(orig_amp * unit_phvect, axis=0)) # approximate to machine precision
        
        z_idxs = np.argsort(z_n, axis=0)
        z_n[z_idxs] = orig_srt
            
        # alternative stop rule (adherent to original paper): stop if the order
        # does not change any longer
        # p_diff = ((shffld_idxs != z_idxs).sum(axis=0))/npt    
        # mismatch = p_diff.sum()                        
        # shffld_idxs = z_idxs
        
        # stop rule
        diff_ = ((shffld_array-z_n)**2).sum(axis=0)
        mismatch = diff_/(z_n**2).sum(axis=0)
        
        shffld_array = z_n
        
        itercount+=1
            
    return shffld_array




