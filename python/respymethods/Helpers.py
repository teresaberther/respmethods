#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set of Helpers functins
"""

import numpy as np

def plv(theta1, theta2):
    """
    The function computes the instantaneous phase of two signals and returns the
    Phase Locking Value (PLV) between the two as:
                    
                        PLV = (1/N) sum(exp(i(theta1-theta2)))    
    
    PLV ranges from 0 (no synchrony) to 1 (perfect synchrony)
    
    Parameters
    ----------
    theta1 : array 1 of M phase angles
    theta2 : array 2 of M phase angles

    Returns
    ----------
    PLV    
    
    References
    
    - Lachaux, J. P., Rodriguez, E., Martinerie, J., & Varela, F. J. (1999). 
      Measuring phase synchrony in brain signals. Human brain mapping, 8(4), 194-208.
      
    """
        
    N = len(theta1)
    phase_diff = theta1-theta2
    PLV = np.abs(np.sum(np.exp(1j*phase_diff)))/N;
    
    return PLV


def prune_nan(sig1, sig2):
    """
    The function selects the longest numeric (aka, non-NAN) segment common 
    between two signals sig1 and sig2.
    """
    
    num_idx1 =  np.where(~np.isnan(sig1))[0]
    num_idx2 =  np.where(~np.isnan(sig2))[0]

    common_onset = np.max([np.min(num_idx1), np.min(num_idx2)])
    common_offset = np.min([np.max(num_idx1), np.max(num_idx2)])

    out_sig1 = sig1[common_onset:common_offset]
    out_sig2 = sig2[common_onset:common_offset]
    
    return out_sig1, out_sig2
    