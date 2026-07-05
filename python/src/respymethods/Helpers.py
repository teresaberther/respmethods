#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set of Helpers functins
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from cmcrameri import cm 

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


def plot_clusters(summary, x_axis, y_axis):

    if len(summary["idxs"])==0:
        print("no cluster found")
        return

    
    # 1. Plot the empirical results
    plt.plot(x_axis, y_axis, label="empirical results", color='black', linewidth=2, zorder=1)

    # 2. Plot the cluster tick lines
    cmap = cm.batlow

    # 3. Create an array of positions from 0.0 to 1.0
    positions = np.linspace(0, 1, len(summary["idxs"]))

    # 4. Pass the positions into the colormap to get your array of RGBA colors
    colors = cmap(positions)
    
    # Calculate a small offset so the ticks don't sit exactly on the line
    # This ensures they are "clearly separated"
    offset = np.abs(np.min(y_axis)) * 0.1 

    for i, idxs in enumerate(summary["idxs"]):
                
        # Determine the Y value for the tick (average Y of the cluster + offset)
        clust_vals = np.zeros(x_axis.shape)*np.nan
        clust_vals[idxs] = np.min(y_axis) - offset
        
        # Plot the horizontal line
        p_val = summary["p"][i]
        plt.plot(x_axis, clust_vals, 
                   color=colors[i], linewidth=3, 
                   label=f"Cluster {i+1} (p={p_val:.4f})", zorder=2)

    plt.xlabel("theta")
    plt.ylabel("data")
    plt.title("Empirical Results with Significant Clusters")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left') # Move legend outside
    plt.tight_layout()
    plt.show()