#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import numpy as np
from scipy.stats import ecdf


def circ_perm(EMPdat, PERMdat, theta, alternative="two_sided"):
    
    
    match alternative:        
        case "two_sided":
            low_p = .025
            up_p  = .975
        case "lesser":
            low_p = .05
            up_p  = 1
        case "greater":
            low_p = 0
            up_p  = .95
        case _:            
            raise TypeError("'{}' not recognized\n".format(alternative) + 
                            "alternatives: 'two_sided', 'lesser', 'greater' ")            
    
        
    PERMEMPdat = np.concatenate((PERMdat, EMPdat[:, np.newaxis]), axis=1)
    nbins = theta.shape[0]
    nperms = PERMdat.shape[1]

    low_bounds, up_bounds = [], []    
    for ibin in range(nbins):
        
        this_bin  = PERMEMPdat[ibin, :]
        ECDF = ecdf(this_bin)
        
        tmp_low = ECDF.cdf.quantiles[ECDF.cdf.probabilities<=low_p][-1]
        low_bounds.append(tmp_low)
        
        tmp_up = ECDF.cdf.quantiles[ECDF.cdf.probabilities>=up_p][0]
        up_bounds.append(tmp_up)

    low_bounds = np.array(low_bounds)
    up_bounds = np.array(up_bounds)

    surrogatedist = []
    for iperm in range(nperms):
        
        thisperm = PERMEMPdat[:, iperm]
        
        tmp_out = circ_clust(theta, thisperm, low_bounds, up_bounds)
        
        all_cluststats = np.array(tmp_out["low_cluststats"] + 
                                  tmp_out["up_cluststats"])
    
        if len(all_cluststats)>0:
            maxabs_idx = np.argmax(np.abs(all_cluststats))
            surr_val = all_cluststats[maxabs_idx]
        else:
            maxabs_idx = np.argmax(np.abs(thisperm))
            surr_val = thisperm[maxabs_idx]
            
        surrogatedist.append(surr_val)

    # clustering stats for the actual data
    out = circ_clust(theta, EMPdat, low_bounds, up_bounds)

    all_cluststats = np.array(out["low_cluststats"] + out["up_cluststats"])
    all_cluststidxs = out["low_clusts_idxs"] + out["up_clusts_idxs"]

    summary = {"idxs" : [], "clustmass_stat" : [], "p_ecdf" : [],
               "p" : []}
    
    acc=0
    for imassstat in all_cluststats:
        
        tmpstats = np.array(surrogatedist + [imassstat] + [np.inf])
        p_ecdf = np.mean(tmpstats <= imassstat)
        p = np.mean(np.abs(tmpstats) > np.abs(imassstat)) # again, prevent p=0 by doing the compoarison against a vector containing inf.
        
        summary["idxs"].append(all_cluststidxs[acc])
        summary["clustmass_stat"].append(imassstat)
        summary["p_ecdf"].append(p_ecdf)
        summary["p"].append(p)
                
        acc+=1
    
    bounds = {"upper": up_bounds,
              "lower": low_bounds}
    
    
    return bounds, summary

def circ_clust(theta, target, low_bounds, up_bounds, threshtype="empirical"):
    
    # classify points below and above the thresholds
    below_ptile = target<low_bounds
    low_clusts = _local_clust(below_ptile, theta)
    above_ptile = target>up_bounds
    up_clusts = _local_clust(above_ptile, theta)
    
    # compute statistics for clusters above the upper threshold
    up_cluststats = []; 
    for iclust in up_clusts:
        up_cluststats.append(np.sum(target[iclust]))

    # compute statistics for clusters below the lower threshold
    low_cluststats = []; 
    for iclust in low_clusts:
        low_cluststats.append(np.sum(target[iclust]))
        
    # output dictionary with clustering results
    out = {"up_clusts_idxs" : up_clusts,
           "low_clusts_idxs" : low_clusts,
           "up_cluststats" : up_cluststats,
           "low_cluststats" : low_cluststats}
    
    return out


def _local_clust(masklength, anglevect):
    
    if not np.any(masklength):
        clusts = []
    else:
        cmplx_repr = masklength * np.exp(anglevect * 1j)        
        
        [X, Y] = np.meshgrid(cmplx_repr, cmplx_repr)        
        dists = np.round(np.abs(X-Y), 5)
        minsigdist = np.unique(dists)
        
        shortlenarc = dists==minsigdist[1];
        [row, col] = np.where(shortlenarc)

        acc = 0
        clusts = []
        
        while col.size > 0:
            
            idxs_acc = [];
            seed_entry = [int(col[0])]
            
            while True:
                
                idxs_acc += seed_entry
                tmp_mask = np.isin(col, seed_entry)
                
                if not np.any(tmp_mask):
                    break
                else:
                    seed_entry = np.unique(row[tmp_mask]).tolist()
                    row = np.delete(row, tmp_mask)
                    col = np.delete(col, tmp_mask)
            
            clusts.append(np.unique(np.array(idxs_acc)).tolist())
            acc += 1
            
    return clusts





