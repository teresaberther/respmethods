#!/usr/bin/env python3
"""
RespStats module for respymethods. It icludes:
- circ_perm
- circ_clust
- omni_perm
- _local_clust
"""

import numpy as np
from scipy.stats import ecdf


def circ_perm(EMPdat, PERMdat, theta, alternative):
    """
    circ_perm performs a cluster permutation test on circular data.

    Args:
        EMPdat (np.ndarray): Empirical data array.
        PERMdat (np.ndarray): Permutation/surrogate data array.
        theta (np.ndarray): The coordinate or binning vector.
        alternative (str): The type of test to perform.
                           Options: 'two_sided', 'less', 'greater'.

    Returns:
        tuple: A tuple containing:
            - bounds (dict): Dictionary with 'upper' and 'lower' boundary arrays.
            - summary (dict): Dictionary containing identified cluster indices ('idxs'),
                              their statistics ('clustmass_stat'), and calculated
                              p-values ('p_ecdf', 'p').

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """
    match alternative:
        case "two_sided":
            low_p = 0.025
            up_p = 0.975
        case "less":
            low_p = 0.05
            up_p = 1
        case "greater":
            low_p = 0
            up_p = 0.95
        case _:
            raise TypeError(
                "'{}' not recognized\n".format(alternative)
                + "alternatives: 'two_sided', 'less', 'greater' "
            )

    PERMEMPdat = np.concatenate((PERMdat, EMPdat[:, np.newaxis]), axis=1)
    nbins = theta.shape[0]
    nperms = PERMdat.shape[1]

    low_bounds, up_bounds = [], []
    for ibin in range(nbins):
        this_bin = PERMEMPdat[ibin, :]
        ECDF = ecdf(this_bin)

        tmp_low = ECDF.cdf.quantiles[ECDF.cdf.probabilities <= low_p][-1]
        low_bounds.append(tmp_low)

        tmp_up = ECDF.cdf.quantiles[ECDF.cdf.probabilities >= up_p][0]
        up_bounds.append(tmp_up)

    low_bounds = np.array(low_bounds)
    up_bounds = np.array(up_bounds)

    surrogatedist = []
    for iperm in range(nperms):
        thisperm = PERMEMPdat[:, iperm]

        tmp_out = circ_clust(theta, thisperm, low_bounds, up_bounds)

        all_cluststats = np.array(tmp_out["low_cluststats"] + tmp_out["up_cluststats"])

        if len(all_cluststats) > 0:
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

    summary = {"idxs": [], "clustmass_stat": [], "p_ecdf": [], "p": []}

    acc = 0
    for imassstat in all_cluststats:
        tmpstats = np.array(surrogatedist + [imassstat] + [np.inf])
        p_ecdf = np.mean(tmpstats <= imassstat)
        p = np.mean(
            np.abs(tmpstats) > np.abs(imassstat)
        )  # again, prevent p=0 by doing the compoarison against a vector containing inf.

        summary["idxs"].append(all_cluststidxs[acc])
        summary["clustmass_stat"].append(imassstat)
        summary["p_ecdf"].append(p_ecdf)
        summary["p"].append(p)

        acc += 1

    bounds = {"upper": up_bounds, "lower": low_bounds}

    return bounds, summary


def circ_clust(theta, target, low_bounds, up_bounds, threshtype="empirical"):
    """
    Identify contiguous clusters of data points that exceed specified upper or lower bounds.

    This function scans the target data to find sequences of adjacent points that fall
    outside the provided thresholds. For each identified cluster, it calculates the
    'cluster mass' (the sum of the values within that cluster).

    Args:
        theta (np.ndarray): The coordinate or binning vector (used by _local_clust
                            to determine adjacency).
        target (np.ndarray): The data array to be tested for clusters.
        low_bounds (np.ndarray): The lower threshold values for each bin.
        up_bounds (np.ndarray): The upper threshold values for each bin.
        threshtype (str, optional): Type of thresholding used. Defaults to "empirical".

    Returns:
        dict: A dictionary containing:
            - "up_clusts_idxs": List of index groups for clusters above the upper bound.
            - "low_clusts_idxs": List of index groups for clusters below the lower bound.
            - "up_cluststats": The sum of values (mass) for each upper cluster.
            - "low_cluststats": The sum of values (mass) for each lower cluster.

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """

    # classify points below and above the thresholds
    below_ptile = target < low_bounds
    low_clusts = _local_clust(below_ptile, theta)
    above_ptile = target > up_bounds
    up_clusts = _local_clust(above_ptile, theta)

    # compute statistics for clusters above the upper threshold
    up_cluststats = []
    for iclust in up_clusts:
        up_cluststats.append(np.sum(target[iclust]))

    # compute statistics for clusters below the lower threshold
    low_cluststats = []
    for iclust in low_clusts:
        low_cluststats.append(np.sum(target[iclust]))

    # output dictionary with clustering results
    out = {
        "up_clusts_idxs": up_clusts,
        "low_clusts_idxs": low_clusts,
        "up_cluststats": up_cluststats,
        "low_cluststats": low_cluststats,
    }

    return out


def _local_clust(masklength, anglevect):
    """
    Identify contiguous clusters of True values in a mask, accounting for circularity.

    This function treats the indices of the mask as points on a circle using the
    provided angle vector. It calculates the distance between all points to
    identify immediate neighbors and then iteratively groups these neighbors
    into distinct clusters.

    Args:
        masklength (np.ndarray): A boolean array (mask) where True indicates a point
                                 of interest to be clustered.
        anglevect (np.ndarray): An array of angles (in radians) corresponding to the
                                positions of the points, used to handle circular
                                wrap-around adjacency.

    Returns:
        list: A list of lists, where each inner list contains the indices of a
              contiguous cluster. Returns an empty list if no True values are found.

    Notes:
        The function converts angles to complex numbers (e^{i*theta}) to calculate
        chordal distances between points, ensuring that the first and last
        elements of the array are correctly identified as neighbors if they
        are spatially adjacent on the circle.

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """
    if not np.any(masklength):
        clusts = []
    else:
        cmplx_repr = masklength * np.exp(anglevect * 1j)

        [X, Y] = np.meshgrid(cmplx_repr, cmplx_repr)
        dists = np.round(np.abs(X - Y), 5)
        minsigdist = np.unique(dists)

        if len(minsigdist) < 3:
            row, col = np.array(masklength), np.array(masklength)
        else:
            shortlenarc = dists == minsigdist[1]
            row, col = np.where(shortlenarc)

        acc = 0
        clusts = []

        while col.size > 0:
            idxs_acc = []
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


def omni_perm(EMPdat, PERMdat, theta):
    """
    Perform an omnibus permutation test to identify individual significant points.

    Unlike cluster-based tests, this function evaluates the significance of individual bins.
    It corrects for multiple comparisons by comparing the global maximum and minimum values across
    all bins and permutations against the 97.5th and 2.5th percentiles, respectively.

    Args:
        EMPdat (np.ndarray): Empirical data array.
        PERMdat (np.ndarray): Permutation/surrogate data array.
        theta (np.ndarray): The coordinate or binning vector.

    Returns:
        tuple: A tuple containing:
            - bounds (dict): Dictionary with 'upper' and 'lower' global thresholds
                              (broadcasted to the shape of theta).
            - summary_omni (dict): Dictionary containing the indices of significant
                                   points ('idxs'), their values ('clustmass_stat'),
                                   and placeholder p-values ('p', 'p_ecdf').

    Notes:
        This function is designed for compatibility with cluster-based output formats,
        treating each significant point as a cluster of size one.

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """
    PERMEMPdat = np.concatenate((PERMdat, EMPdat[:, np.newaxis]), axis=1)

    up_dist = PERMEMPdat.max(axis=0)
    low_dist = PERMEMPdat.min(axis=0)

    up_ECDF = ecdf(up_dist)
    low_ECDF = ecdf(low_dist)

    low_bound = low_ECDF.cdf.quantiles[low_ECDF.cdf.probabilities <= 0.025][-1]
    up_bound = up_ECDF.cdf.quantiles[up_ECDF.cdf.probabilities >= 0.975][0]

    up_idxs = list(np.where(EMPdat > up_bound)[0])
    low_idxs = list(np.where(EMPdat < low_bound)[0])
    all_idxs = up_idxs + low_idxs

    # for compatibility
    summary_omni = {
        "clustmass_stat": [],
        "idxs": [],
        "p": [],
        "p_ecdf": [],
    }

    for idx in all_idxs:
        summary_omni["clustmass_stat"].append(EMPdat[idx])
        summary_omni["idxs"].append([idx])
        summary_omni["p"].append(np.float64(0.001))
        summary_omni["p_ecdf"].append(np.float64(0.001))

    bounds = {
        "upper": np.ones(theta.shape) * up_bound,
        "lower": np.ones(theta.shape) * low_bound,
    }

    return bounds, summary_omni
