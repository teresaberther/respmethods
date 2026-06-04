import tomllib as toml

import numba as nb
import numpy as np


@nb.vectorize
def nb_binom(n, p):
    # numba wrapper for binomial
    return np.random.binomial(n, p)


@nb.njit
def model_resp(n_trials, effect_mag, phase_event, phase_lag, trl_noise):
    """
    Models the response probability based on a cosine-modulated (respiration) phase effect.

    The function calculates a probability of response by combining a baseline (0.5),
    a periodic effect (cosine), and random Gaussian noise.

    Args:
        n_trials (int): Number of trials to simulate.
        effect_mag (float): The amplitude/strength of the effect.
        phase_event (float): The phase at which the event happens (radians).
        phase_lag (float): The shift or lag applied to the phase, in radians. Used for expressing inter subjects variability.
        trl_noise (float): The standard deviation of the random noise added to each trial.

    Returns:
        ndarray: An array of response probabilities for each trial.

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """

    p_resp = (
        0.5
        + effect_mag * np.cos(phase_event + phase_lag)
        + np.random.randn(n_trials) * trl_noise
    )

    return p_resp


@nb.njit
def observer(n_trials, effect_mag, phase_lag, trl_noise):
    """
    Simulates the behavior of a single observer across multiple trials.

    The function generates random event phases, calculates the probability of
    response based on those phases, ensures probabilities stay within [0, 1],
    and samples the final binary responses.

    Args:
        n_trials (int): Number of trials to simulate.
        effect_mag (float): Strength of the effect on the response.
        phase_lag (float): The phase shift/lag of the observer.
        trl_noise (float): Amount of random noise per trial.

    Returns:
        ndarray: A 2D array of shape (n_trials, 2) where:
                 Column 0 = the random phase of the event.
                 Column 1 = the observer's binary response (0 or 1).

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """

    # 1. Generate random phases for each trial (uniformly distributed between -π and π)
    phase_event = np.random.uniform(-np.pi, np.pi, n_trials)

    # 2. Calculate the probability of response using the modulation model
    p_resp = model_resp(n_trials, effect_mag, phase_event, phase_lag, trl_noise)

    # 3. Clip probabilities to prevent values outside the valid [0, 1] range
    p_resp[p_resp <= 0] = 0.01
    p_resp[p_resp >= 1] = 0.99

    # 4. Sample binary outcomes (0 or 1) based on the calculated probabilities
    resps = nb_binom(1, p_resp)
    resps = resps.astype(np.float64)

    # 5. Combine phases and responses into a single table (trials x 2)
    beh_events = np.stack((phase_event, resps)).T

    return beh_events


@nb.njit
def create_phase_bins(phase_vect, center_bin, range_bin):
    """
    Groups phase values into circular bins based on their proximity to center points.

    This function converts phases to Cartesian coordinates (unit circle) and
    determines which points fall within a specified angular range around
    the provided center bins.

    Args:
        phase_vect (ndarray): Array of phase values (in radians).
        center_bin (array_like): The center phase(s) for the bin(s) to be created.
        range_bin (float): The angular half-width of the bin (the distance from
                           center to edge).

    Returns:
        list: A list of arrays, where each array contains the indices of the
              phase_vect that fall within the corresponding bin.

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """

    x_vect = np.cos(phase_vect)
    y_vect = np.sin(phase_vect)

    plus_min_range = np.exp(np.array([range_bin, -range_bin]) * 1j)
    cmplx_cntr = np.exp(center_bin * 1j)

    bin_idxs = []
    idx_acc = 0
    for ibin in cmplx_cntr:
        x_ranges = np.real(ibin * plus_min_range)
        y_ranges = np.imag(ibin * plus_min_range)

        # correct for arc/axis intersection
        if (0 >= x_ranges.min()) & (0 <= x_ranges.max()):
            y_ranges = np.append(y_ranges, 1 * np.sign(y_ranges[0]))
        if (0 >= y_ranges.min()) & (0 <= y_ranges.max()):
            x_ranges = np.append(x_ranges, 1 * np.sign(x_ranges[0]))

        lgcl = (
            (x_vect >= x_ranges.min())
            & (x_vect < x_ranges.max())
            & (y_vect >= y_ranges.min())
            & (y_vect < y_ranges.max())
        )

        bin_idxs.append(np.where(lgcl)[0])
        idx_acc += 1

    return bin_idxs


@nb.njit
def pool_beh_resp(beh_events, bin_idxs):
    """
    Calculates the average response probability for each defined phase bin.

    This function iterates through the indices of each bin, extracts the
    corresponding responses from the behavior data, and computes the mean
    response for that specific angular segment.

    Args:
        beh_events (ndarray): A 2D array where column 0 is the phase
                              and column 1 is the binary response (0 or 1).
        bin_idxs (list): A list of arrays containing the indices of
                         trials that fall into each bin.

    Returns:
        ndarray: An array of mean responses (probabilities) for each bin.

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """

    mv_avg = []
    for idxs in bin_idxs:
        binned_events = beh_events[idxs, 1]
        mv_avg.append(np.mean(binned_events))

    return np.array(mv_avg)


@nb.njit
def create_sample(
    n_subjects,
    intersubj_noise,
    n_trials,
    trl_noise,
    effect_mag,
    range_bin,
    n_phase_bins,
    n_perms,
):
    """
    Simulates a complete dataset of multiple subjects and generates a
    permutation-based null distribution for each subject.

    The function simulates observers with varying phase lags, calculates their
    response rates across phase bins, and then shuffles the phase data
    (permutations) to determine a null hypothesis distribution.

    Args:
        n_subjects (int): Number of simulated participants.
        intersubj_noise (float): Log-scale noise for inter-subject variability.
        n_trials (int): Number of trials per subject.
        trl_noise (float): Log-scale noise for trial-by-trial variability.
        effect_mag (float): The strength of the modulation effect.
        range_bin (float): The bin width as a multiplier of pi.
        n_phase_bins (int): Number of angular bins to divide the circle into.
        n_perms (int): Number of permutations to run for the null distribution.

    Returns:
        smpl_hr (ndarray): Observed response rates [bins x subjects].
        perms_hr (ndarray): Permuted response rates [bins x subjects x permutations].
        center_bin (ndarray): The center phase of each bin.

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """

    trl_noise = np.exp(trl_noise)
    intersubj_noise = np.exp(intersubj_noise)
    center_bin = np.linspace(-np.pi, np.pi, n_phase_bins)
    range_bin *= 2 * np.pi  # range_bin is a pi multiplier
    smpl_hr = np.empty((n_phase_bins, n_subjects))
    perms_hr = np.empty((n_phase_bins, n_subjects, n_perms))

    for isubj in range(n_subjects):
        subj_lag = np.random.randn() * intersubj_noise * np.pi
        exp_obs = observer(n_trials, effect_mag, subj_lag, trl_noise)
        binned = create_phase_bins(exp_obs[:, 0], center_bin, range_bin)
        hr_obs = pool_beh_resp(exp_obs, binned)
        smpl_hr[:, isubj] = hr_obs

        tmp_bins = exp_obs[:, 0]
        for iperm in range(n_perms):
            np.random.shuffle(tmp_bins)
            perm_binned = create_phase_bins(tmp_bins, center_bin, range_bin)
            ps_hr = pool_beh_resp(exp_obs, perm_binned)
            perms_hr[:, isubj, iperm] = ps_hr

    return smpl_hr, perms_hr, center_bin


def create_dict_from_toml(fname):
    """
    Loads simulation parameters from a TOML file and expands them into arrays.

    The function supports two ways of defining parameters in the TOML file:
    1. 'fixed_range': A literal list of values to be used.
    2. 'strtval', 'endval', 'n_entries': A range that will be expanded into
       a linear space (np.linspace).

    Args:
        fname (str): Path to the .toml configuration file.

    Returns:
        dict: A dictionary where keys are parameter names and values are
              NumPy arrays of the values to be tested.

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """

    cfg = {}
    ncombs_counter = 1
    with open(fname, "rb") as f:
        fc = toml.load(f)
        for cfg_name, cfg_values in fc.items():
            if "fixed_range" in cfg_values.keys():
                cfg.update({cfg_name: np.array(cfg_values["fixed_range"])})
            else:
                fr = np.linspace(
                    cfg_values["strtval"],
                    cfg_values["endval"],
                    num=cfg_values["n_entries"],
                )
                cfg.update({cfg_name: fr})
            ncombs_counter *= len(cfg[cfg_name])
        print(f"\n{ncombs_counter} combinations will be created. \n")
        print("They will be based on the following configuration: \n")
        print(cfg)

    return cfg


def sdt(summary, effect_mag, n_phase_bins):
    """
    Performs Signal Detection Theory (SDT) analysis to evaluate the
    detection of a theoretical effect.

    The function compares the statistically significant clusters found in
    the 'summary' data against a theoretical response curve, determined by model_resp.
    A 'Hit' (H) is recorded if a significant cluster matches the sign of the theoretical
    effect at the same phase location.

    Args:
        summary (dict): Results of the cluster analysis, containing:
            - "p": p-values for each cluster.
            - "idxs": Indices of the bins belonging to each cluster.
            - "clustmass_stat": The test statistic for each cluster.
        effect_mag (float): The magnitude of the effect used in the simulation.
        n_phase_bins (int): Number of bins used in the phase analysis.

    Returns:
        dict: A dictionary containing:
            - "H": Hit (1 if the effect was correctly detected, 0 otherwise).
            - "FA": False Alarm (currently initialized to 0).

    Copyright (C) 2026, Elio Balestrieri & Teresa Berther, University of Münster, Germany
    """

    center_bin = np.linspace(-np.pi, np.pi, n_phase_bins)

    trl_noise = 0
    phase_lag = 0
    n_trials = n_phase_bins
    p_resp = model_resp(n_trials, effect_mag, center_bin, phase_lag, trl_noise)
    # demean p_resp to avoid false signs
    p_resp = p_resp - p_resp.mean()

    H, FA = 0, 0
    for idx_effect, iclust_p in enumerate(summary["p"]):
        if iclust_p < 0.05:
            sign_theory = np.sign(p_resp[summary["idxs"][idx_effect]].sum())
            sign_clustmass = np.sign(summary["clustmass_stat"][idx_effect])

            if effect_mag:
                if sign_clustmass == sign_theory:
                    H = 1
            # EB: the following condition, while being theoretically justified, invariably inflates the FA rate based
            # on boundary points. In other words, the cluster does exist, but because of maybe extending one
            # point above what is considered the ground truth cluster -an arbitrary definition anyway- leads
            # to irrealistic FAR. Commented out, but with the memorandum that could be a valuable condition, if
            # better criteria for ground truth are defined.
            #    else:
            #        FA = 1
            else:
                FA = 1

    sdt_dict = {
        "H": H,
        "FA": FA,
    }

    return sdt_dict
