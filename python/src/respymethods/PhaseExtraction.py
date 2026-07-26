import numpy as np
from numba import njit
from scipy.signal import find_peaks


def two_point_interp(ts):
    # SciPy call stays in pure Python (fast enough, and fully supported)
    peak_idxs, _ = find_peaks(ts, prominence=1)
    peak_idxs = peak_idxs.astype(np.int64)
    return _two_point_interp_core(ts, peak_idxs)


@njit(cache=True)
def _two_point_interp_core(ts, peak_idxs):
    troughs = np.empty(peak_idxs.shape[0]-1, dtype=np.int64)
    for k in range(1, peak_idxs.shape[0]):
        this_idx = peak_idxs[k]
        prev_idx = peak_idxs[k - 1]
        tmp = ts[prev_idx:this_idx]
        troughs[k - 1] = prev_idx + np.argmin(tmp)

    phasevec = np.nan * np.ones(ts.shape)
    phasevec[peak_idxs] = 0.0
    phasevec[troughs] = np.pi

    for k in range(peak_idxs.shape[0] - 1):
        this_peak_idx = peak_idxs[k]
        n_p2t = troughs[k] - this_peak_idx
        n_t2p = peak_idxs[k + 1] - troughs[k]

        s1 = np.linspace(np.pi / n_p2t, np.pi, n_p2t)
        s2 = np.linspace(-np.pi + np.pi / n_t2p, 0.0, n_t2p)

        phasevec[this_peak_idx : troughs[k]] = s1
        phasevec[troughs[k] : peak_idxs[k + 1]] = s2

    return phasevec


# def two_point_interp(ts):
#
#    # Initialize variables for types Numba expects
#    peak_idxs = np.empty(0, dtype=np.int64)
#
#    # Temporarily switch to object mode to run SciPy
#    with objmode(peak_idxs="int64[:]"):
#        peak_idxs, _ = find_peaks(ts, prominence=1)
#    troughs = np.ones(peak_idxs.shape, dtype=np.int64)
#    for k, this_idx in enumerate(peak_idxs[1:], start=1):
#        tmp = ts[peak_idxs[k - 1] : this_idx]
#        tmp_idxs = np.arange(peak_idxs[k - 1], this_idx)
#        troughs[k - 1] = tmp_idxs[np.argmin(tmp)]
#
#    phasevec = np.nan * np.ones(ts.shape)
#    phasevec[peak_idxs] = 0
#    phasevec[troughs] = np.pi
#
#    for k, this_peak_idx in enumerate(peak_idxs[:-1]):
#        n_p2t = troughs[k] - this_peak_idx  # n points between peak1 and trough1
#        n_t2p = peak_idxs[k + 1] - troughs[k]  # n points between trough1 and peak2
#        s1 = np.linspace(
#            0 + np.pi / n_p2t, np.pi, n_p2t
#        )  # linear interpolation peak2trough
#        s2 = np.linspace(
#            -np.pi + np.pi / n_t2p, 0, n_t2p
#        )  # linear interpolation trough2peak
#        phasevec[this_peak_idx : troughs[k]] = (
#            s1  # substitute with phase angles (peak2trough)
#        )
#        phasevec[troughs[k] : peak_idxs[k + 1]] = (
#            s2  # substitute with phase angles (peak2trough)
#        )
#    return phasevec
#
