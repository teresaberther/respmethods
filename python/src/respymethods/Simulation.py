import numpy as np
import numba as nb

import matplotlib.pyplot as plt

@nb.vectorize
def nb_binom(n, p):
    return np.random.binomial(n, p)

@nb.njit
def observer(ntrials, effect_mag, phase_lag, noise):
    
    phase_event = np.random.uniform(-np.pi, np.pi, ntrials)
    p_resp = (.5 + effect_mag*np.cos(phase_event + phase_lag) + 
              np.random.randn(ntrials)*noise)
    resps = nb_binom(1, p_resp)
    resps = resps.astype(np.float64)
    
    beh_events = np.stack((phase_event, resps)).T
    
    return beh_events

@nb.njit
def create_phase_bins(phase_vect, center_bin, range_bin):
    
    x_vect = np.cos(phase_vect)
    y_vect = np.sin(phase_vect)    
    
    plus_min_range = np.exp(np.array([range_bin, -range_bin])*1j)
    cmplx_cntr = np.exp(center_bin*1j)

    bin_idxs = []
    idx_acc = 0
    for ibin in cmplx_cntr:
        
        x_ranges = np.real(ibin*plus_min_range) 
        y_ranges = np.imag(ibin*plus_min_range) 
                
        # correct for arc/axis intersection
        if (0>=x_ranges.min()) & (0<=x_ranges.max()):
            y_ranges = np.append(y_ranges, 1*np.sign(y_ranges[0]))
        if (0>=y_ranges.min()) & (0<=y_ranges.max()):
            x_ranges = np.append(x_ranges, 1*np.sign(x_ranges[0]))
        
        lgcl = ((x_vect>=x_ranges.min()) & (x_vect<x_ranges.max()) & 
                (y_vect>=y_ranges.min()) & (y_vect<y_ranges.max()))
        
        bin_idxs.append(np.where(lgcl)[0])
        idx_acc += 1
        
    return bin_idxs


def pool_beh_resp(beh_events, bin_idxs):
    
    mv_avg = []
    for idxs in bin_idxs:
        
        binned_events = beh_events[idxs, 1]
        mv_avg.append(np.mean(binned_events))
        
    return np.array(mv_avg)

def create_sample(nsubj, phase_lag_sd, noise):
    return

#%%  common

ntrials = 500
effect_mag = .2
phase_lag = .2
center_bin = np.linspace(-np.pi, np.pi, 51)
range_bin = np.pi/10
noise = .001

exp_obs = observer(ntrials, effect_mag, phase_lag, noise)
binned = create_phase_bins(exp_obs[:, 0], center_bin, range_bin)
hr_obs = pool_beh_resp(exp_obs, binned)

plt.figure()
plt.plot(center_bin, hr_obs)
