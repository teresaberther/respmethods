import numpy as np
import numba as nb
import tomllib as toml 
import matplotlib.pyplot as plt

@nb.vectorize
def nb_binom(n, p):
    return np.random.binomial(n, p)

@nb.njit
def observer(n_trials, effect_mag, phase_lag, trl_noise):
    
    phase_event = np.random.uniform(-np.pi, np.pi, n_trials)
    p_resp = (.5 + effect_mag*np.cos(phase_event + phase_lag) + 
              np.random.randn(n_trials)*trl_noise)
    ## maintain the presp within plausible range 
    # ...
    
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

@nb.njit
def pool_beh_resp(beh_events, bin_idxs):
    
    mv_avg = []
    for idxs in bin_idxs:
        
        binned_events = beh_events[idxs, 1]
        mv_avg.append(np.mean(binned_events))
        
    return np.array(mv_avg)

@nb.njit
def create_sample(n_subjects, intersubj_noise, n_trials, trl_noise,
                  effect_mag, range_bin, n_phase_bins, n_perms):
    
    center_bin = np.linspace(-np.pi, np.pi, n_phase_bins)
    range_bin *= np.pi # range_bin is a pi multiplier
    smpl_hr = np.empty((n_phase_bins, n_subjects))
    perms_hr = np.empty((n_phase_bins, n_subjects, n_perms))    
    
    for isubj in range(n_subjects):
        
        subj_lag = np.random.randn()*intersubj_noise
        exp_obs = observer(n_trials, effect_mag, subj_lag, trl_noise)
        binned = create_phase_bins(exp_obs[:, 0], center_bin, range_bin)
        hr_obs = pool_beh_resp(exp_obs, binned)
        smpl_hr[:,isubj] = hr_obs
        
        tmp_bins = exp_obs[:, 0]
        for iperm in range(n_perms):
            
            np.random.shuffle(tmp_bins)
            perm_binned = create_phase_bins(tmp_bins, center_bin, range_bin)
            ps_hr = pool_beh_resp(exp_obs, perm_binned)
            perms_hr[:, isubj, iperm] = ps_hr
        
    return smpl_hr, perms_hr, center_bin





def create_dict_from_toml(fname):
    
    cfg = {}
    ncombs_counter = 1
    with open(fname, 'rb') as f:
        fc = toml.load(f)        
        for cfg_name, cfg_values in fc.items():            
            if "fixed_range" in cfg_values.keys():            
                cfg.update({cfg_name : np.array(cfg_values["fixed_range"])})
            else:
                fr = np.linspace(cfg_values["strtval"],
                                 cfg_values["endval"],
                                 num=cfg_values["n_entries"])
                cfg.update({cfg_name : fr})
            ncombs_counter *= len(cfg[cfg_name])
        print(f"\n{ncombs_counter} combinations will be created. \n")
    return cfg


#%%  common

# ntrials = 500
# effect_mag = .2
# phase_lag = .2
# center_bin = np.linspace(-np.pi, np.pi, 51)
# range_bin = np.pi/10
# noise = .001


# plt.figure()
# plt.plot(center_bin, hr_obs)
