import numpy as np
import numba as nb
import tomllib as toml 
from scipy.stats import zscore

# pack configs combinations
from itertools import product

# visualize
import matplotlib.pyplot as plt

# respmethods
import Simulation as sim
from RespStats import circ_perm


#%% load toml
fc = sim.create_dict_from_toml("sim-cfg.toml")

# generation of all the parameters combination
combinations = product(*fc.values())
results = [dict(zip(fc.keys(), c)) for c in combinations]

tmp_dict = results[0]
smpl, perms, theta = sim.create_sample(**tmp_dict)
aggr_dat = np.concatenate((perms, smpl[:, :, np.newaxis]), axis=2)
tmp_z = zscore(aggr_dat, axis=2)

permemp_t = np.sqrt(tmp_dict["n_subjects"])*tmp_z.mean(axis=1)/tmp_z.std(axis=1)

emp_t = permemp_t[:, -1]
perm_t = permemp_t[:, 0:-1]

[bounds, summary] = circ_perm(emp_t, perm_t, theta)


#%%

gavg = smpl.mean(axis=1)
plt.figure()
plt.plot(theta, emp_t)

#%%
out = sim.sdt(summary, tmp_dict["effect_mag"], tmp_dict["n_phase_bins"])

