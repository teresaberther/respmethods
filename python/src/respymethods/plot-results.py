#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 15:34:32 2026

@author: ebalestr
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np

fname = "palma-217"

DF = pd.read_csv(f"{fname}.csv")
DFpr = DF.loc[DF["effect_mag"]>0, :]
DFabs = DF.loc[DF["effect_mag"]==0, :]


#%%

plt.figure()

plt.subplot(221)
sns.lineplot(data=DFpr, x="effect_mag", y="H", hue="n_trials")
plt.subplot(222)
sns.lineplot(data=DFpr, x="effect_mag", y="H", hue="n_subjects")
plt.subplot(223)
sns.lineplot(data=DFpr, x="effect_mag", y="H", hue="trl_noise")
plt.subplot(224)
sns.lineplot(data=DFpr, x="effect_mag", y="H", hue="intersubj_noise")

plt.suptitle("Prob reporting existing effect", fontsize=18)

#%%

fig = plt.figure()

acc=1
for i_intersubjnoise in np.sort(DFpr["intersubj_noise"].unique()):
    
    plt.subplot(2, 4, acc)
    tmpDF = DFpr.loc[DFpr["intersubj_noise"]==i_intersubjnoise, :]
    piv1 = tmpDF.pivot_table(index="n_trials", columns="n_subjects", 
                             values="H", aggfunc="mean")
    
    sns.heatmap(data=piv1, vmin=0, vmax=1, cbar=False, cmap="rocket_r")
    plt.title(f"intersubj noise\nexp({round(i_intersubjnoise, 3)})")
    acc+=1

# Colorbar in the 8th subplot
cbar_ax = plt.subplot(2, 4, 8)
norm = mcolors.Normalize(vmin=0, vmax=1)
cmap = cm.get_cmap("rocket_r")
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)

plt.suptitle("Prob reporting existing effect\n N subjects, N trials, and  inter-subject noise", 
             fontsize=18)

plt.tight_layout()

#%% 

fig = plt.figure()

acc=1
for i_intersubjnoise in np.sort(DFabs["intersubj_noise"].unique()):
    
    plt.subplot(2, 4, acc)
    tmpDF = DFabs.loc[DFabs["intersubj_noise"]==i_intersubjnoise, :]
    piv1 = tmpDF.pivot_table(index="n_trials", columns="n_subjects", 
                             values="FA", aggfunc="mean")
    
    sns.heatmap(data=piv1, vmin=0, vmax=1, cbar=False, cmap="rocket_r")
    plt.title(f"intersubj noise\nexp({round(i_intersubjnoise, 3)})")
    acc+=1

# Colorbar in the 8th subplot
cbar_ax = plt.subplot(2, 4, 8)
norm = mcolors.Normalize(vmin=0, vmax=.05)
cmap = cm.get_cmap("rocket_r")
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)

plt.suptitle("Prob reporting NON-existing effect\n N subjects, N trials, and  inter-subject noise", 
             fontsize=18)

plt.tight_layout()



#%%

    
    
# plt.subplot(222)
# sns.lineplot(data=DFpr, x="effect_mag", y="H", hue="n_subjects")
# plt.subplot(223)
# sns.lineplot(data=DFpr, x="effect_mag", y="H", hue="trl_noise")
# plt.subplot(224)
# sns.lineplot(data=DFpr, x="effect_mag", y="H", hue="intersubj_noise")

# plt.suptitle("Prob reporting existing effect", fontsize=18)



#%%

plt.figure()

plt.subplot(223)
sns.lineplot(data=DF, x="effect_mag", y="FA", hue="trl_noise")
plt.subplot(224)
sns.lineplot(data=DF, x="effect_mag", y="FA", hue="intersubj_noise")

#%%

plt.figure()

plt.subplot(121)
sns.lineplot(data=DF, x="effect_mag", y="H", hue="n_phase_bins")
plt.subplot(122)
sns.lineplot(data=DF, x="effect_mag", y="FA", hue="n_phase_bins")


