#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 15:34:32 2026

@author: ebalestr
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rc('axes', titlesize=18)
plt.rc('axes', labelsize=15)
plt.rcParams.update({'legend.fontsize': 12})

fname = "palma-217"

DF = pd.read_csv(f"{fname}.csv")


#%%

plt.figure()

plt.subplot(221)
sns.lineplot(data=DF, x="effect_mag", y="H", hue="n_trials")
plt.subplot(222)
sns.lineplot(data=DF, x="effect_mag", y="H", hue="n_subjects")
plt.subplot(223)
sns.lineplot(data=DF, x="effect_mag", y="FA", hue="n_trials")
plt.subplot(224)
sns.lineplot(data=DF, x="effect_mag", y="FA", hue="n_subjects")

#%%

plt.figure()

plt.subplot(221)
sns.lineplot(data=DF, x="effect_mag", y="H", hue="trl_noise")
plt.subplot(222)
sns.lineplot(data=DF, x="effect_mag", y="H", hue="intersubj_noise")
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


