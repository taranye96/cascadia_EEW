#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 09:35:25 2023

@author: tnye
"""

# Imports
import numpy as np
import pandas as pd
from obspy import read
from glob import glob
from mtspec import mtspec
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, MultipleLocator, ScalarFormatter, FormatStrFormatter
import signal_average_fns as avg
import IM_fns


velmod = np.genfromtxt('/Users/tnye/ONC/files/cascadia.mod')


d = []

for i in range(len(velmod)):
    if i == 0:
        d.append(0)
    else:
        d.append(np.sum(velmod[:i,0])) 



fig, ax = plt.subplots(1,1, figsize=(3,6))

ax.plot(velmod[:,2],d,c='mediumpurple',lw=1,label='Vp')
ax.plot(velmod[:,3],d,c='goldenrod',lw=1,label='Vs')
ax.grid(alpha=0.25)
ax.invert_yaxis()
# ax.set_xlim(xmax=3.5)
# ax.set_ylim(1,0)
ax.set_xlabel('V (km/s)')
ax.set_ylabel('Depth (km)')
ax.legend()

plt.subplots_adjust(top=0.95, bottom=0.1, right=0.95, left=0.2, wspace=0.35, hspace=0.3)

plt.savefig('/Users/tnye/ONC/manuscript/figures/S_velmod.png',dpi=300)
