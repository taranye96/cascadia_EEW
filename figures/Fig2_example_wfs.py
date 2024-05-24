#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 17:28:07 2023

@author: tnye
"""

# Imports
import numpy as np
from obspy import read
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt

lf_st = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/output/waveforms/cascadia.000050/AL2H.LYE.sac')
hf_st = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/output/waveforms/cascadia.000050/AL2H.HNE.mpi.sac')
bb_st = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/output/waveforms/cascadia.000050/AL2H.bb.HNE.sac')

# Make figure
fig, axs = plt.subplots(3,1, figsize=(4,3))
axs[0].plot(lf_st[0].times()[200:400],lf_st[0].data[200:400],c='k',lw=0.8)
axs[0].set_ylabel('Disp (m)')
axs[0].text(0.98,0.95,'Low frequency',transform=axs[0].transAxes,va='top',ha='right')
axs[0].xaxis.set_major_locator(ticker.NullLocator())
axs[0].yaxis.set_major_locator(MultipleLocator(0.02))
axs[1].plot(hf_st[0].times()[10000:20000],hf_st[0].data[10000:20000],c='k',lw=0.8)
axs[1].text(0.98,0.95,'High frequency',transform=axs[1].transAxes,va='top',ha='right')
axs[1].xaxis.set_major_locator(ticker.NullLocator())
axs[1].set_ylabel(r'Acc (m/s$^2$)')
axs[1].yaxis.set_major_locator(MultipleLocator(0.01))
axs[1].set_ylim(-0.0125,0.0125)
axs[2].plot(bb_st[0].times()[10000:20000],bb_st[0].data[10000:20000],c='k',lw=0.8)
axs[2].text(0.98,0.95,'Broadband',transform=axs[2].transAxes,va='top',ha='right')
axs[2].set_xticks([100,120,140,160,180,200])
axs[2].set_xticklabels([0,20,40,60,80,100])
axs[2].yaxis.set_major_locator(MultipleLocator(0.01))
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel(r'Acc (m/s$^2$)')
axs[2].set_ylim(-0.0125,0.0125)
plt.subplots_adjust(top=0.96, bottom=0.15, right=0.98, left=0.21, hspace=0.1)
plt.savefig('/Users/tnye/ONC/manuscript/figures/unannotated/example_wfs.png',dpi=300)