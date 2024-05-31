#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 09:47:00 2023

@author: tnye
"""

# Imports
import numpy as np
import pandas as pd
from obspy import read
from mtspec import mtspec
import tsueqs_main_fns as tmf
import matplotlib.pyplot as plt
import seaborn as sns

bb_file = '/Users/tnye/ONC/simulations/cascadia/raw_waveforms/ONC-strong-motion/eews.000095/AL2H.bb.HNE.sac'
hf_file = '/Users/tnye/ONC/simulations/cascadia/AL2H.HNE.mpi.sac'
lf_file = '/Users/tnye/ONC/simulations/cascadia/raw_waveforms/ONC-GNSS/eews.000095/AL2H.LYE.sac'

bb_st = read(bb_file)
lf_st = read(lf_file)
hf_st = read(hf_file)

lf_acc = np.diff(np.diff(lf_st[0].data))
hf_data = hf_st[0].data
bb_data = bb_st[0].data

lf_filt = np.diff(np.diff(tmf.lowpass(lf_st,0.998, lf_st[0].stats.sampling_rate, 4)[0].data))
hf_filt = tmf.highpass(hf_st, 0.01, bb_st[0].stats.sampling_rate, 4)[0].data

lf_amp2, freq_lf = mtspec(lf_acc, lf_st[0].stats.delta, time_bandwidth=4, 
                          number_of_tapers=5, nfft=lf_st[0].stats.npts, quadratic=True)
hf_amp2, freq_hf = mtspec(hf_st[0].data, hf_st[0].stats.delta, time_bandwidth=4, 
                          number_of_tapers=5, nfft=hf_st[0].stats.npts, quadratic=True)
bb_amp2, freq_hf = mtspec(bb_st[0].data, bb_st[0].stats.delta, time_bandwidth=4, 
                          number_of_tapers=5, nfft=bb_st[0].stats.npts, quadratic=True)

lf_amp = np.sqrt(lf_amp2)
hf_amp = np.sqrt(hf_amp2)
bb_amp = np.sqrt(bb_amp2)

lf_amp2_filt, freq_lf = mtspec(lf_filt, lf_st[0].stats.delta, time_bandwidth=4, 
                          number_of_tapers=5, nfft=lf_st[0].stats.npts, quadratic=True)
hf_amp2_filt, freq_hf = mtspec(hf_filt, hf_st[0].stats.delta, time_bandwidth=4, 
                          number_of_tapers=5, nfft=hf_st[0].stats.npts, quadratic=True)

lf_amp_filt = np.sqrt(lf_amp2_filt)
hf_amp_filt = np.sqrt(hf_amp2_filt)


plt.subplots(1,1,figsize=(5,1.5))
plt.plot(lf_st[0].times()[:-2],lf_acc,c='k',lw=1,label='AL2H')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()
plt.subplots_adjust(bottom=0.3, right=0.95, left=0.15)
plt.savefig('/Users/tnye/ONC/plots/misc/example_wfs/lf_acc.png',dpi=300)

plt.subplots(1,1,figsize=(5,1.5))
plt.plot(hf_st[0].times(),hf_st[0].data,c='k',lw=1,label='AL2H')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()
plt.subplots_adjust(bottom=0.3, right=0.95, left=0.15)
plt.savefig('/Users/tnye/ONC/plots/misc/example_wfs/hf.png',dpi=300)

plt.subplots(1,1,figsize=(5,1.5))
plt.plot(bb_st[0].times(),bb_st[0].data,c='k',lw=1,label='AL2H')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()
plt.subplots_adjust(bottom=0.3, right=0.95, left=0.15)
plt.savefig('/Users/tnye/ONC/plots/misc/example_wfs/bb.png',dpi=300)

# spectra
sns.set_style('darkgrid')
plt.subplots(1,1,figsize=(4,2.75))
plt.loglog(freq_lf, lf_amp_filt, alpha=0.7, lw=1, label='low frequency')
plt.loglog(freq_hf, hf_amp_filt, alpha=0.7, lw=1, label='high frequency')
plt.loglog(freq_hf, bb_amp, lw=1, alpha=0.7, label='full spectrum')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.legend()
plt.subplots_adjust(bottom=0.2, right=0.95, left=0.15)
plt.savefig('/Users/tnye/ONC/plots/misc/mathced_filter.png',dpi=300)






sns.set_style('white')

st = read('/Users/tnye/ONC/simulations/cascadia/raw_waveforms/ONC-strong-motion/eews.000095/AL2H.bb.HNE.sac')
fig, ax = plt.subplots(1,1,figsize=(5.5,1.5))
ax.plot(st[0].times(),st[0].data,c='k',lw=1)
ax.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
plt.subplots_adjust(bottom=0.3, right=0.95, left=0.05)
plt.xticks([]) 
plt.yticks([]) 
plt.savefig('/Users/tnye/ONC/plots/misc/example_wfs/misc1.png',dpi=300)



st = read('/Users/tnye/ONC/simulations/cascadia/raw_waveforms/ONC-strong-motion/eews.000065/FAHB.bb.HNE.sac')
fig, ax = plt.subplots(1,1,figsize=(5.5,1.5))
ax.plot(st[0].times(),st[0].data,c='k',lw=1)
ax.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
plt.subplots_adjust(bottom=0.3, right=0.95, left=0.05)
plt.xticks([]) 
plt.yticks([]) 
plt.savefig('/Users/tnye/ONC/plots/misc/example_wfs/misc2.png',dpi=300)



st = read('/Users/tnye/ONC/simulations/cascadia/raw_waveforms/ONC-strong-motion/eews.000100/QUAD.bb.HNE.sac')
fig, ax = plt.subplots(1,1,figsize=(5.5,1.5))
ax.plot(st[0].times(),st[0].data,c='k',lw=1)
ax.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
plt.subplots_adjust(bottom=0.3, right=0.95, left=0.05)
plt.xticks([]) 
plt.yticks([]) 
plt.savefig('/Users/tnye/ONC/plots/misc/example_wfs/misc3.png',dpi=300)








