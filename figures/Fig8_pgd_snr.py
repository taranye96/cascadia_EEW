#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 13:05:09 2023

@author: tnye
"""

###############################################################################
# This script makes Figure 8, which shows the relationships between noise
# levels, SNR, and PGD.
###############################################################################

# Imports
import numpy as np
import pandas as pd
from obspy import read
from glob import glob
import matplotlib.pyplot as plt

# Read in SNR files for the different noise levels
snr_10p_df = pd.read_csv('/Users/tnye/ONC/groundmotion_analysis/snr_files/snr_10p.csv')
snr_50p_df = pd.read_csv('/Users/tnye/ONC/groundmotion_analysis/snr_files/snr_50p.csv')
snr_90p_df = pd.read_csv('/Users/tnye/ONC/groundmotion_analysis/snr_files/snr_90p.csv')

# Get SNR
snr_10p = snr_10p_df.SNR.values
snr_50p = snr_50p_df.SNR.values
snr_90p = snr_90p_df.SNR.values

# Get magnitude
M = snr_10p_df.Mag.values

# Gather intensity measure files for the different noise levels
IM_files_10p = sorted(glob('/Users/tnye/ONC/groundmotion_analysis/IMs/cascadia/*_10p.csv'))
IM_files_50p = sorted(glob('/Users/tnye/ONC/groundmotion_analysis/IMs/cascadia/*_50p.csv'))
IM_files_90p = sorted(glob('/Users/tnye/ONC/groundmotion_analysis/IMs/cascadia/*_90p.csv'))

# Initialize lists
pgd_res_melgar_90p = np.array([])
pgd_res_melgar_50p = np.array([])
pgd_res_melgar_10p = np.array([])
pgd_res_gold_90p = np.array([])
pgd_res_gold_50p = np.array([])
pgd_res_gold_10p = np.array([])
rhyp_list = np.array([])
rrup_list = np.array([])
events = np.array([])
stns = np.array([])

# Loop over records in intensity measure files
for i in range(len(IM_files_10p)):
    
    # Get PGD residuals
    pgd_res_melgar_90p = np.append(pgd_res_melgar_90p, pd.read_csv(IM_files_90p[i])['lnPGD_noise_res_MA15'].values)
    pgd_res_melgar_50p = np.append(pgd_res_melgar_50p, pd.read_csv(IM_files_50p[i])['lnPGD_noise_res_MA15'].values)
    pgd_res_melgar_10p = np.append(pgd_res_melgar_10p, pd.read_csv(IM_files_10p[i])['lnPGD_noise_res_MA15'].values)
    pgd_res_gold_90p = np.append(pgd_res_gold_90p, pd.read_csv(IM_files_90p[i])['lnPGD_noise_res_GA21'].values)
    pgd_res_gold_50p = np.append(pgd_res_gold_50p, pd.read_csv(IM_files_50p[i])['lnPGD_noise_res_GA21'].values)
    pgd_res_gold_10p = np.append(pgd_res_gold_10p, pd.read_csv(IM_files_10p[i])['lnPGD_noise_res_GA21'].values)
    
    # Get distances
    rhyp_list = np.append(rhyp_list,pd.read_csv(IM_files_10p[i])['Rhyp(km)'])
    rrup_list = np.append(rrup_list,pd.read_csv(IM_files_10p[i])['Rrup(km)'])
    
    # Get run and station names
    events = np.append(events,[IM_files_10p[i].split('/')[-1].split('_')[0]]*len(pd.read_csv(IM_files_10p[i])))
    stns = np.append(stns,pd.read_csv(IM_files_10p[i])['Station'])

# Read in waveforms for a smaller event
small_raw = read('/Users/tnye/ONC/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_Signal/Cas-ONC-On-GNSS-Sig_cascadia-000022/Cas-ONC-On-GNSS-Sig_cascadia-000022_AL2H-LYE.mseed')
small_10p = read('/Users/tnye/ONC/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_10p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000022/Cas-ONC-On-GNSS-SigNoise_cascadia-000022_AL2H-LYE.mseed')
small_50p = read('/Users/tnye/ONC/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000022/Cas-ONC-On-GNSS-SigNoise_cascadia-000022_AL2H-LYE.mseed')
small_90p = read('/Users/tnye/ONC/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_90p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000022/Cas-ONC-On-GNSS-SigNoise_cascadia-000022_AL2H-LYE.mseed')

# Read in waveforms for a larger event
large_raw = read('/Users/tnye/ONC/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_Signal/Cas-ONC-On-GNSS-Sig_cascadia-000083/Cas-ONC-On-GNSS-Sig_cascadia-000083_AL2H-LYE.mseed')
large_10p = read('/Users/tnye/ONC/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_10p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000083/Cas-ONC-On-GNSS-SigNoise_cascadia-000083_AL2H-LYE.mseed')
large_50p = read('/Users/tnye/ONC/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000083/Cas-ONC-On-GNSS-SigNoise_cascadia-000083_AL2H-LYE.mseed')
large_90p = read('/Users/tnye/ONC/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_90p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000083/Cas-ONC-On-GNSS-SigNoise_cascadia-000083_AL2H-LYE.mseed')

# Get indicies for records within 400 km
dist_ind = np.where(rhyp_list<=400)[0]

#%% Make figure

# Set up mosaic layout
layout = [
    ["A", "B"],
    ["null", "null"],
    ["C", "C"],
    ["D", "D"]]

fig, axs = plt.subplot_mosaic(layout, figsize=(6.3,6), gridspec_kw={'height_ratios':[1.25,0.2,2,2]})

# Smaller event timeseries with different noise levels
axs['A'].plot(small_raw[0].times(),small_90p[0].data,lw=1,alpha=0.8,c='goldenrod',label='90th perc')
axs['A'].plot(small_raw[0].times(),small_50p[0].data,lw=1,alpha=0.8,c='mediumturquoise',label='50th perc')
axs['A'].plot(small_raw[0].times(),small_10p[0].data,lw=1,alpha=0.8,c='darkorchid',label='10th perc')
axs['A'].plot(small_raw[0].times(),small_raw[0].data,lw=1,alpha=0.8,c='k',label='Signal')
axs['A'].text(0.97, 0.95, 'M7.3', ha='right', va='top', transform=axs['A'].transAxes)
axs['A'].set_xlim(0,600)
axs['A'].set_ylim(-0.35,0.35)
axs['A'].set_xlabel('Time (s)')
axs['A'].xaxis.set_label_coords(1.1,-0.35)
axs['A'].set_ylabel('Amplitude')
handles, labels = axs['A'].get_legend_handles_labels()
handles = [handles[i] for i in [3,2,1,0]]
labels = [labels[i] for i in [3,2,1,0]]
axs['A'].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.11,-0.6),facecolor='white',
                frameon=True,fontsize=9,ncol=4,markerscale=2)
axs['A'].text(-0.35, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs['A'].transAxes)

# Larger event timeseries with different noise 
axs['B'].plot(small_raw[0].times(),large_90p[0].data,lw=1,alpha=0.8,c='goldenrod',label='90p Noise')
axs['B'].plot(small_raw[0].times(),large_50p[0].data,lw=1,alpha=0.8,c='mediumturquoise',label='50p Noise')
axs['B'].plot(small_raw[0].times(),large_10p[0].data,lw=1,alpha=0.8,c='darkorchid',label='10p Noise')
axs['B'].plot(small_raw[0].times(),large_raw[0].data,lw=1,alpha=0.8,c='k',label='Signal')
axs['B'].text(0.97, 0.95, 'M9.0', ha='right', va='top', transform=axs['B'].transAxes)
axs['B'].set_xlim(0,600)
axs['B'].text(-0.2, 1, r'$\bf{(b)}$', ha='right', va='top', transform=axs['B'].transAxes)

# Magnitude vs SNR
axs['null'].remove()
axs['C'].scatter(snr_10p[dist_ind], M[dist_ind], s=10, marker='o', lw=0.4, facecolors='none', edgecolors='darkorchid', alpha=0.6, label='10p Noise')
axs['C'].scatter(snr_50p[dist_ind], M[dist_ind], s=10, marker='^', lw=0.4, facecolors='none', edgecolors='mediumturquoise', alpha=0.6, label='50p Noise')
axs['C'].scatter(snr_90p[dist_ind], M[dist_ind], s=10, marker='s', lw=0.4, facecolors='none', edgecolors='goldenrod', alpha=0.6, label='90p Noise')
axs['C'].set_xscale('log')
axs['C'].set_ylabel('Magnitude')
axs['C'].grid(alpha=0.25)
axs['C'].text(-0.15, 1, r'$\bf{(c)}$', ha='right', va='top', transform=axs['C'].transAxes)

# Magnitude vs PGD residual
axs['D'].scatter(snr_10p[dist_ind], pgd_res_gold_10p[dist_ind], s=10, marker='o', lw=0.4, facecolors='none', edgecolors='darkorchid', alpha=0.6, label='10th perc')
axs['D'].scatter(snr_50p[dist_ind], pgd_res_gold_50p[dist_ind], s=10, marker='^', lw=0.4, facecolors='none', edgecolors='mediumturquoise', alpha=0.6, label='50th perc')
axs['D'].scatter(snr_90p[dist_ind], pgd_res_gold_90p[dist_ind], s=10, marker='s', lw=0.4, facecolors='none', edgecolors='goldenrod', alpha=0.6, label='90th perc')
axs['D'].axhline(np.log(0.255),lw=0.7,ls='--',c='k')
axs['D'].axhline(-np.log(0.255),lw=0.7,ls='--',c='k')
axs['D'].axhline(0, ls='-', c='k', lw=1)
axs['D'].set_xscale('log')
# axs['D'].set_xlim(1,1300)
axs['D'].set_xlabel('SNR')
axs['D'].set_ylabel(r'$\delta_{PGD,GA21}$')
axs['D'].grid(alpha=0.25)
axs['D'].legend(loc='upper center', bbox_to_anchor=(0.5,-0.4),facecolor='white',
                frameon=True,fontsize=9,ncol=4,markerscale=2)
axs['D'].text(-0.15, 1, r'$\bf{(d)}$', ha='right', va='top', transform=axs['D'].transAxes)

plt.subplots_adjust(hspace=0.5,wspace=0.35,left=0.165,right=0.975,top=0.975,bottom=0.175)
plt.savefig('/Users/tnye/ONC/manuscript/figures/Fig8_PGD_SNR.png',dpi=300)
