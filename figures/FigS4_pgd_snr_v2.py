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
from scipy.stats import binned_statistic_2d
import matplotlib as mpl
from matplotlib import cm

home = '/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW'

# Read in SNR files for the different noise levels
snr_10p_df = pd.read_csv(f'{home}/groundmotion_analysis/snr_files/snr_10p.csv')
snr_50p_df = pd.read_csv(f'{home}/groundmotion_analysis/snr_files/snr_50p.csv')
snr_90p_df = pd.read_csv(f'{home}/groundmotion_analysis/snr_files/snr_90p.csv')

# Get SNR
snr_10p = snr_10p_df.SNR.values
snr_50p = snr_50p_df.SNR.values
snr_90p = snr_90p_df.SNR.values

# Get magnitude
M = snr_10p_df.Mag.values

# Gather intensity measure files for the different noise levels
IM_files_10p = sorted(glob(f'{home}/groundmotion_analysis/IMs/cascadia/*_10p.csv'))
IM_files_50p = sorted(glob(f'{home}/groundmotion_analysis/IMs/cascadia/*_50p.csv'))
IM_files_90p = sorted(glob(f'{home}/groundmotion_analysis/IMs/cascadia/*_90p.csv'))

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

# # Read in waveforms for a smaller event
# small_raw = read(f'{home}/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_Signal/Cas-ONC-On-GNSS-Sig_cascadia-000022/Cas-ONC-On-GNSS-Sig_cascadia-000022_AL2H-LYE.mseed')
# small_10p = read(f'{home}/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_10p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000022/Cas-ONC-On-GNSS-SigNoise_cascadia-000022_AL2H-LYE.mseed')
# small_50p = read(f'{home}/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000022/Cas-ONC-On-GNSS-SigNoise_cascadia-000022_AL2H-LYE.mseed')
# small_90p = read(f'{home}/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_90p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000022/Cas-ONC-On-GNSS-SigNoise_cascadia-000022_AL2H-LYE.mseed')

# # Read in waveforms for a larger event
# large_raw = read(f'{home}/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_Signal/Cas-ONC-On-GNSS-Sig_cascadia-000083/Cas-ONC-On-GNSS-Sig_cascadia-000083_AL2H-LYE.mseed')
# large_10p = read(f'{home}/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_10p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000083/Cas-ONC-On-GNSS-SigNoise_cascadia-000083_AL2H-LYE.mseed')
# large_50p = read(f'{home}/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000083/Cas-ONC-On-GNSS-SigNoise_cascadia-000083_AL2H-LYE.mseed')
# large_90p = read(f'{home}/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_90p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000083/Cas-ONC-On-GNSS-SigNoise_cascadia-000083_AL2H-LYE.mseed')

# Get indicies for records within 400 km
dist_ind = np.where(rhyp_list<=400)[0]

#%% Make figure

cmap = plt.cm.seismic.copy()
cmap.set_bad(color='lightgray')  # Options: 'gray', 'lightgray', 'black', etc.

# Set up mosaic layout
layout = [
    ["A", "B"],
    ["null", "null"],
    ["C", "D"],
    ["E", "null2"]]

fig, axs = plt.subplot_mosaic(layout, figsize=(6.3,6.5), gridspec_kw={'height_ratios':[1.25,0.2,2,2]})

# # Smaller event timeseries with different noise levels
# axs['A'].plot(small_raw[0].times(),small_90p[0].data,lw=1,alpha=0.8,c='goldenrod',label='90th perc')
# axs['A'].plot(small_raw[0].times(),small_50p[0].data,lw=1,alpha=0.8,c='mediumturquoise',label='50th perc')
# axs['A'].plot(small_raw[0].times(),small_10p[0].data,lw=1,alpha=0.8,c='darkorchid',label='10th perc')
# axs['A'].plot(small_raw[0].times(),small_raw[0].data,lw=1,alpha=0.8,c='k',label='Signal')
# axs['A'].text(0.97, 0.95, 'M7.3', ha='right', va='top', transform=axs['A'].transAxes)
# axs['A'].set_xlim(0,600)
# axs['A'].set_ylim(-0.35,0.35)
# axs['A'].set_xlabel('Time (s)')
# axs['A'].xaxis.set_label_coords(1.1,-0.35)
# axs['A'].set_ylabel('Amplitude')
# handles, labels = axs['A'].get_legend_handles_labels()
# handles = [handles[i] for i in [3,2,1,0]]
# labels = [labels[i] for i in [3,2,1,0]]
# axs['A'].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.11,-0.6),facecolor='white',
#                 frameon=True,fontsize=9,ncol=4,markerscale=2)
# axs['A'].text(-0.35, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs['A'].transAxes)

# # Larger event timeseries with different noise 
# axs['B'].plot(small_raw[0].times(),large_90p[0].data,lw=1,alpha=0.8,c='goldenrod',label='90p Noise')
# axs['B'].plot(small_raw[0].times(),large_50p[0].data,lw=1,alpha=0.8,c='mediumturquoise',label='50p Noise')
# axs['B'].plot(small_raw[0].times(),large_10p[0].data,lw=1,alpha=0.8,c='darkorchid',label='10p Noise')
# axs['B'].plot(small_raw[0].times(),large_raw[0].data,lw=1,alpha=0.8,c='k',label='Signal')
# axs['B'].text(0.97, 0.95, 'M9.0', ha='right', va='top', transform=axs['B'].transAxes)
# axs['B'].set_xlim(0,600)
# axs['B'].text(-0.2, 1, r'$\bf{(b)}$', ha='right', va='top', transform=axs['B'].transAxes)

x_bins = np.linspace(6.5, 9.5, 10)
y_bins = np.logspace(np.log10(0.4), np.log10(700), 10)

# 10p noise
# statistic, xedges, yedges, binnumber = binned_statistic_2d(M, snr_10p, pgd_res_gold_10p, statistic='median', bins=[x_bins, y_bins])
# masked_statistic = np.ma.masked_invalid(statistic.T)
# X, Y = np.meshgrid(xedges, yedges)
# pcm10 = axs['C'].pcolormesh(X, Y, masked_statistic, vmin=-4.75, vmax=4.75, cmap=cmap, shading='auto') 
# axs['C'].set_yscale('log')
# axs['C'].set_xlabel('Magnitude')
# axs['C'].set_ylabel('SNR')

axs['C'].scatter(M, snr_10p, marker='o', fc='none', ec=pgd_res_gold_10p, cmap=cmap, vmin=-4.75, vmax=4.75, alpha=0.65, lw=0.4)
axs['C'].set_yscale('log')

# 50p noise
statistic, xedges, yedges, binnumber = binned_statistic_2d(M, snr_50p, pgd_res_gold_50p, statistic='median', bins=[x_bins, y_bins])
masked_statistic = np.ma.masked_invalid(statistic.T)
X, Y = np.meshgrid(xedges, yedges)
pcm50 = axs['D'].pcolormesh(X, Y, masked_statistic, vmin=-4.75, vmax=4.75, cmap=cmap, shading='auto') 
axs['D'].set_yscale('log')
axs['D'].set_xlabel('Magnitude')
axs['D'].set_ylabel('SNR')

# 90p noise
statistic, xedges, yedges, binnumber = binned_statistic_2d(M, snr_90p, pgd_res_gold_90p, statistic='median', bins=[x_bins, y_bins])
masked_statistic = np.ma.masked_invalid(statistic.T)
X, Y = np.meshgrid(xedges, yedges)
pcm90 = axs['E'].pcolormesh(X, Y, masked_statistic, vmin=-4.75, vmax=4.75, cmap=cmap, shading='auto') 
axs['E'].set_yscale('log')
axs['E'].set_xlabel('Magnitude')
axs['E'].set_ylabel('SNR')

axs['null'].remove()
axs['null2'].remove()

# Colorbar
cb_ax = fig.add_axes([0.65, 0.35, 0.31, 0.02])  # [left, bottom, width, height] in figure fraction

fig.colorbar(
    pcm90,
    cax=cb_ax,
    orientation='horizontal',
    label=r'$\delta_{PGD,GA21}$'
)

plt.subplots_adjust(hspace=0.5,wspace=0.35,left=0.165,right=0.975,top=0.975,bottom=0.175)
# plt.savefig(f'{home}/manuscript/figures/Fig8_PGD_SNR.png',dpi=300)
