#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 13:05:09 2023

@author: tnye
"""

###############################################################################
# Script used to make the PGD residual SNR analysis figure.
###############################################################################

# Imports
from os import path, makedirs
import numpy as np
import pandas as pd
from obspy import read
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

snr_10p_df = pd.read_csv('/Users/tnye/ONC/groundmotion_analysis/snr_files/snr_10p.csv')
snr_50p_df = pd.read_csv('/Users/tnye/ONC/groundmotion_analysis/snr_files/snr_50p.csv')
snr_90p_df = pd.read_csv('/Users/tnye/ONC/groundmotion_analysis/snr_files/snr_90p.csv')

M = snr_10p_df.Mag.values
snr_10p = snr_10p_df.SNR.values
snr_50p = snr_50p_df.SNR.values
snr_90p = snr_90p_df.SNR.values

IM_files_10p = sorted(glob('/Users/tnye/ONC/flatfiles/IMs/cascadia_longer_wfs/*_10p.csv'))
IM_files_50p = sorted(glob('/Users/tnye/ONC/flatfiles/IMs/cascadia_longer_wfs/*_50p.csv'))
IM_files_90p = sorted(glob('/Users/tnye/ONC/flatfiles/IMs/cascadia_longer_wfs/*_90p.csv'))

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

for i in range(len(IM_files_10p)):
    pgd_res_melgar_90p = np.append(pgd_res_melgar_90p, pd.read_csv(IM_files_90p[i])['lnPGD_noise_res_MA15'].values)
    pgd_res_melgar_50p = np.append(pgd_res_melgar_50p, pd.read_csv(IM_files_50p[i])['lnPGD_noise_res_MA15'].values)
    pgd_res_melgar_10p = np.append(pgd_res_melgar_10p, pd.read_csv(IM_files_10p[i])['lnPGD_noise_res_MA15'].values)
    pgd_res_gold_90p = np.append(pgd_res_gold_90p, pd.read_csv(IM_files_90p[i])['lnPGD_noise_res_GA21'].values)
    pgd_res_gold_50p = np.append(pgd_res_gold_50p, pd.read_csv(IM_files_50p[i])['lnPGD_noise_res_GA21'].values)
    pgd_res_gold_10p = np.append(pgd_res_gold_10p, pd.read_csv(IM_files_10p[i])['lnPGD_noise_res_GA21'].values)
    rhyp_list = np.append(rhyp_list,pd.read_csv(IM_files_10p[i])['Rhyp(km)'])
    rrup_list = np.append(rrup_list,pd.read_csv(IM_files_10p[i])['Rrup(km)'])
    events = np.append(events,[IM_files_10p[i].split('/')[-1].split('_')[0]]*len(pd.read_csv(IM_files_10p[i])))
    stns = np.append(stns,pd.read_csv(IM_files_10p[i])['Station'])

small_raw = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_Signal/Cas-ONC-On-GNSS-Sig_cascadia-000022/Cas-ONC-On-GNSS-Sig_cascadia-000022_AL2H-LYE.mseed')
small_10p = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-ONC-Onshore-GNSS_10p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000022/Cas-ONC-On-GNSS-SigNoise_cascadia-000022_AL2H-LYE.mseed')
small_50p = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000022/Cas-ONC-On-GNSS-SigNoise_cascadia-000022_AL2H-LYE.mseed')
small_90p = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-ONC-Onshore-GNSS_90p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000022/Cas-ONC-On-GNSS-SigNoise_cascadia-000022_AL2H-LYE.mseed')

large_raw = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_Signal/Cas-ONC-On-GNSS-Sig_cascadia-000083/Cas-ONC-On-GNSS-Sig_cascadia-000083_AL2H-LYE.mseed')
large_10p = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-ONC-Onshore-GNSS_10p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000083/Cas-ONC-On-GNSS-SigNoise_cascadia-000083_AL2H-LYE.mseed')
large_50p = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-ONC-Onshore-GNSS_50p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000083/Cas-ONC-On-GNSS-SigNoise_cascadia-000083_AL2H-LYE.mseed')
large_90p = read('/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-ONC-Onshore-GNSS_90p/Cas-ONC-On-GNSS_SignalwithNoise/Cas-ONC-On-GNSS-SigNoise_cascadia-000083/Cas-ONC-On-GNSS-SigNoise_cascadia-000083_AL2H-LYE.mseed')

dist_ind = np.where(rhyp_list<=400)[0]

#%%

# mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['font.family'] = 'Helvetica'

layout = [
    ["A", "B"],
    ["null", "null"],
    ["C", "C"],
    ["D", "D"]]

fig, axs = plt.subplot_mosaic(layout, figsize=(6.3,6), gridspec_kw={'height_ratios':[1.25,0.2,2,2]})

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
# axs['A'].set_title('cascadia.000022, AL2H.HNE', fontsize=11)
handles, labels = axs['A'].get_legend_handles_labels()
handles = [handles[i] for i in [3,2,1,0]]
labels = [labels[i] for i in [3,2,1,0]]
axs['A'].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.11,-0.6),facecolor='white',
                frameon=True,fontsize=9,ncol=4,markerscale=2)
axs['A'].text(-0.35, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs['A'].transAxes)

axs['B'].plot(small_raw[0].times(),large_90p[0].data,lw=1,alpha=0.8,c='goldenrod',label='90p Noise')
axs['B'].plot(small_raw[0].times(),large_50p[0].data,lw=1,alpha=0.8,c='mediumturquoise',label='50p Noise')
axs['B'].plot(small_raw[0].times(),large_10p[0].data,lw=1,alpha=0.8,c='darkorchid',label='10p Noise')
axs['B'].plot(small_raw[0].times(),large_raw[0].data,lw=1,alpha=0.8,c='k',label='Signal')
axs['B'].text(0.97, 0.95, 'M9.0', ha='right', va='top', transform=axs['B'].transAxes)
axs['B'].set_xlim(0,600)
# axs['B'].set_ylim(-0.35,0.35)
axs['B'].text(-0.2, 1, r'$\bf{(b)}$', ha='right', va='top', transform=axs['B'].transAxes)

axs['null'].remove()
axs['C'].scatter(snr_10p[dist_ind], M[dist_ind], s=10, marker='o', lw=0.4, facecolors='none', edgecolors='darkorchid', alpha=0.6, label='10p Noise')
axs['C'].scatter(snr_50p[dist_ind], M[dist_ind], s=10, marker='^', lw=0.4, facecolors='none', edgecolors='mediumturquoise', alpha=0.6, label='50p Noise')
axs['C'].scatter(snr_90p[dist_ind], M[dist_ind], s=10, marker='s', lw=0.4, facecolors='none', edgecolors='goldenrod', alpha=0.6, label='90p Noise')
axs['C'].set_xscale('log')
# axs['C'].set_xlim(1,1300)
axs['C'].set_ylabel('Magnitude')
axs['C'].grid(alpha=0.25)
axs['C'].text(-0.15, 1, r'$\bf{(c)}$', ha='right', va='top', transform=axs['C'].transAxes)

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
plt.savefig('/Users/tnye/ONC/manuscript/figures/PGD_SNR.png',dpi=300)
