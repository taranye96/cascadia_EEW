#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 26 13:27:29 2025

@author: tnye
"""

# Imports 
from os import path, makedirs
from obspy import read, Stream
from glob import glob
import numpy as np
import pandas as pd
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.trigger import plot_trigger
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from datetime import datetime
import matplotlib.dates as mdates
from scipy.stats import gaussian_kde
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator

# Local Imports
import cascadiaEEW_main_fns as emf


home = '/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW'

# Name of batch for simulations
batch = 'cascadia'

# Home directory where waveforms are stored
# home_dir = f'{home}/simulations/{batch}/waveforms_data_curation'
home_dir = '/Users/tnye/cascadia_simulations'

# Read in list of ruptures
ruptures = np.genfromtxt(f'{home}/simulations/{batch}/data/ruptures.list',dtype=str)

# Read in arrival time dataframe
arriv_df = pd.read_csv(f'{home}/event_detection/p_waves/{batch}_parrivals.csv')

# Are we using the polarization filter for the ONC parameter test?
polarize = True

run = 'cascadia.000015'

stns = ['COOS','WYLD']

# Extract event info from rupture .log file
with open(f"{home}/simulations/{batch}/output/ruptures/{run.replace('.','-')}.log") as f:
    lines = f.readlines()
    
    # Get hypocenter coordinates
    hyp_line = lines[16]
    hyp_line = hyp_line.replace(' ', '')
    hyp_line = hyp_line.replace('\n', '')
    hyp_line = hyp_line.replace('(', '')
    hyp_line = hyp_line.replace(')', '')
    hyp_line = hyp_line.split(':')[1].split(',')
    hyplon, hyplat, hypdepth = [float(i) for i in hyp_line]
    
    # Get magnitude
    mag_line = lines[15]
    mag_line = mag_line.replace('\n', '')
    mag = float(mag_line.split(':')[1].split(' ')[2])

fig, axs = plt.subplots(2,2, figsize=(5.75,5))

for istn, stn in enumerate(stns):

    # Gather waveforms
    bb_file_Z = glob(f"{home_dir}/*/*SignalwithNoise/*{run.replace('.','-')}/*{stn}-HNZ.mseed")[0]
    bb_file_N = glob(f"{home_dir}/*/*SignalwithNoise/*{run.replace('.','-')}/*{stn}-HNN.mseed")[0]
    bb_file_E = glob(f"{home_dir}/*/*SignalwithNoise/*{run.replace('.','-')}/*{stn}-HNE.mseed")[0]
    
    noise_file_Z = glob(f"{home_dir}/*/*_Noise/*{run.replace('.','-')}/*HNZ.mseed")[0]
    noise_file_N = glob(f"{home_dir}/*/*_Noise/*{run.replace('.','-')}/*HNN.mseed")[0]
    noise_file_E = glob(f"{home_dir}/*/*_Noise/*{run.replace('.','-')}/*HNE.mseed")[0]
        
    # Read in waveforms
    st_acc_unfilt_E = read(bb_file_E)
    st_acc_unfilt_N = read(bb_file_N)
    st_acc_unfilt_Z = read(bb_file_Z)
    dt = st_acc_unfilt_E[0].stats.delta
    
    # Append values to lists
    rhyp = arriv_df['Rhyp(km)'].iloc[np.where((arriv_df['Station']==stn)&(arriv_df['Batch']==batch)&(arriv_df['Run']==run))[0][0]]
    
    # Get theoretical p-wave arrival
    true_arriv = arriv_df['P-wave arriv'].iloc[np.where((arriv_df['Station']==stn)&(arriv_df['Batch']==batch)&(arriv_df['Run']==run))[0][0]]
    yyyy,mth,dd=true_arriv.split('T')[0].split('-')
    yyyy = int(yyyy)
    mth = int(mth)
    dd = int(dd)
    hh,mm,full_sec = true_arriv.split('T')[1].split(':')
    hh = int(hh)
    mm = int(mm)
    sec = int(full_sec.split('.')[0])
    micsec = int(full_sec.split('.')[1][:-1])
    
    st = Stream()
    st.append(st_acc_unfilt_N[0])
    
    true_p = UTCDateTime(true_arriv)
    sig_ind = np.where(st[0].times('UTCDateTime') >= true_p)[0][0]
    
    ## Test using Obspy with EPIC parameters
    sta, lta, thresh = [0.05,5,20]
    
    try:
        st_acc_filt = emf.highpass(st_acc_unfilt_N,1/15,1/dt,4,zerophase=False)
        st_vel_unfilt = emf.accel_to_veloc(st_acc_filt)
        st_vel_filt = emf.highpass(st_vel_unfilt,1/15,1/dt,4,zerophase=False)
        cft = classic_sta_lta(st_vel_filt[0].data, int(sta * 1/dt), int(lta * 1/dt))
        trig_ind = np.array([], dtype=int)
        obspy_pick_times = np.array([])
        thresh_ind = np.where(np.array(cft) >= thresh)[0]
        thresh_sig = thresh_ind[np.where(thresh_ind >= sig_ind)[0]] 
        pick = st_vel_filt[0].times('UTCDateTime')[thresh_sig[0]]
    except:
        try:
            st_acc_filt = emf.highpass(st_acc_unfilt_E,1/15,1/dt,4,zerophase=False)
            st_vel_unfilt = emf.accel_to_veloc(st_acc_filt)
            st_vel_filt = emf.highpass(st_vel_unfilt,1/15,1/dt,4,zerophase=False)
            cft = classic_sta_lta(st_vel_filt[0].data, int(sta * 1/dt), int(lta * 1/dt))
            trig_ind = np.array([], dtype=int)
            obspy_pick_times = np.array([])
            thresh_ind = np.where(np.array(cft) >= thresh)[0]
            thresh_sig = thresh_ind[np.where(thresh_ind >= sig_ind)[0]] 
            pick = st_vel_filt[0].times('UTCDateTime')[thresh_sig[0]]
        except:
            try:
                st_acc_filt = emf.highpass(st_acc_unfilt_Z,1/15,1/dt,4,zerophase=False)
                st_vel_unfilt = emf.accel_to_veloc(st_acc_filt)
                st_vel_filt = emf.highpass(st_vel_unfilt,1/15,1/dt,4,zerophase=False)
                cft = classic_sta_lta(st_vel_filt[0].data, int(sta * 1/dt), int(lta * 1/dt))
                trig_ind = np.array([], dtype=int)
                obspy_pick_times = np.array([])
                thresh_ind = np.where(np.array(cft) >= thresh)[0]
                thresh_sig = thresh_ind[np.where(thresh_ind >= sig_ind)[0]] 
                pick = st_vel_filt[0].times('UTCDateTime')[thresh_sig[0]]
            except:
                break
        
    ax = axs[0, istn]

    # --- main waveform ---
    ax.plot(st_vel_filt[0].times(),
            st_vel_filt[0].data,
            color='k',
            zorder=3)
    ax.set_ylim(-0.00175,0.0019)
    
    ax.minorticks_on()
    ax.tick_params(axis='x', which='minor', bottom=True, length=2)
    ax.tick_params(axis='x', which='major', bottom=True, length=4)

    ax.xaxis.set_major_locator(MultipleLocator(200))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(0.001))
    ax.yaxis.set_minor_locator(MultipleLocator(0.0005))
    ax.set_title(f'Scenario 15, {stn}')
    
    if istn == 1:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel('Amplitude (cm/s)')
        
    # --- true arrival line on top ---
    ax.vlines(UTCDateTime(true_arriv) - UTCDateTime(st_vel_filt[0].stats.starttime),
              -0.0005,
              0.0005,
              color='blue',
              lw=1,
              label='True arrival',
              zorder=5)
    
    # --- twin axis (CFT) ---
    twin = ax.twinx()
    twin.plot(st_vel_filt[0].times(),
              cft,
              lw=0.5,
              color='darkgray',
              zorder=1)
    twin.set_ylim(0,24)
    twin.yaxis.set_major_locator(MultipleLocator(10))
    twin.yaxis.set_minor_locator(MultipleLocator(5))
    twin.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    
    if len(np.where(cft > thresh)[0]) > 0:
        p_ind = np.where(cft > thresh)[0][0]
        ax.vlines(st_vel_filt[0].times()[p_ind],
                  -0.0005,
                  0.0005,
                  color='red',
                  lw=1,
                  label='STA/LTA pick',
                  zorder=4)
        
    twin.axhline(20, c='red', ls='--', lw=0.5, label='Threshold')
    if istn == 0:
        twin.set_yticklabels([])
    
    # Move twin axis behind
    twin.set_zorder(0)
    ax.set_zorder(1)
    
    # Make foreground axis transparent so twin is visible
    ax.patch.set_visible(False)
    
    # Right-side y-label
    if istn == 1:
        twin.set_ylabel('STA(t) / LTA(t)')
    twin.yaxis.set_label_position('right')
    twin.yaxis.tick_right()
    twin.tick_params(axis='y', which='both', colors='darkgray')
    twin.yaxis.label.set_color('darkgray')
    # twin.spines['right'].set_color('darkgray')
    
    # Remove twin grid
    twin.grid(False)
    
    # --- labels & annotations ---
    ax.set_xlabel('Time (s)')
    ax.text(0.02, 0.98, f'M{round(mag, 1)}',
            ha='left', va='top', transform=ax.transAxes)
    ax.text(0.98, 0.98, r'$R_{hyp}$'+ f' = {round(rhyp)} (km)',
            ha='right', va='top', transform=ax.transAxes)
    
    ax.grid(False)
    
    for spine in ax.spines.values():
        spine.set_color('black')

    for spine in twin.spines.values():
        spine.set_color('black')

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = twin.get_legend_handles_labels()
    
axs[0,0].legend(h1+h2, l1+l2, loc=('upper center'), bbox_to_anchor=(1.05,-0.3), fontsize=9, ncol=3)


# Are we wanting the  'obspy' or 'elarmS' picker for the STA/LTA results using
    # EPIC/ShakeAlert parameters?
shakealert_dtype = 'obspy'

# Read in dataframes
az_df = pd.read_csv(f'{home}/event_detection/azimuth/cascadia_azimuth-amplification.csv')
pick_df = pd.read_csv(f'{home}/event_detection/sta_lta/{batch}_sta_lta_polar.csv')
arriv_df = pd.read_csv(f'{home}/event_detection/p_waves/{batch}_parrivals.csv')

# Merge dataframes
df = pd.merge(arriv_df, pick_df, on=['Batch','Run','Station'], how='outer')
df = pd.merge(df, az_df, on=['Run', 'Station'], how='outer')

# Get azimuth and P-wave amplification factor
az = df['Azimuth']
P = df['P']

# Get STA/LTA triggers (Yes or No) 
trig_OncEW = df['1/5/4.9_triggered']
trig_ShakeAlert = df[f'0.05/5/20_{shakealert_dtype}_triggered']

# Get successful and unsuccessful trigger indices
OncEW_nan_ind = np.where(trig_OncEW=='No')[0]
OncEW_trig_ind = np.where(trig_OncEW=='Yes')[0]
OncEW_only_nan_ind = np.where((trig_OncEW=='No')&(trig_ShakeAlert=='Yes'))[0]
OncEW_only_trig_ind = np.where((trig_OncEW=='Yes')&(trig_ShakeAlert=='No'))[0]
ShakeAlert_nan_ind = np.where(trig_ShakeAlert=='No')[0]
ShakeAlert_trig_ind = np.where(trig_ShakeAlert=='Yes')[0]
ShakeAlert_only_nan_ind = np.where((trig_ShakeAlert=='No')&(trig_OncEW=='Yes'))[0]
ShakeAlert_only_trig_ind = np.where((trig_ShakeAlert=='Yes')&(trig_OncEW=='No'))[0]
both_nan_ind = np.where((trig_ShakeAlert=='No')&(trig_OncEW=='No'))[0]
both_trig_ind = np.where((trig_ShakeAlert=='Yes')&(trig_OncEW=='Yes'))[0]

# Get P-wave pick time residuals
ppick_res_OncEW = df['1/5/4.9_P-arrival Residual (s)']
ppick_res_ShakeAlert = df[f'0.05/5/20_{shakealert_dtype}_P-arrival Residual (s)']

# Get absolute value of pick time residuals
OncEW_res_abs = np.abs(ppick_res_OncEW)
ShakeAlert_res_abs = np.abs(ppick_res_ShakeAlert)

# Get magnitude and distance metrics
# repi = arriv_df['Repi(km)'].values
# rrup = arriv_df['Rrup(km)'].values
rhyp = df['Rhyp(km)'].values
mag = df['Magnitude'].values

ind = np.where(np.isnan(ppick_res_ShakeAlert)==False)[0]
xy = np.vstack([rhyp[ind],ppick_res_ShakeAlert[ind]])
z = gaussian_kde(xy)(xy)
z_norm = (z - z.min()) / (z.max() - z.min())
norm = mcolors.Normalize(vmin=z_norm.min(), vmax=z_norm.max())
cmap = plt.cm.Blues
facecolors = cmap(norm(z_norm))
axs[1,0].minorticks_on()
axs[1,0].tick_params(which='minor', bottom=True, left=True)
axs[1,0].scatter(rhyp[ind],ppick_res_ShakeAlert[ind],lw=0.1,s=15,marker='^',fc=facecolors,ec='k',alpha=1)
axs[1,0].axhline(0,0,ls='--',lw=1,color='k')
axs[1,0].xaxis.set_major_locator(MultipleLocator(500))
axs[1,0].xaxis.set_minor_locator(MultipleLocator(250))
axs[1,0].yaxis.set_major_locator(MultipleLocator(100))
axs[1,0].yaxis.set_minor_locator(MultipleLocator(50))
axs[1,0].grid(alpha=0)
axs[1,0].tick_params(axis='x', which='minor')
axs[1,0].tick_params(axis='y', which='minor')
axs[1,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[1,0].set_xlim(0,1300)
axs[1,0].set_ylim(-425,20)
axs[1,0].set_ylabel(r'$\delta_{arr}$ (s)', fontsize=10)
axs[1,0].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[1,0].set_title('EPIC')

ind = np.where(np.isnan(ppick_res_OncEW)==False)[0]
xy = np.vstack([rhyp[ind],ppick_res_OncEW[ind]])
z = gaussian_kde(xy)(xy)
z_norm = (z - z.min()) / (z.max() - z.min())
norm = mcolors.Normalize(vmin=z_norm.min(), vmax=z_norm.max())
cmap = plt.cm.Blues
facecolors = cmap(norm(z_norm))
axs[1,1].minorticks_on()
axs[1,1].tick_params(which='minor', bottom=True, left=True)
axs[1,1].scatter(rhyp[ind],ppick_res_OncEW[ind],lw=0.1,s=15,marker='^',fc=facecolors,ec='k',alpha=1)
axs[1,1].axhline(0,0,ls='--',lw=1,color='k')
axs[1,1].xaxis.set_major_locator(MultipleLocator(500))
axs[1,1].xaxis.set_minor_locator(MultipleLocator(250))
axs[1,1].yaxis.set_major_locator(MultipleLocator(100))
axs[1,1].yaxis.set_minor_locator(MultipleLocator(50))
axs[1,1].grid(alpha=0)
axs[1,1].tick_params(axis='x', which='minor')
axs[1,1].tick_params(axis='y', which='minor')
axs[1,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[1,1].set_xlim(0,1300)
axs[1,1].set_ylim(-425, 20)
axs[1,1].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[1,1].set_yticklabels([])
axs[1,1].set_title('ONC-EW')

sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cax = fig.add_axes([0.875, 0.08, 0.02, 0.35])  # left, bottom, width, height
cbar = fig.colorbar(sm, cax,
                    shrink=1,   
                    pad=0.15,
                    anchor=(0.0, 1))   
cbar.set_label('Relative point density', labelpad=-10)   
cbar.set_ticks([0, 1])
cbar.set_ticklabels(['Low', 'High'])


plt.subplots_adjust(left=0.15, right=0.85, top=0.925, wspace=0.1, hspace=0.75)

plt.savefig('/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW/manuscript/figures/subfigs/FigS3_STA-LTA.png',dpi=300)