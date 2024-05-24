#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 13:25:48 2023

@author: tnye
"""

##############################################################################
# This script calculates the SNR of the simulated GNSS waveforms that have been
# combined with noise. The results are used in Figure 8.
##############################################################################

# Imports
import numpy as np
import pandas as pd
from glob import glob
from obspy import read, UTCDateTime

# Noise percentile (10, 50, 90)
noise_perc = 90

# Read in rupture list
ruptures = np.genfromtxt(f'/Users/tnye/ONC/simulations/cascadia/data/ruptures.list',dtype=str)

# Read in df with P- and S-wave arrival times
arriv_df = pd.read_csv('/Users/tnye/ONC/magnitude_estimation/cascadia_arrivals_taupy.csv')

# Initialize lists
event_list = []
stn_list = []
mag_list = []
snr_list = []

# Loop over ruptures 
for rupture in ruptures:
    
    # Get run number
    run = rupture.replace('.rupt','')
    print(run)
    
    # Read in rupture log file
    mod_run = run.replace('.','-')
    logfile = f"/Users/tnye/ONC/simulations/cascadia/output/ruptures/{rupture.replace('.rupt','.log')}"
    f = open(logfile, 'r')
    lines = f.readlines()
    
    # Get event magnitude
    mag = float(lines[15].split(' ')[-1].split('\n')[0])
    
    # Gather mseed files
    mseed_files = sorted(glob(f'/Users/tnye/ONC/simulations/cascadia/waveforms_data_curation/Cas-ONC-Onshore-GNSS_{noise_perc}p/Cas-ONC-On-GNSS_SignalwithNoise/*{mod_run}/*.mseed'))
    
    # Group channel files by station
    N = 3
    stn_files = [mseed_files[n:n+N] for n in range(0, len(mseed_files), N)]
    
    # Loop over stations
    for group in stn_files:
        
        # Read in components
        st_E = read(group[0])
        st_N = read(group[1])
        st_Z = read(group[2])
        stn = st_E[0].stats.station
        
        # Append values to lists
        mag_list.append(mag)
        event_list.append(run)
        stn_list.append(stn)
        
        # Find arrival time df row for event-station pair
        arriv_ind = np.where((arriv_df['Run']==run) & (arriv_df['Station']==stn))[0][0]
        
        # Get S-wave arrival times
        s_arr = arriv_df['S-wave arriv'].values[arriv_ind]
        s_pick = (np.abs(np.array(st_E[0].times('UTCDateTime'))-UTCDateTime(s_arr))).argmin()
        
        # Get signal and noise segmetns
        sig_data_E = st_E[0].data[s_pick-5:s_pick+60]
        noise_data_E = st_E[0].data[:120]
        sig_data_N = st_N[0].data[s_pick-5:s_pick+60]
        noise_data_N = st_N[0].data[:120]
        sig_data_Z = st_Z[0].data[s_pick-5:s_pick+60]
        noise_data_Z = st_Z[0].data[:120]
        
        # Compute SNR
        snr_E = np.std(sig_data_E)/np.std(noise_data_E) 
        snr_N = np.std(sig_data_N)/np.std(noise_data_N) 
        snr_Z = np.std(sig_data_Z)/np.std(noise_data_Z) 
        
        snr_all = np.sqrt(snr_E**2 + snr_N**2 + snr_Z**2)
        
        # Append SNR value to list
        snr_list.append(snr_all)

# Save SNR values to df
data = {'Event':event_list, 'Station':stn_list, 'Mag':mag_list, 'SNR':snr_list}
df = pd.DataFrame(data)
df.to_csv(f'/Users/tnye/ONC/groundmotion_analysis/snr_files/snr_{noise_perc}p.csv')
        
        