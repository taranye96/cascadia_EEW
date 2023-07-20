#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 14:51:58 2022

@author: tnye
"""

###############################################################################
# Script that adds noise downloaded from the Precise Point Positioning (PPP)
# correction streams from the ONC database to the low frequency synthetic 
# waveforms.
###############################################################################

# Imports
import numpy as np
import pandas as pd
from os import path, makedirs
from glob import glob
from obspy import read
import random
from mpi4py import MPI
from gnss_noise_fns import add_synthetic_gnss_noise

###################### Set up parallelization parameters ######################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

ncpus = size

batch = 'cascadia_longer_wfs'
stn_type = 'onc-onshore-gnss'
run_dirs = sorted(glob(f'/Users/tnye/ONC/simulations/{batch}/output/waveforms/*'))
# run_dirs = ['/Users/tnye/ONC/simulations/leechriver/output/waveforms/leechriver.000009', '/Users/tnye/ONC/simulations/leechriver/output/waveforms/leechriver.000012']

# Read in the noise dfs downloaded from PPP
    # Curently we have the stations AL2H, TFNO, and CMBR
noise_dfs = sorted(glob('/Users/tnye/ONC/data/station_noise/GNSS/PPP/*.csv'))

stn_list = np.genfromtxt(f"/Users/tnye/ONC/data/station_info/{stn_type.split('-')[0]}_list.txt",delimiter=',',dtype=str)[:,0]

############################# Start Parallelization ###########################

# Set up full array of data on main process 
if rank == 0:
    fulldata = np.arange(len(run_dirs), dtype=int)
else:
    fulldata=None

# Number of items on each process
count = len(run_dirs)//ncpus

# Set up empty array for each process to receive data
subdata = np.empty(count, dtype=int)

# Scatter data
comm.Scatter(fulldata,subdata,root=0)

############################### Do Calculations ###############################

for index in subdata:
    
    run = run_dirs[index].split('/')[-1]
    print(f'Rank {rank} working on {run}')

    # Read in the pre-noise low frequency synthetic waveforms 
    raw_wfs = sorted(glob(run_dirs[index]+'/*LY*.sac'))
    
    # Create directories to save mseed files
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/noise/{run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/noise/{run}')
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal_with_noise/{run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal_with_noise/{run}')
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal/{run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal/{run}')
    
    # Group channel files by station
    N = 3
    stn_files = [raw_wfs[n:n+N] for n in range(0, len(raw_wfs), N)]
    
    # Loop through clean waveforms
    for stn_group in stn_files:
        
        # Read in waveform stream and get data
        st_E = read(stn_group[0])
        st_N = read(stn_group[1])
        st_Z = read(stn_group[2])
        
        stn = st_E[0].stats.station
        
        if stn in stn_list:
        
            # # Randomly assign a station's noise
            # w = np.random.randint(0,3)
            # noise_df = pd.read_csv(noise_dfs[w])
            
            # # Get noise for specific channel
            # if channel == 'LYE':
            #     noise_array = np.array(noise_df['east_offset'])
            # elif channel == 'LYN':
            #     noise_array = np.array(noise_df['north_offset'])
            # elif channel == 'LYZ':
            #     noise_array = np.array(noise_df['vertical_offset'])
            
            # if w == 0:
            #     noise_min = 1
            #     noise_max = 227000
            # if w == 1:
            #     noise_min = 1
            #     noise_max = len(noise_array)
            # if w == 2:
            #     noise_min = 100000
            #     noise_max = 180000

            # noise_array = noise_array[noise_min:noise_max]
            
            # # Get length of time series
            # max_t = int(st[0].stats.npts+(120*st[0].stats.sampling_rate))
            
            # # Get random subset of noise to add to timeseries
            # while True:
            #     try:
            #         noise_start = random.choice(range(len(noise_array)-max_t))
            #         noise_end = noise_start + max_t
            #         noise = noise_array[noise_start:noise_end]
            #         True not in np.isnan(noise)
            #     except:
            #         continue
            #     break

            # # Save noise waveform
            # tr_noise = st[0].copy()
            # tr_noise.data = noise - np.mean(noise[:61])
            # # tr_noise.times = st[0].times()[start:end]
            
            # Add 2 minutes of padding to timeseries
            st_E[0].data = np.pad(st_E[0].data, (int(120*st_E[0].stats.sampling_rate),0))
            st_N[0].data = np.pad(st_N[0].data, (int(120*st_N[0].stats.sampling_rate),0))
            st_Z[0].data = np.pad(st_Z[0].data, (int(120*st_Z[0].stats.sampling_rate),0))
            
            # Get synthetic noise
            st_E_noisy,st_N_noisy,st_Z_noisy,E_noise,N_noise,Z_noise = add_synthetic_gnss_noise(st_E,st_N,st_Z,percentile=50)
            
            # Make traces for the noise
            tr_E_noise = st_E[0].copy()
            tr_N_noise = st_N[0].copy()
            tr_Z_noise = st_Z[0].copy()
            
            # Add noise to traces
            tr_E_noise.data = E_noise
            tr_N_noise.data = N_noise
            tr_Z_noise.data = Z_noise
            
            # Change starttime
            st_E[0].stats.starttime='2023-03-22T23:58:00.000000Z'
            st_N[0].stats.starttime='2023-03-22T23:58:00.000000Z'
            st_Z[0].stats.starttime='2023-03-22T23:58:00.000000Z'
            
            # Save padded waveform
            st_E[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal/{run}/{stn}.LYE.mseed', format='MSEED')
            st_N[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal/{run}/{stn}.LYN.mseed', format='MSEED')
            st_Z[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal/{run}/{stn}.LYZ.mseed', format='MSEED')
            
            # Add noise to timeseries
            st_E[0].data = st_E[0].data + E_noise
            st_N[0].data = st_N[0].data + N_noise
            st_Z[0].data = st_Z[0].data + Z_noise
            
            # Save noisy waveform
            st_E[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal_with_noise/{run}/{stn}.LYE.mseed', format='MSEED')
            st_N[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal_with_noise/{run}/{stn}.LYN.mseed', format='MSEED')
            st_Z[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal_with_noise/{run}/{stn}.LYZ.mseed', format='MSEED')
            
            # Save noise waveform
            tr_E_noise.write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/noise/{run}/{stn}.LYE.mseed', format='MSEED')
            tr_N_noise.write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/noise/{run}/{stn}.LYN.mseed', format='MSEED')
            tr_Z_noise.write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/noise/{run}/{stn}.LYZ.mseed', format='MSEED')
            
            
