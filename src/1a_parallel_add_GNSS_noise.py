#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 14:51:58 2022

@author: tnye
"""

###############################################################################
# This script adds synthetic noise to the GNSS synthetic waveforms.
###############################################################################

# Imports
import numpy as np
from os import path, makedirs
from glob import glob
from obspy import read
from mpi4py import MPI
from cascadiaEEW_main_fns import add_synthetic_gnss_noise

###################### Set up parallelization parameters ######################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

ncpus = size

# Name of batch for simulations
batch = 'cascadia_longer_wfs'

# Station type. Only option is:
    # 'onc-onshore-gnss'
stn_type = 'onc-onshore-gnss'

# Gather subdirectory names for the different runs
run_dirs = sorted(glob(f'/Users/tnye/ONC/simulations/{batch}/output/waveforms/*'))

# Get list of stations
stn_list = np.genfromtxt(f"/Users/tnye/ONC/data/station_info/{stn_type.split('-')[0]}_list.txt",delimiter=',',dtype=str)[:,0]

# Percentile of noise to be added. Values used for this study are 10, 50, and 90.
gnss_noise_perc = 90


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
    mod_run = run.replace('.','-')
    print(f'Rank {rank} working on {run}')
    
    # Read in the raw (i.e., no noise added) broadband synthetic waveforms 
    raw_wfs = sorted(glob(run_dirs[index]+'/*LY*.sac'))
    
    # Define names of data and run directories based on station type 
    data_dir = f'Cas-ONC-Onshore-GNSS_{gnss_noise_perc}p'
    run_dir = 'Cas-ONC-On-GNSS'
    
    # Create directories to save mseed files
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Noise/{run_dir}-Noise_{mod_run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Noise/{run_dir}-Noise_{mod_run}')
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Signal/{run_dir}-Sig_{mod_run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Signal/{run_dir}-Sig_{mod_run}')
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_SignalwithNoise/{run_dir}-SigNoise_{mod_run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_SignalwithNoise/{run_dir}-SigNoise_{mod_run}')
    
    # Group channel files by station
    N = 3
    stn_files = [raw_wfs[n:n+N] for n in range(0, len(raw_wfs), N)]
    
    # Loop through clean waveforms
    for stn_group in stn_files:
        
        # Read in waveform stream and get data
        st_E = read(stn_group[0])
        st_N = read(stn_group[1])
        st_Z = read(stn_group[2])
        
        # Get station name
        stn = st_E[0].stats.station
        
        # Loop over stations
        if stn in stn_list:
            
            # Add 2 minutes of padding to timeseries
            st_E[0].data = np.pad(st_E[0].data, (int(120*st_E[0].stats.sampling_rate),0))
            st_N[0].data = np.pad(st_N[0].data, (int(120*st_N[0].stats.sampling_rate),0))
            st_Z[0].data = np.pad(st_Z[0].data, (int(120*st_Z[0].stats.sampling_rate),0))
            
            # Change starttime
            st_E[0].stats.starttime='2023-03-22T23:58:00.000000Z'
            st_N[0].stats.starttime='2023-03-22T23:58:00.000000Z'
            st_Z[0].stats.starttime='2023-03-22T23:58:00.000000Z'
            
            # Get synthetic noise
            st_E_noisy,st_N_noisy,st_Z_noisy,E_noise,N_noise,Z_noise = add_synthetic_gnss_noise(st_E,st_N,st_Z,percentile=gnss_noise_perc)
            
            # Make traces for the noise
            tr_E_noise = st_E[0].copy()
            tr_N_noise = st_N[0].copy()
            tr_Z_noise = st_Z[0].copy()
            
            # Add noise to traces
            tr_E_noise.data = E_noise
            tr_N_noise.data = N_noise
            tr_Z_noise.data = Z_noise
            
            # Add channels to trace stats
            st_E[0].stats.channel = 'LYE'
            st_N[0].stats.channel = 'LYN'
            st_Z[0].stats.channel = 'LYZ'
            st_E_noisy[0].stats.channel = 'LYE'
            st_N_noisy[0].stats.channel = 'LYN'
            st_Z_noisy[0].stats.channel = 'LYZ'
            
            # Save padded waveform
            st_E[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Signal/{run_dir}-Sig_{mod_run}/{run_dir}-Sig_{mod_run}_{stn}-LYE.mseed', format='MSEED')
            st_N[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Signal/{run_dir}-Sig_{mod_run}/{run_dir}-Sig_{mod_run}_{stn}-LYN.mseed', format='MSEED')
            st_Z[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Signal/{run_dir}-Sig_{mod_run}/{run_dir}-Sig_{mod_run}_{stn}-LYZ.mseed', format='MSEED')
            
            # Save noisy waveform
            st_E_noisy[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_SignalwithNoise/{run_dir}-SigNoise_{mod_run}/{run_dir}-SigNoise_{mod_run}_{stn}-LYE.mseed', format='MSEED')
            st_N_noisy[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_SignalwithNoise/{run_dir}-SigNoise_{mod_run}/{run_dir}-SigNoise_{mod_run}_{stn}-LYN.mseed', format='MSEED')
            st_Z_noisy[0].write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_SignalwithNoise/{run_dir}-SigNoise_{mod_run}/{run_dir}-SigNoise_{mod_run}_{stn}-LYZ.mseed', format='MSEED')
            
            # Save noise waveform
            tr_E_noise.write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Noise/{run_dir}-Noise_{mod_run}/{run_dir}-Noise_{mod_run}_{stn}-LYE.mseed', format='MSEED')
            tr_N_noise.write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Noise/{run_dir}-Noise_{mod_run}/{run_dir}-Noise_{mod_run}_{stn}-LYN.mseed', format='MSEED')
            tr_Z_noise.write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Noise/{run_dir}-Noise_{mod_run}/{run_dir}-Noise_{mod_run}_{stn}-LYZ.mseed', format='MSEED')
            
          