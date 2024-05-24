#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 11:39:21 2022

@author: tnye
"""

###############################################################################
# This parallelized script adds real noise to the broadband synthetic
# waveforms.
###############################################################################

# Imports
import numpy as np
from os import path, makedirs
from glob import glob
from obspy import read
import random
from mpi4py import MPI

###################### Set up parallelization parameters ######################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

ncpus = size

# Name of batch for simulations
batch = 'cascadia'

# Station type. Options are:
    # 'onc-onshore-strongmotion'
    # 'onc-offshore-strongmotion'
    # 'pnsn-strongmotion'
stn_type = 'pnsn-strongmotion'

# Home directory where waveforms are stored
home_dir = f'/Users/tnye/ONC/simulations/{batch}/output/waveforms'

# Gather subdirectory names for the different runs
run_dirs = sorted(glob(f'/Users/tnye/ONC/simulations/{batch}/output/waveforms/*'))

# Read in the noise wfs processed from the raw station data
        # Curently we have the stations AL2H, CMBR, and PHRB
onshore_noise_wfs = glob('/Users/tnye/ONC/data/station_noise/broadband/corrected_noise/*.mseed')

# Get list of stations
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
    mod_run = run.replace('.','-')
    print(f'Rank {rank} working on {run}')

    # Read in the raw (i.e., no noise added) broadband synthetic waveforms 
    raw_wfs = sorted(glob(f'{home_dir}/{run}/*bb.HN*.sac'))
    
    # Define names of data and run directories based on station type 
    if stn_type == 'onc-onshore-strongmotion':
        data_dir = 'Cas-ONC-Onshore-StrongMotion'
    elif stn_type == 'onc-offshore-strongmotion':
        data_dir = 'Cas-ONC-Offshore-StrongMotion'
    else:
        data_dir = 'Cas-PNSN-Onshore-StrongMotion'
    
    if stn_type == 'onc-onshore-strongmotion':
        run_dir = 'Cas-ONC-On-SM'
    elif stn_type == 'onc-offshore-strongmotion':
        run_dir = 'Cas-ONC-Off'
    else:
        run_dir = 'Cas-PNSN'
    
    # Create directories to save mseed files
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Noise/{run_dir}-Noise_{mod_run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Noise/{run_dir}-Noise_{mod_run}')
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Signal/{run_dir}-Sig_{mod_run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Signal/{run_dir}-Sig_{mod_run}')
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_SignalwithNoise/{run_dir}-SigNoise_{mod_run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_SignalwithNoise/{run_dir}-SigNoise_{mod_run}')
    
    # Loop through raw waveforms
    for wf in raw_wfs:
        
        # Read in waveform stream and get data
        st = read(wf)
        stn = st[0].stats.station
        channel = wf.split('/')[-1].split('.')[-2]
        
        # Check that station is one of your desired station types
        if stn.split('.')[0] in stn_list:
            
            # If station is offshore, use corresponding offshore noise
            if 'offshore' in stn_type:
                # Read in noise stream
                noise_st = read(glob(f'/Users/tnye/ONC/data/station_noise/offshore/trimmed/*{stn}*2018-10-31.mseed')[0])
                
                # Randomly assign a segment of day (even for 1st part, odd for 2nd)
                w = np.random.randint(0,2)
                if w == 0:
                    # Get noise for specific channel
                    for i, noise_tr in enumerate(noise_st):
                        if i % 2 == 0:
                            if noise_tr.stats.channel == channel:
                                break
                elif w == 1:
                    # Get noise for specific channel
                    for i, noise_tr in enumerate(noise_st):
                        if i % 2 != 0:
                            if noise_tr.stats.channel == channel:
                                break
            
            # If station is not offshore, assign a random noise segment
            else:
                # Randomly assign a station's noise
                w = np.random.randint(0,3)
                noise_st = read(onshore_noise_wfs[w])
            
                # Get noise for specific channel
                for noise_tr in noise_st:
                    if noise_tr.stats.channel == channel:
                        break
            
            # Get ratio of sampling rate between noise and signal
            samp_ratio = int(noise_tr.stats.sampling_rate/st[0].stats.sampling_rate)
            
            # Get length of time series with 2 minutes padding added
            max_t = int(st[0].stats.npts+(120*st[0].stats.sampling_rate))
            
            # Multiply max length by sampling ratio to account for difference
                # in sampling rates
            length_noise = max_t * samp_ratio
            
            # Get random subset of noise to add to timeseries
            start = random.choice(range(len(noise_tr.data)-length_noise)) # start has to be early enough for noise to be correct length for timeseries
            end = start + length_noise
            noise = noise_tr.data[start:end][::samp_ratio] # ::samp_ratio accounts for difference in sampling rates between noise and data
            
            # De-mean noise
            noise -= np.mean(noise)
            
            # Add 2 minutes of padding to timeseries
            st[0].data = np.pad(st[0].data, (int(120*st[0].stats.sampling_rate),0))
            
            # Make sure starttimes are consistent
            st[0].stats.starttime='2023-03-22T23:58:00.000000Z'
            
            # Remove periods from file names for uploading to Borealis
            if '.' in stn:
                stn = stn.replace('.','-')
            
            # Save padded waveform
            filename = f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Signal/{run_dir}-Sig_{mod_run}/{run_dir}-Sig_{mod_run}_{stn}-{channel}.mseed'
            st[0].stats.channel = channel
            st[0].write(filename, format='MSEED')
            
            # Save noise
            st_noise  = st.copy()
            st_noise[0].data = noise
            st_noise[0].times = st[0].times()[start:end][::samp_ratio]
            st_noise[0].stats.channel = channel
            filename = f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_Noise/{run_dir}-Noise_{mod_run}/{run_dir}-Noise_{mod_run}_{stn}-{channel}.mseed'
            st_noise.write(filename, format='MSEED')
            
            # Save noisy waveform
            st_noisy = st.copy()
            st_noisy[0].data = st[0].data + noise
            st_noisy[0].stats.channel = channel
            filename = f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{data_dir}/{run_dir}_SignalwithNoise/{run_dir}-SigNoise_{mod_run}/{run_dir}-SigNoise_{mod_run}_{stn}-{channel}.mseed'
            st_noisy[0].write(filename, format='MSEED')
            
