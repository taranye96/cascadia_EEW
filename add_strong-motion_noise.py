#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 11:39:21 2022

@author: tnye
"""

###############################################################################
# Script that adds real noise obtained from Angela to the broadband synthetic 
# waveforms.
###############################################################################

# Imports
import numpy as np
import pandas as pd
from os import path, makedirs
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

batch = 'cascadia_longer_wfs'
stn_type = 'onc-offshore-strongmotion'
# home_dir = f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/raw_waveforms'
home_dir = f'/Users/tnye/ONC/simulations/{batch}/output/waveforms'

run_dirs = sorted(glob(f'/Users/tnye/ONC/simulations/{batch}/output/waveforms/*'))
# run_dirs = ['/Users/tnye/ONC/simulations/leechriver/output/waveforms/leechriver.000009', '/Users/tnye/ONC/simulations/leechriver/output/waveforms/leechriver.000012']

# Read in the noise wfs processed from the raw station data Angela provided
        # Curently we have the stations AL2H, CMBR, and PHRB
noise_wfs = glob('/Users/tnye/ONC/data/station_noise/broadband/corrected_noise/*.mseed')

stn_list = np.genfromtxt('/Users/tnye/ONC/data/station_info/onc-offshore_list.txt',delimiter=',',dtype=str)[:,0]

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

# print('Rank: '+str(rank)+' received data='+str(subdata))

############################### Do Calculations ###############################

for index in subdata:
    
    run = run_dirs[index].split('/')[-1]
    print(f'Rank {rank} working on {run}')

    # Read in the pre-noise broadband synthetic waveforms 
    raw_wfs = sorted(glob(f'{home_dir}/{run}/*bb.HN*.sac'))
    
    # Create directories to save mseed files
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/noise/{run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/noise/{run}')
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal_with_noise/{run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal_with_noise/{run}')
    if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal/{run}'):
        makedirs(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal/{run}')
    
    # Loop through clean waveforms
    for wf in raw_wfs:
        
        # Read in waveform stream and get data
        st = read(wf)
        stn = st[0].stats.station
        channel = wf.split('/')[-1].split('.')[-2]
        
        # # Search for stn files
        # files = glob(f'/Users/tnye/ONC/simulations/{batch}/SignalwithNoise/{stn_type}/{run}/{stn}*')
        
        # if len(files) == 0:
        
        if stn.split('.')[0] in stn_list:
        
            # Randomly assign a station's noise
            w = np.random.randint(0,3)
            noise_st = read(noise_wfs[w])
            
            # Get noise for specific channel
            for noise_tr in noise_st:
                if noise_tr.stats.channel == channel:
                    break
            
            samp_ratio = int(noise_tr.stats.sampling_rate/st[0].stats.sampling_rate)
            
            # Get length of time series
            max_t = int(st[0].stats.npts+(120*st[0].stats.sampling_rate))
            
            # Multiply max length by sampling ratio to account for difference in sampling rates
            length_noise = max_t * samp_ratio
            
            # Get random subset of noise to add to timeseries
            start = random.choice(range(len(noise_tr.data)-length_noise)) # start has to be early enough for noise to be correct length for timeseries
            end = start + length_noise
            noise = noise_tr.data[start:end][::samp_ratio] # ::samp_ratio accounts for difference in sampling rates between noise and data
            
            # Save noise waveform
            tr_noise = noise_tr.copy()
            tr_noise.data = noise
            tr_noise.times = st[0].times()[start:end][::samp_ratio]
            
            # Add 2 minutes of padding to timeseries
            st[0].data = np.pad(st[0].data, (int(120*st[0].stats.sampling_rate),0))
            
            # Change starttime
            st[0].stats.starttime='2023-03-22T23:58:00.000000Z'
            
            # Save padded waveform
            filename = f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal/{run}/{stn}.{channel}.mseed'
            st[0].write(filename, format='MSEED')
            
            # Add noise to timeseries
            st[0].data = st[0].data + noise
            
            # # Make sure output folder exists for noise
            # if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/noise_wfs/strong-motion/{run}/'):
            #     makedirs(f'/Users/tnye/ONC/simulations/{batch}/noise_wfs/strong-motion/{run}/')
            
            tr_noise.write(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/noise/{run}/{stn}.{channel}.mseed', format='MSEED')
           
            # # Make sure output folder exists for noisy waveforms
            # if not path.exists(f'/Users/tnye/ONC/simulations/{batch}/wfs_with_noise/strong-motion/pnsn/{run}'):
            #     makedirs(f'/Users/tnye/ONC/simulations/{batch}/wfs_with_noise/strong-motion/{run}')
            
            # Save noisy waveform
            filename = f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/{stn_type}/signal_with_noise/{run}/{stn}.{channel}.mseed'
            st[0].write(filename, format='MSEED')
            
