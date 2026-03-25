#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 14:27:58 2022

@author: tnye
"""

###############################################################################
# This script runs the simulated waveforms thorugh the Obspy STA/LTA picker.
###############################################################################

# Imports 
from os import path, makedirs
from obspy import read, Stream
from glob import glob
import numpy as np
import pandas as pd
from mpi4py import MPI
from obspy.signal.trigger import classic_sta_lta
from obspy.core import UTCDateTime

# Local Imports
import cascadiaEEW_main_fns as emf

###################### Set up parallelization parameters ######################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

ncpus = size

home = '/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW'

# Name of batch for simulations
batch = 'cascadia'

# Home directory where waveforms are stored
# home_dir = f'{home}/simulations/{batch}/waveforms_data_curation'
home_dir = '/Users/tnye/cascadia_simulations'

# Names of station folders in the waveforms directory
stn_types = ['Cas-ONC-Onshore-StrongMotion','Cas-ONC-Offshore-StrongMotion',
             'Cas-PNSN-Onshore-StrongMotion']

# Read in list of ruptures
ruptures = np.genfromtxt(f'{home}/simulations/{batch}/data/ruptures.list',dtype=str)

# Read in arrival time dataframe
arriv_df = pd.read_csv(f'{home}/event_detection/p_waves/{batch}_arrivals.csv')

pad = 1

# Set up folder to save sta_lta csv files for different ranks (parallelized)
outdir = f'{home}/event_detection/sta_lta/sta_lta_mpi_{pad}s'
if not path.exists(outdir):
    makedirs(outdir, exist_ok=True)

# Set up folder for figures of outlier pick times (large discrepancy between
# STA/LTA pick and the true arrival time) 
figdir = f'{home}/event_detection/plots/outliers/polarized'
if not path.exists(f'{figdir}/positive'):
    makedirs(f'{figdir}/positive', exist_ok=True)
if not path.exists(f'{figdir}/negative'):
    makedirs(f'{figdir}/negative', exist_ok=True)

# Are we using the polarization filter for the ONC parameter test?
polarize = True

############################# Start Parallelization ###########################

if rank == 0:
    sendbuf = np.arange(len(ruptures), dtype=np.int32)

    # count: the size of each sub-task
    ave, res = divmod(sendbuf.size, ncpus)
    count = [ave + 1 if p < res else ave for p in range(ncpus)]
    count = np.array(count)

    # displacement: the starting index of each sub-task
    displ = [sum(count[:p]) for p in range(ncpus)]
    displ = np.array(displ)
else:
    sendbuf = None
    # initialize count on worker processes
    count = np.zeros(ncpus, dtype=int)
    displ = None

# broadcast count
comm.Bcast(count, root=0)

# initialize recvbuf on all processes
recvbuf = np.zeros(count[rank],dtype=np.int32)

comm.Scatterv([sendbuf, count, displ, MPI.INT], recvbuf, root=0)
print('Rank: '+str(rank)+' received data='+str(recvbuf))


############################### Do Calculations ###############################

# Initialize lists
stns = []
batches = []
rupts = []
triggered = []
p_pick = []
sta_list = []
lta_list = []
thresh_list = []
p_res = []

onc_triggered = []
onc_res = []
onc_pick = []
epic_triggered = []
epic_res = []
epic_pick = []
mags = []
rhyp_list = []
    
for index in recvbuf:
    
    krup = int(index)
    rupture = ruptures[krup]
    
    run = rupture.replace('.rupt','').replace('cascadia.','cascadia-')
    
    print(run)
    
    # Loop over station type
    for stn_type in stn_types:

        # Waveform directory
        wf_dir = f'{home_dir}/{stn_type}'
        
        # Gather waveforms
        bb_files_Z = np.array(sorted(glob(f"{wf_dir}/*SignalwithNoise/*{run}/*HNZ.mseed")))
        bb_files_N = np.array(sorted(glob(f"{wf_dir}/*SignalwithNoise/*{run}/*HNN.mseed")))
        bb_files_E = np.array(sorted(glob(f"{wf_dir}/*SignalwithNoise/*{run}/*HNE.mseed")))
        
        noise_files_Z = np.array(sorted(glob(f"{wf_dir}/*_Noise/*{run}/*HNZ.mseed")))
        noise_files_N = np.array(sorted(glob(f"{wf_dir}/*_Noise/*{run}/*HNN.mseed")))
        noise_files_E = np.array(sorted(glob(f"{wf_dir}/*_Noise/*{run}/*HNE.mseed")))
        
        # Loop over waveforms
        for i_file in range(len(bb_files_Z)):
            
            # Read in waveforms
            st_acc_unfilt_E = read(bb_files_E[i_file])
            st_acc_unfilt_N = read(bb_files_N[i_file])
            st_acc_unfilt_Z = read(bb_files_Z[i_file])
            dt = st_acc_unfilt_E[0].stats.delta
            stn =  st_acc_unfilt_E[0].stats.station.split('.')[0]
            
            # Extract event info from rupture .log file
            with open(f"{home}/simulations/{batch}/output/ruptures/{run}.log") as f:
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
            
            # Append values to lists
            rhyp = arriv_df['Rhyp(km)'].iloc[np.where((arriv_df['Station']==stn)&(arriv_df['Batch']==batch)&(arriv_df['Run']==run))[0][0]]
            rhyp_list.append(rhyp)
            stns.append(stn)
            batches.append(batch)
            rupts.append(run)
            mags.append(mag)
            
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
            

            ###################################################################
            
            # Testing usng Obspy with EPIC parameters
            
            ## Vertical component only
            
            st = Stream()
            st.append(st_acc_unfilt_Z[0])
            st[0].stats.channel = 'HNZ'
            
            true_p = UTCDateTime(true_arriv)
            sig_ind = np.where(st[0].times('UTCDateTime') >= true_p-pad)[0][0]
            
            ## Test using Obspy with EPIC parameters
            sta, lta, thresh = [0.05,5,20]
            
            # Highpass filter acceleration at 0.075 Hz
            st_acc_filt_Z = emf.highpass(st_acc_unfilt_Z,1/15,1/dt,4,zerophase=False)
            
            # Integrate to velocity and highpass filter again (EPIC only)
            st_vel_unfilt_Z = emf.accel_to_veloc(st_acc_filt_Z)
            st_vel_filt_Z = emf.highpass(st_vel_unfilt_Z,1/15,1/dt,4,zerophase=False)
        
            # Run obspy STA/LTA algorithm
            cft_Z = classic_sta_lta(st_vel_filt_Z[0].data, int(sta * 1/dt), int(lta * 1/dt))
            
            # Get indices of where STA/LTA passed threshold 
            thresh_ind_Z = np.where(np.array(cft_Z) >= thresh)[0]
            thresh_sig_Z = thresh_ind_Z[np.where(thresh_ind_Z >= sig_ind)[0]]
                                    
            
            # Check for each component if algorithm was triggered
            if len(thresh_sig_Z)>0:
                pick = st_vel_filt_Z[0].times('UTCDateTime')[thresh_sig_Z[0]]
                epic_triggered.append('Yes')
                epic_pick.append(pick)
                epic_res.append(UTCDateTime(true_arriv) - pick)
            else:
                epic_triggered.append('No')
                epic_pick.append(np.nan)
                epic_res.append(np.nan)
            
            
            ###################################################################
            
            ## Test with ONC parameters
            sta, lta, thresh = [1,5,4.9]
            
            # Read in 3 components of noise (used in polarization)
            noise_Z = read(noise_files_Z[i_file])
            noise_N = read(noise_files_N[i_file])
            noise_E = read(noise_files_E[i_file])
            
            # Highpass filter acceleration at 0.075 Hz
            st_Z_filt = emf.highpass(st_acc_unfilt_Z,1/15,1/dt,4,zerophase=False)
            st_N_filt = emf.highpass(st_acc_unfilt_N,1/15,1/dt,4,zerophase=False)
            st_E_filt = emf.highpass(st_acc_unfilt_E,1/15,1/dt,4,zerophase=False)
            
            # Lowpass filter acceleration at 10 Hz for ONC parameters
            st_Z_filtfilt = emf.lowpass(st_Z_filt,10,1/dt,4,zerophase=False)
            st_N_filtfilt = emf.lowpass(st_N_filt,10,1/dt,4,zerophase=False)
            st_E_filtfilt = emf.lowpass(st_E_filt,10,1/dt,4,zerophase=False)
            
            if polarize == True:
                st_P, st_S = emf.polarization(st_Z_filtfilt,st_N_filtfilt,st_E_filtfilt,
                                   noise_Z,noise_N,noise_E)
                tr_Z, tr_N, tr_E = st_P
                data = tr_Z.data**2 + tr_N.data**2 +tr_E.data**2
            else:
                data = st_Z_filtfilt[0]
            
            
            # Pick P-wave using Obspy
            cft = classic_sta_lta(data, int(sta * 1/dt), int(lta * 1/dt))
            
            thresh_ind = np.where(np.array(cft) >= thresh)[0]
            thresh_sig = thresh_ind[np.where(thresh_ind >= sig_ind)[0]]
            
            
            
            # Determine if algorithm was triggered
            if len(thresh_sig)>0:
                pick = st_Z_filt[0].times('UTCDateTime')[thresh_sig[0]]
                onc_triggered.append('Yes')
                onc_pick.append(pick)
                onc_res.append(UTCDateTime(true_arriv) - pick)
                
            else:
                onc_triggered.append('No')
                onc_pick.append(np.nan)
                onc_res.append(np.nan)
            
    
# Make dataframe
data = {'Batch':batches,'Run':rupts,'Magnitude':mags,'Station':stns,
        '0.05/5/20_triggered':epic_triggered,
        '0.05/5/20_P-picktime':epic_pick,'0.05/5/20_P-arrival Residual (s)':epic_res,
        '1/5/4.9_triggered':onc_triggered, '1/5/4.9_P-picktime':onc_pick,
        '1/5/4.9_P-arrival Residual (s)':onc_res}

df = pd.DataFrame(data)
df.to_csv(f'{outdir}/{batch}_sta_lta_{rank}.csv', index=False)


# Ensure all file I/O is finished before any rank exits
comm.Barrier()

# Explicit MPI shutdown (defensive; avoids "exited improperly")
MPI.Finalize()


    
