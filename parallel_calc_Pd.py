#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 14:13:07 2022

@author: tnye
"""

###############################################################################
# Parallelized script that calculates the amplitude of P-waves (Pd) using a
# specified window length.
###############################################################################

# Imports
from mpi4py import MPI
import numpy as np
import pandas as pd
from os import path, makedirs
from glob import glob
from obspy import read
from scipy.integrate import cumtrapz
from obspy.core import UTCDateTime
import filtering_fns as filt

###################### Set up parallelization parameters ######################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

ncpus = size

batch = 'cascadia'
wf_folder = 'waveforms'
ruptures = sorted(glob(f'/Users/tnye/ONC/simulations/{batch}/{wf_folder}/*'))

windows = [1,2,4,6,12,20]
arriv_df = pd.read_csv('/Users/tnye/ONC/outfiles/p_waves/'+batch.split('/')[-1]+'_parrivals.csv')

if not path.exists(f'/Users/tnye/ONC/outfiles/Pd/mpi_files/{batch}'):
    makedirs(f'/Users/tnye/ONC/outfiles/Pd/mpi_files/{batch}')
if not path.exists(f'/Users/tnye/ONC/outfiles/Pd/collected_files/'):
    makedirs(f'/Users/tnye/ONC/outfiles/Pd/collected_files/')


############################# Start Parallelization ###########################

# Set up full array of data on main process 
if rank == 0:
    fulldata = np.arange(len(ruptures), dtype=int)
else:
    fulldata=None

# Number of items on each process
ecount = len(ruptures)//ncpus

# Set up empty array for each process to receive data
subdata = np.empty(count, dtype=int)

# Scatter data
comm.Scatter(fulldata,subdata,root=0)

print('Rank: '+str(rank)+' received data='+str(subdata))


############################# Start Computations ##############################

mag_list = []
Pd_1 = []
Pd_2 = []
Pd_4 = []
Pd_6 = []
Pd_12 = []
Pd_20 = []

for index in subdata:
    
    rupture = ruptures[index]

    flatfile_path = open(f'/Users/tnye/ONC/outfiles/Pd/mpi_files/{batch}/Pd_mpi_{rank}.csv', 'w')
        
    # Get batch and run number
    run = rupture.split('/')[-1].replace('.rupt','')
    
    if path.exists(f'/Users/tnye/ONC/flatfiles/IMs/{batch}/{run}.csv'):
    
        # Get mseed files for this batch and run
        files = glob(f'/Users/tnye/ONC/simulations/{batch}/{wf_folder}/{run}/*.bb.HNZ*')
        
        # Get earthquake info for this batch and run
        f = open(f'/Users/tnye/ONC/simulations/{batch}/ruptures/{run}.log')
        lines = f.readlines()
        mag_line = lines[15]
        mag = float(mag_line.split(' ')[3].split('\n')[0])
        
        # Loop through mseed files
        for file in files:
            
            acc_st = read(file)
            
            stn = acc_st[0].stats.station.split('.')[0]
            stsamprate = acc_st[0].stats.sampling_rate
            
            ## Demean data
            idx = np.where((np.array(arriv_df['Station']==stn)) &
                            (np.array(arriv_df['Batch'])==batch) & 
                            (np.array(arriv_df['Run'])==run))[0][0]
            p_arriv = UTCDateTime(arriv_df['P-wave arriv'].iloc[idx])
            idx_p = (np.abs(np.array(acc_st[0].times('UTCDateTime'))-p_arriv)).argmin()
            if int(idx_p-(10*stsamprate)) >=0:
                mean = np.mean(acc_st[0].data[int(idx_p-(10*stsamprate)):idx_p])
            else:
                mean = np.mean(acc_st[0].data[:idx_p+1])
            acc_st[0].data = acc_st[0].data - mean
            
            ## Highpass filter to 0.075 Hz
            acc_filt = filt.highpass(acc_st,0.075,stsamprate,4,zerophase=True)
            
            ## Get demeaned amplitudes:
            acc_amplitude = acc_filt[0].data
            
            ## And times:
            acc_times = acc_filt[0].times()
            
            ## Double integrate the acceration time series to get displacement:
            vel_amplitude = cumtrapz(acc_amplitude,x=acc_times,initial=0)
            disp_amplitude = cumtrapz(vel_amplitude,x=acc_times,initial=0)
            
            ## Make a copy of the old stream object:
            disp_st = acc_st.copy()
            
            ## Put the integrated displacement in there in place of accel:
            disp_st[0].data = disp_amplitude
        
            ## Bandpass filter between 0.075 and 3.0 Hz
            disp_filt = filt.bandpass(disp_st,0.075,3.0,stsamprate,4,zerophase=True)
            
            ## Distance correction
            rhyp = arriv_df['Rhyp(km)'].iloc[np.where((np.array(arriv_df['Station']==stn))&
                                              (np.array(arriv_df['Batch']==batch))&
                                                (np.array(arriv_df['Run']==run)))[0][0]]
            
            ## Get start and end time of 4 second p-wave interval
            idx_start = (np.abs(np.array(disp_st[0].times('UTCDateTime'))-p_arriv)).argmin()
            
            Pds = []
            for window in windows:
                idx_end = int(idx_start + window*stsamprate)
            
                Pd = np.max(np.abs(disp_filt[0].data[idx_start:idx_end+1]))
                
                ## Get log and correct to 10 km
                Pds.append(np.log10(Pd*100)+np.log10(rhyp/10))
            
            mag_list.append(mag)
            Pd_1.append(Pds[0])
            Pd_2.append(Pds[1])
            Pd_4.append(Pds[2])
            Pd_6.append(Pds[3])
            Pd_12.append(Pds[4])
            Pd_20.append(Pds[5])
        
        f.close()
        
dataset_dict = {'Magnitude':mag_list,'1s_Pd_logcm':Pd_1,'2s_Pd_logcm':Pd_2,
                '4s_Pd_logcm':Pd_4,'6s_Pd_logcm':Pd_6,'12s_Pd_logcm':Pd_12,
                '20s_Pd_logcm':Pd_20}

df = pd.DataFrame(data=dataset_dict)
df.to_csv(flatfile_path,index=False)


############################# End Parallelization #############################

#Set up empty array to gather data on
recvbuf=None
if rank == 0:
    recvbuf = np.empty(count*size, dtype=int)

# comm.Gather(subdata, recvbuf, root=0)


############################## Combine mpi files ##############################

# Gather mpi files
mpi_files = glob(f'/Users/tnye/ONC/outfiles/Pd/mpi_files/{batch}/*.csv')

# Combine mpi files into one dataframe
df = pd.concat(map(pd.read_csv, mpi_files), ignore_index=True)    

# Sort dataframe by magnitude
df = df.sort_values(by=['Magnitude'])

# Save dataframe   
df.to_csv('/Users/tnye/ONC/outfiles/Pd/collected_files/'+batch.split('/')[-1]+'.csv')





