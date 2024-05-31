#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 14:24:57 2022

@author: tnye
"""

###############################################################################
# Script that calculates the amplitude of P-waves (Pd) using a specified
# window length.
###############################################################################

# Imports
from os import path, makedirs
import numpy as np
import pandas as pd
from glob import glob
from mpi4py import MPI
from obspy import read
from scipy.integrate import cumtrapz
from obspy.core import UTCDateTime

# Local imports
import cascadiaEEW_main_fns as emf

###################### Set up parallelization parameters ######################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

ncpus = size

# Name of batch for simulations
batch = 'cascadia' 

# Station name descriptors (used in names of waveform folders)
stn_types = ['ONC-Onshore','ONC-Offshore','PNSN']

# Read in list of ruptures
ruptures = np.genfromtxt(f'/Users/tnye/ONC/simulations/{batch}/data/ruptures.list',dtype=str)

# Time windoes used to calculate Pd
windows = [1,2,4,6,12,20]

# Read in arrival time df
arriv_df = pd.read_csv('/Users/tnye/ONC/event_detection/p_waves/cascadia_parrivals.csv')

############################# Start Parallelization ###########################

if rank == 0:
    sendbuf = np.arange(float(len(ruptures)))

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
recvbuf = np.zeros(count[rank])

comm.Scatterv([sendbuf, count, displ, MPI.DOUBLE], recvbuf, root=0)
print('Rank: '+str(rank)+' received data='+str(recvbuf))


############################### Do Calculations ###############################

# Make directory to save Pd csv files (parallelize)
if not path.exists(f'/Users/tnye/ONC/magnitude_estimation/Pd_files/Pd_mpi_files/{batch}/cascadia_noisy'):
    makedirs(f'/Users/tnye/ONC/magnitude_estimation/Pd_files/Pd_mpi_files/{batch}/cascadia_noisy')

# Initialize a file to save Pd values to
outfile = open(f'/Users/tnye/ONC/magnitude_estimation/Pd_files/Pd_mpi_files/{batch}/cascadia_noisy/Pd_{rank}.csv', 'w')
outfile.write('Event,Station,Station Type,Magnitude,Rhyp,Repi,1s_Pd_logcm,2s_Pd_logcm,4s_Pd_logcm,6s_Pd_logcm,12s_Pd_logcm,20s_Pd_logcm\n')

for index in recvbuf:

    rupture = ruptures[int(index)]
    run = rupture.strip('.rupt')
    modified_run = run.replace('.','-')
    
    print(run)
    
    # Loop over station types
    for stations in stn_types:
    
        # Gather mseed files
        files = sorted(glob(f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation/*{stations}*/*SignalwithNoise/*{modified_run}/*HNZ*'))
        
        # Get earthquake info for this batch and run
        f = open(f'/home/tnye/onc/simulations/{batch}/output/ruptures/{run}.log')
        lines = f.readlines()
        mag_line = lines[15]
        mag = float(mag_line.split(' ')[3].split('\n')[0])
        f.close()
        
        # Loop over mseed files
        for file in files:
            
            ## Read in waveform and get station name and sampling rate
            acc_st = read(file)
            stn = acc_st[0].stats.station.strip('.')
            stsamprate = acc_st[0].stats.sampling_rate
            
            print(f'{stn}, {run}')

            # Get P- and S-arrivals
            idx = np.where((np.array(arriv_df['Station'])==stn) &
                            (np.array(arriv_df['Batch'])==batch) & 
                            (np.array(arriv_df['Run'])==run))[0][0]
       
            p_arriv = UTCDateTime(arriv_df['P-wave arriv'].iloc[idx])
            s_arriv = UTCDateTime(arriv_df['S-wave arriv'].iloc[idx])
            
            ## Demean data
            idx_p = (np.abs(np.array(acc_st[0].times('UTCDateTime'))-p_arriv)).argmin()
            if int(idx_p-(10*stsamprate)) >=0:
                mean = np.mean(acc_st[0].data[int(idx_p-(10*stsamprate)):idx_p])
            else:
                mean = np.mean(acc_st[0].data[:idx_p+1])
            acc_st[0].data = acc_st[0].data - mean
            
            ## Highpass filter to 0.075 Hz
            acc_filt = emf.highpass(acc_st,0.075,stsamprate,4,zerophase=False)
            
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
            disp_filt = emf.bandpass(disp_st,0.075,3.0,stsamprate,4,zerophase=False)
            
            ## Distance correction
            rhyp = arriv_df['Rhyp(km)'].iloc[idx]
            repi = arriv_df['Repi(km)'].iloc[idx]
            
            ## Get start time of time window
            idx_start = (np.abs(np.array(disp_st[0].times('UTCDateTime'))-p_arriv)).argmin()
            
            # Initialize Pd list
            Pd_list = []
            
            # Loop over time windows
            for window in windows:
                if window <= 0.95*(s_arriv-p_arriv):
                    
                    ## Get end time of time window
                    idx_end = int(idx_start + window*stsamprate)
                    
                    # Calculate Pd
                    Pd = np.max(np.abs(disp_filt[0].data[idx_start:idx_end+1]))
                
                    ## Get log and correct to 10 km
                    corrected_Pd = np.log10(Pd*100)+np.log10(rhyp/10)
                    Pd_list.append(corrected_Pd)
                else:
                    Pd_list.append(np.nan)
                
            # Add Pd values to file
            line = f'{run},{stn},{stations},{mag},{rhyp},{repi},{Pd_list[0]},{Pd_list[1]},{Pd_list[2]},{Pd_list[3]},{Pd_list[4]},{Pd_list[5]}\n'
            outfile.write(line)
                
# Close Pd file
outfile.close()
        
