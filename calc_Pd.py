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
import numpy as np
import pandas as pd
from glob import glob
from obspy import read
from scipy.integrate import cumtrapz
from obspy.core import UTCDateTime
import tsueqs_main_fns as tmf

# batch = 'small_events'
batch = 'cascadia_longer_wfs' 
stations = 'onc_onshore-strongmotion'

# files = glob('/Users/tnye/ONC/simulations/batch_09/waveforms/eews.000000/*.bb.HNZ*')
# files = glob('/Users/tnye/ONC/simulations/*/waveforms/*/*.bb.HNZ*')
ruptures = sorted(glob(f'/Users/tnye/ONC/simulations/{batch}/SignalwithNoise/{stations}/*'))

windows = [1,2,4,6,12,20]

# arriv_df = pd.read_csv(f'/Users/tnye/ONC/outfiles/p_waves/{batch}_parrivals.csv')
arriv_df = pd.read_csv(f'/Users/tnye/ONC/outfiles/p_waves/{batch}_parrivals_{stations}_vel.csv')

# outfile = open(f'/Users/tnye/ONC/files/Pd_4s_{batch}.txt', 'w')
outfile = open(f'/Users/tnye/ONC/files/Pd_{batch}_{stations}_withnoise.txt', 'w')
outfile.write('#Magnitude \t 1s_Pd_logcm \t 2s_Pd_logcm \t 4s_Pd_logcm \t 6s_Pd_logcm \t 12s_Pd_logcm \t 20s_Pd_logcm\n')

for folder in ruptures:
    
    # Get batch and run number
    batch = folder.split('/')[5]
    run = folder.split('/')[-1]
    
    print(run)
    
    # Get mseed files for this batch and run
    # files = glob(f'/Users/tnye/ONC/simulations/{batch}/waveforms/{run}/*.bb.HNZ*')
    files = glob(f'/Users/tnye/ONC/simulations/{batch}/Signal/{stations}/{run}/*.HNZ*')
    
    # Get earthquake info for this batch and run
    f = open(f'/Users/tnye/ONC/simulations/{batch}/ruptures/{run}.log')
    lines = f.readlines()
    mag_line = lines[15]
    mag = float(mag_line.split(' ')[3].split('\n')[0])
   
    # x = open(f'/Users/tnye/ONC/files/Pd_4s_{batch}.txt', 'w')
    # x.write('#Magnitude,Pd_logcm\n')
    
    # Loop through mseed files
    for file in files:
        
        acc_st = read(file)
        
        stn = acc_st[0].stats.station
        stsamprate = acc_st[0].stats.sampling_rate
        
        ## Demean data
        idx = np.where((np.array(arriv_df['Station'])==stn) &
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
        acc_filt = tmf.highpass(acc_st,0.075,stsamprate,4,zerophase=True)
        
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
        disp_filt = tmf.bandpass(disp_st,0.075,3.0,stsamprate,4,zerophase=True)
        
        ## Distance correction
        rhyp = arriv_df['Rhyp(km)'].iloc[np.where((np.array(arriv_df['Station']==stn))&
                                          (np.array(arriv_df['Batch']==batch))&
                                            (np.array(arriv_df['Run']==run)))[0][0]]
        
        ## Get start and end time of 4 second p-wave interval
        idx_start = (np.abs(np.array(disp_st[0].times('UTCDateTime'))-p_arriv)).argmin()
        
        file_data = [mag]
        for window in windows:
            idx_end = int(idx_start + window*stsamprate)
        
            Pd = np.max(np.abs(disp_filt[0].data[idx_start:idx_end+1]))
            
            ## Get log and correct to 10 km
            file_data.append(np.log10(Pd*100)+np.log10(rhyp/10))
        
        file_data = np.array(file_data,dtype=object)
        np.savetxt(outfile, [file_data], fmt='%E', delimiter='\t')
    
    f.close()
        
outfile.close()
        