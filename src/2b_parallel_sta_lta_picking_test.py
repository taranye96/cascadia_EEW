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
from obspy.signal.trigger import plot_trigger
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from datetime import datetime
import matplotlib.dates as mdates

# Local Imports
import cascadiaEEW_main_fns as emf
import elarmS_picker as picker

###################### Set up parallelization parameters ######################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

ncpus = size

# Name of batch for simulations
batch = 'cascadia'

# Home directory where waveforms are stored
home_dir = f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation'

# Names of station folders in the waveforms directory
stn_types = ['Cas-ONC-Onshore-StrongMotion','Cas-ONC-Offshore-StrongMotion',
             'Cas-PNSN-Onshore-StrongMotion']

# Read in list of ruptures
ruptures = np.genfromtxt(f'/Users/tnye/ONC/simulations/{batch}/data/ruptures.list',dtype=str)

# Read in arrival time dataframe
arriv_df = pd.read_csv(f'/Users/tnye/ONC/event_detection/p_waves/{batch}_parrivals.csv')

# Set up folder to save sta_lta csv files for different ranks (parallelized)
outdir = '/Users/tnye/ONC/event_detection/sta_lta/sta_lta_mpi_polar_test'
if not path.exists(outdir):
    makedirs(outdir)

# Set up folder for figures of outlier pick times (large discrepancy between
# STA/LTA pick and the true arrival time) 
figdir = '/Users/tnye/ONC/event_detection/plots/outliers/polarized'
if not path.exists(f'{figdir}/positive'):
    makedirs(f'{figdir}/positive')
if not path.exists(f'{figdir}/negative'):
    makedirs(f'{figdir}/negative')

# Are we using the polarization filter for the ONC parameter test?
polarize = True

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
obspy_triggered = []
obspy_res = []
obspy_pick = []
obspy_channel = []
elarmS_triggered = []
elarmS_res = []
elarmS_pick = []
elarmS_channel = []
mags = []
rhyp_list = []
    
for index in recvbuf:
    
    krup = int(index)
    rupture = ruptures[krup]
    
    run = rupture.replace('.rupt','')
    
    print(run)
    
    # Loop over station type
    for stn_type in stn_types:

        # Waveform directory
        wf_dir = f'{home_dir}/{stn_type}'
        
        # Gather waveforms
        bb_files_Z = np.array(sorted(glob(f"{wf_dir}/*SignalwithNoise/*{run.replace('.','-')}/*HNZ.mseed")))
        bb_files_N = np.array(sorted(glob(f"{wf_dir}/*SignalwithNoise/*{run.replace('.','-')}/*HNN.mseed")))
        bb_files_E = np.array(sorted(glob(f"{wf_dir}/*SignalwithNoise/*{run.replace('.','-')}/*HNE.mseed")))
        
        noise_files_Z = np.array(sorted(glob(f"{wf_dir}/*_Noise/*{run.replace('.','-')}/*HNZ.mseed")))
        noise_files_N = np.array(sorted(glob(f"{wf_dir}/*_Noise/*{run.replace('.','-')}/*HNN.mseed")))
        noise_files_E = np.array(sorted(glob(f"{wf_dir}/*_Noise/*{run.replace('.','-')}/*HNE.mseed")))
        
        # Loop over waveforms
        for i_file in range(len(bb_files_Z)):
            
            # Read in waveforms
            st_acc_unfilt_E = read(bb_files_E[i_file])
            st_acc_unfilt_N = read(bb_files_N[i_file])
            st_acc_unfilt_Z = read(bb_files_Z[i_file])
            dt = st_acc_unfilt_E[0].stats.delta
            stn =  st_acc_unfilt_E[0].stats.station.split('.')[0]
            
            # Extract event info from rupture .log file
            with open(f'/Users/tnye/ONC/simulations/{batch}/output/ruptures/{run}.log') as f:
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
            
            ## Test using elarmS (EPIC) picker
            components = ['N','E','Z']
            
            st = Stream()
            st.append(st_acc_unfilt_N[0])
            st.append(st_acc_unfilt_E[0])
            st.append(st_acc_unfilt_Z[0])
            
            st[0].stats.channel = 'HNN'
            st[1].stats.channel = 'HNE'
            st[2].stats.channel = 'HNZ'
            
            true_p = UTCDateTime(true_arriv)
            sig_ind = np.where(st[0].times('UTCDateTime') >= true_p)[0]
            
            elarms_pick_times = []
            trig_ind = []
            
            for i, tr in enumerate(st):
                D_series,V_series,A_series,N,S = picker.SUBROUTINE_TO_GROUND(tr)
                _,_,_,stalta = picker.SUBROUTINE_PICKER(V_series,tr)
                thresh_ind = np.where(np.array(stalta) >= 20)[0]
                
                try:
                    t_trig = tr.times('UTCDateTime')[np.intersect1d(sig_ind, thresh_ind)[0]]
                    elarms_pick_times.append(t_trig)
                    trig_ind.append(i)
                except:
                    continue
            
            if len(trig_ind) > 0:
                pick_ind = np.where(np.array(elarms_pick_times)==np.min(np.array(elarms_pick_times)))[0][0]
                elarmS_triggered.append('Yes')
                elarmS_pick.append(elarms_pick_times[pick_ind])
                elarmS_res.append(true_p - UTCDateTime(elarms_pick_times[pick_ind]))
                elarmS_channel.append(components[pick_ind])
            else:
                elarmS_triggered.append('No')
                elarmS_pick.append(np.nan)
                elarmS_res.append(np.nan)
                elarmS_channel.append(np.nan)
            
            # fig, axs = plt.subplots(2,1)
            # axs[0].plot(tr.times('UTCDateTime'),V_series,color='k',lw=1)
            # axs[0].vlines(UTCDateTime(elarms_pick_times[0]),np.min(V_series),np.max(V_series),color='red',label='P-arrival pick',lw=1)
            # axs[0].grid(linestyle='-',alpha=0.5)
            # axs[0].vlines(UTCDateTime(true_p),np.min(V_series),np.max(V_series),color='blue',label='True P-arrival',lw=1)
            # axs[0].set_xlabel('Time (s)')
            # axs[0].set_ylabel('Amplitude (m/s)')
            # axs[0].legend(facecolor='white')
            # axs[0].text(0.98,0.1,f'M{round(mag,1)}',ha='right',transform = axs[0].transAxes)
            # axs[0].text(0.98,0.02,f'Rhyp = {round(rhyp)} (km)',ha='right',transform = axs[0].transAxes)
            
            # axs[1].plot(st_vel_filt[0].times(),cft_vel,color='k')
            # axs[1].grid(linestyle='-',alpha=0.5)
            # axs[1].hlines(thresh, np.min(st_vel_filt[0].times()), np.max(st_vel_filt[0].times()),
            #               ls='--', lw=1, color='red',label='triggering threshold')
            # axs[1].set_xlabel('Time (s)')
            # axs[1].set_ylabel('STA/LTA')
            # axs[1].legend(facecolor='white')
            # fig.suptitle(f'{stn}')
            # plt.subplots_adjust(left=0.15,hspace=0.3,top=0.925)

            ###################################################################
            
            ## Test using Obspy with EPIC parameters
            sta, lta, thresh = [0.05,5,20]
            components = ['N','E','Z']
            
            # Highpass filter acceleration at 0.075 Hz
            st_acc_filt_E = emf.highpass(st_acc_unfilt_E,1/15,1/dt,4,zerophase=False)
            st_acc_filt_N = emf.highpass(st_acc_unfilt_N,1/15,1/dt,4,zerophase=False)
            st_acc_filt_Z = emf.highpass(st_acc_unfilt_Z,1/15,1/dt,4,zerophase=False)
            
            # Integrate to velocity and highpass filter again (EPIC only)
            st_vel_unfilt_E = emf.accel_to_veloc(st_acc_filt_E)
            st_vel_filt_E = emf.highpass(st_vel_unfilt_E,1/15,1/dt,4,zerophase=False)
            st_vel_unfilt_N = emf.accel_to_veloc(st_acc_filt_N)
            st_vel_filt_N = emf.highpass(st_vel_unfilt_N,1/15,1/dt,4,zerophase=False)
            st_vel_unfilt_Z = emf.accel_to_veloc(st_acc_filt_Z)
            st_vel_filt_Z = emf.highpass(st_vel_unfilt_Z,1/15,1/dt,4,zerophase=False)
        
            # Run obspy STA/LTA algorithm
            cft_E = classic_sta_lta(st_vel_filt_E[0].data, int(sta * 1/dt), int(lta * 1/dt))
            cft_N = classic_sta_lta(st_vel_filt_N[0].data, int(sta * 1/dt), int(lta * 1/dt))
            cft_Z = classic_sta_lta(st_vel_filt_Z[0].data, int(sta * 1/dt), int(lta * 1/dt))
            
            # Initialize lists
            trig_ind = []
            obspy_pick_times = []
            
            # Get indices of where STA/LTA passed threshold 
            thresh_ind_N = np.where((np.array(cft_N) >= thresh) & (st_acc_filt_E[0].times('matplotlib') > (true_p - 5)))[0]
            thresh_ind_E = np.where((np.array(cft_E) >= thresh) & (st_acc_filt_E[0].times('matplotlib') > (true_p - 5)))[0]
            thresh_ind_Z = np.where((np.array(cft_Z) >= thresh) & (st_acc_filt_E[0].times('matplotlib') > (true_p - 5)))[0]
            
            # Check for each component if algorithm was triggered
            if len(np.where(sig_ind == thresh_ind_N)[0])>0:
                trig_ind.append(0)
                obspy_pick_times.append(st_vel_filt_N[0].times('UTCDateTime')[np.where(sig_ind == thresh_ind_N)[0][0]])
            if len(np.where(sig_ind == thresh_ind_E)[0])>0:
                trig_ind.append(1)
                obspy_pick_times.append(st_vel_filt_E[0].times('UTCDateTime')[np.where(sig_ind == thresh_ind_E)[0][0]])
            if len(np.where(sig_ind == thresh_ind_Z)[0])>0:
                trig_ind.append(2)
                obspy_pick_times.append(st_vel_filt_Z[0].times('UTCDateTime')[np.where(sig_ind == thresh_ind_Z)[0][0]])
            
            # If triggered, add pick time and residuals to lists
            if len(trig_ind) > 0:
                pick_ind = np.where(np.array(obspy_pick_times)==np.min(np.array(obspy_pick_times)))[0][0]
                obspy_triggered.append('Yes')
                obspy_pick.append(obspy_pick_times[pick_ind])
                obspy_res.append(UTCDateTime(true_arriv) - obspy_pick_times[pick_ind])
                obspy_channel.append(components[trig_ind[pick_ind]])
            else:
                obspy_triggered.append('No')
                obspy_pick.append(np.nan)
                obspy_res.append(np.nan)
                obspy_channel.append(np.nan)

                
            # # Plot outlier
            # if np.abs(UTCDateTime(true_arriv) - pick)>=75:
                
            #     if (UTCDateTime(true_arriv) - pick)>0:
            #         folder = 'positive'
            #     else:
            #         folder = 'negative'
            #     fig, axs = plt.subplots(2,1)
            #     axs[0].plot(st_acc_filt[0].times(),st_vel_filt[0].data,color='k')
            #     if len(np.where(cft_vel>thresh)[0]>0):
            #         p_ind = np.where(cft_vel>thresh)[0][0]
            #         axs[0].vlines(st_vel_filt[0].times()[p_ind],np.min(st_vel_filt[0].data),np.max(st_vel_filt[0].data),color='red',label='P-arrival pick',lw=1)
            #     axs[0].grid(linestyle='-',alpha=0.5)
            #     axs[0].vlines(UTCDateTime(true_arriv)-UTCDateTime(st_vel_filt[0].stats.starttime),np.min(st_vel_filt[0].data),np.max(st_vel_filt[0].data),color='blue',label='True P-arrival',lw=1)
            #     axs[0].set_xlabel('Time (s)')
            #     axs[0].set_ylabel('Amplitude (m/s)')
            #     axs[0].legend(facecolor='white')
            #     axs[0].text(0.98,0.1,f'M{round(mag,1)}',ha='right',transform = axs[0].transAxes)
            #     axs[0].text(0.98,0.02,f'Rhyp = {round(rhyp)} (km)',ha='right',transform = axs[0].transAxes)
                
            #     axs[1].plot(st_vel_filt[0].times(),cft_vel,color='k')
            #     axs[1].grid(linestyle='-',alpha=0.5)
            #     axs[1].hlines(thresh, np.min(st_vel_filt[0].times()), np.max(st_vel_filt[0].times()),
            #                   ls='--', lw=1, color='red',label='triggering threshold')
            #     axs[1].set_xlabel('Time (s)')
            #     axs[1].set_ylabel('STA/LTA')
            #     axs[1].legend(facecolor='white')
            #     fig.suptitle(f'{stn}')
            #     plt.subplots_adjust(left=0.15,hspace=0.3,top=0.925)
               
                
            #     plt.savefig(f'{figdir}/{folder}/{run}_{stn}.{channel}_{sta}-{lta}-{thresh}_vel.png',dpi=300)
            #     plt.close()
            
            
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
            
            # Determine if algorithm was triggered
            if len(np.where(sig_ind == thresh_ind)[0])>0:
                pick = st_Z_filt[0].times('UTCDateTime')[np.where(sig_ind == thresh_ind)[0][0]]
                onc_triggered.append('Yes')
                onc_pick.append(pick)
                onc_res.append(UTCDateTime(true_arriv) - pick)
                
                # # Plot outlier
                # if np.abs(UTCDateTime(true_arriv) - pick)>=75:
                    
                #     if (UTCDateTime(true_arriv) - pick)>0:
                #         folder = 'positive'
                #     else:
                #         folder = 'negative'
                #     fig, axs = plt.subplots(2,1)
                #     axs[0].plot(st_Z_filtfilt[0].times(),data,color='k')
                #     if len(np.where(cft>thresh)[0]>0):
                #         p_ind = np.where(cft>thresh)[0][0]
                #         axs[0].vlines(st_Z_filtfilt[0].times()[p_ind],np.min(data),np.max(data),color='red',label='P-arrival pick',lw=1)
                #     axs[0].grid(linestyle='-',alpha=0.5)
                #     axs[0].vlines(UTCDateTime(true_arriv)-UTCDateTime(st_Z_filtfilt[0].stats.starttime),np.min(data),np.max(data),color='blue',label='True P-arrival',lw=1)
                #     axs[0].set_xlabel('Time (s)')
                #     axs[0].set_ylabel(r'Amplitude (m/s$^2$)')
                #     axs[0].legend(facecolor='white')
                #     axs[0].text(0.98,0.1,f'M{round(mag,1)}',ha='right',transform = axs[0].transAxes)
                #     axs[0].text(0.98,0.02,f'Rhyp = {round(rhyp)} (km)',ha='right',transform = axs[0].transAxes)
                    
                #     axs[1].plot(st_acc_filt[0].times(),cft,color='k')
                #     axs[1].grid(linestyle='-',alpha=0.5)
                #     axs[1].hlines(thresh, np.min(st_acc_filt_filt[0].times()), np.max(st_acc_filt_filt[0].times()),
                #                   ls='--', lw=1, color='red',label='triggering threshold')
                #     axs[1].set_xlabel('Time (s)')
                #     axs[1].set_ylabel('STA/LTA')
                #     axs[1].legend(facecolor='white')
                #     fig.suptitle(f'{stn}')
                #     plt.subplots_adjust(left=0.15,hspace=0.3,top=0.925)
                   
                #     plt.savefig(f'{figdir}/{folder}/{run}_{stn}.{channel}_{sta}-{lta}-{thresh}_acc.png',dpi=300)
                #     plt.close()
            else:
                onc_triggered.append('No')
                onc_pick.append(np.nan)
                onc_res.append(np.nan)
            
    
# Make dataframe
data = {'Batch':batches,'Run':rupts,'Magnitude':mags,'Station':stns,
        '1/5/4.9_triggered':onc_triggered,'1/5/4.9_P-pick time':onc_pick,
        '1/5/4.9_P-arrival Residual (s)':onc_res,'0.05/5/20_obspy_triggered':obspy_triggered,
        '0.05/5/20_obspy_P-pick time':obspy_pick,'0.05/5/20_obspy_P-arrival Residual (s)':obspy_res,
        '0.05/5/20_obspy_channel':obspy_channel,'0.05/5/20_elarmS_triggered':elarmS_triggered,
        '0.05/5/20_elarmS_P-pick time':elarmS_pick,'0.05/5/20_elarmS_P-arrival Residual (s)':elarmS_res,
        '0.05/5/20_elarmS_channel':elarmS_channel,}

df = pd.DataFrame(data)
df.to_csv(f'{outdir}/{batch}_sta_lta_{rank}.csv')


    
