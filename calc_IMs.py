#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 04:30:50 2022

@author: tnye
"""

###############################################################################
# Script that goes through synethic waveforms generated by FakeQuakes,
# calculates IMs, and computes the residuals between GMM IM and synthetic wf IM.
###############################################################################

# Standard library imports 
from glob import glob
from os import makedirs, path
import sys
import numpy as np
import pandas as pd
from obspy import read
from mpi4py import MPI

# Local imports
import single_shakemaps
sys.path.append('/Users/tnye/tsuquakes/code')
import tsueqs_main_fns as tmf
sys.path.append('/Users/tnye/tsuquakes/code/processing')
import signal_average_fns as avg
import valid_fns as valid
sys.path.append('/Users/tnye/kappa/code/github_code/processing')
from rotd50 import compute_rotd50

min_mag = 0
min_lat = 0
max_lat = 90


###################### Set up parallelization parameters ######################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

ncpus = size

# Names of batch and run
batch = 'cascadia_longer_wfs'
# stn_type = 'onc-onshore'
# wf_folder = 'onc'
# stn_type = 'pnsn'
# wf_folder = 'pnsn'
stn_type = 'onc-offshore'
wf_folder = 'onc'
# stn_types = ['lf','hf']
stn_types = ['hf']

# Read in station df (contains station names and coordinates)
stn_df = pd.read_csv(f'/Users/tnye/ONC/data/station_info/{stn_type}_list.txt', delimiter=',')

vs30_csv = f'/Users/tnye/ONC/data/files/vs30_coords_{stn_type}.txt'

ruptures = np.genfromtxt(f'/Users/tnye/ONC/simulations/{batch}/data/ruptures.list',dtype=str)

fcorner_high = 1/15.   
order = 2

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

for index in recvbuf:

    rupture = ruptures[int(index)]
    run = rupture.replace('.rupt','')
    
    logfile = f"/Users/tnye/ONC/simulations/{batch}/output/ruptures/{rupture.replace('.rupt','.log')}"
    
    f = open(logfile, 'r')
    lines = f.readlines()
    mag = float(lines[15].split(' ')[-1].split('\n')[0])
    lat = float(lines[16].split('(')[-1].split(',')[1])
    
    if mag >= min_mag and lat <= max_lat and lat >= min_lat:
    
        print(f'Working on rupture {run}')
    
        # Path to save output csv
        flatfile_path = f'/Users/tnye/ONC/flatfiles/IMs/{batch}/{run}_{stn_type}.csv'
        
        # Synthetic miniseed directories
        wf_dir = f'/Users/tnye/ONC/simulations/{batch}/waveforms_data_curation'
        
        # Gather waveforms
        if stn_type == 'onc-offshore':
            bb_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-offshore-strongmotion/signal/{run}/*.mseed')))
            bb_noise_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-offshore-strongmotion/signal_with_noise/{run}/*.mseed')))
            # lf_files = np.array(sorted(glob(f'{wf_dir}/ONC-offshore/{run}/' + '*.W1.LY*.sac')))
        elif stn_type == 'onc-onshore':
            bb_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-onshore-strongmotion/signal/{run}/*.mseed')))
            bb_noise_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-onshore-strongmotion/signal_with_noise/{run}/*.mseed')))
            lf_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-onshore-gnss/signal/{run}/' + '*.LY*.mseed')))
            lf_noise_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-onshore-gnss/signal_with_noise/{run}/' + '*.mseed')))
        else:
            bb_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-onshore-strongmotion/signal/{run}/' + '*.mseed')))
            bb_noise_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-onshore-strongmotion/signal_with_noise/{run}/' + '*.mseed')))
            lf_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-onshore-lowfreq/signal/{run}/' + '*.mseed')))
            lf_noise_files = np.array(sorted(glob(f'{wf_dir}/{wf_folder}-onshore-lowfreq/signal_with_noise/{run}/' + '*.mseed')))
            
        # Get rupture files
        rupt_file = f'/Users/tnye/ONC/simulations/{batch}/output/ruptures/{run}.rupt'
        log_file = f'/Users/tnye/ONC/simulations/{batch}/output/ruptures/{run}.log'
        
        # Home directory 
        home_dir = '/Users/tnye/ONC'
        
        # Read in station metadata file
        stn_metadata = pd.read_csv('/Users/tnye/ONC/data/station_info/all_stations_metadata.csv')
        # stn_metadata = pd.read_csv('/Users/tnye/ONC/data/station_info/pnsn_metadata.csv')
        
        # Create shakefiles
        single_shakemaps.run_shakemaps(batch, rupt_file, vs30_csv, stn_type)
    
        # Read in dfs for BCHydro and NGA-Sub and extract IM info
        shake_file = f'/Users/tnye/ONC/shakemaps/shakefiles/{batch}/{stn_type}/{run}.shake'
        shake_df = pd.read_csv(shake_file,delimiter=',')
        
        # Get GMM predictions
        BCHydro_pga = np.array(shake_df['PGA_bchydro2018_g'])
        BCHydro_sd = np.exp(np.array(shake_df['PGAsd_bchydro2018_lng']))
        BCHydro_MMI_pga = np.array(shake_df['MMI_bchydro_wgrw12_pga'])
        NGASub_pga = np.array(shake_df['PGA_NGASub_g'])
        NGASub_pga_sd = np.exp(np.array(shake_df['PGAsd_NGASub_lng']))
        NGASub_MMI_pga = np.array(shake_df['MMI_NGASub_wgrw12_pga'])
        NGASub_MMI_pgv = np.array(shake_df['MMI_NGASub_wgrw12_pgv'])
        NGASub_pgv = np.array(shake_df['PGV_NGASub_cm/s'])
        NGASub_pgv_sd = np.exp(np.array(shake_df['PGVsd_NGASub_lncm/s']))
        
        # Set up folder for flatfile
        if not path.exists(f'{home_dir}/flatfiles/IMs/{batch}'):
            makedirs(f'{home_dir}/flatfiles/IMs/{batch}')
        
        
        ########################## Set up initial parameters ##########################
        
        # Extract event info from rupture .log file
        with open(log_file) as f:
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
            mw = float(mag_line.split(':')[1].split(' ')[2])
            m0 = 10**(mw*(3/2.) + 9.1)
            
            # Get origintime
            orig_line = lines[17]
            orig_line = orig_line.replace('\n', '')
            origintime = orig_line.split(': ')[1].split('.')[0]
        
        # Create lists for of the event and station info for the df
        eventnames = np.array([])
        origintimes = np.array([])
        hyplons = np.array([])
        hyplats = np.array([])
        hypdepths = np.array([])
        mws = np.array([])
        m0s = np.array([])
        networks = np.array([])
        stations = np.array([])
        stlons = np.array([])
        stlats = np.array([])
        stelevs = np.array([])
        hypdists = np.array([])
        epidists = np.array([])
        instrument_codes = np.array([])
        pga_list = np.array([])
        pga_noise_list = np.array([])
        pgv_list = np.array([])
        pgv_noise_list = np.array([])
        pgd_list = np.array([])
        pgd_noise_list = np.array([])
        melgar_pgd = np.array([])
        melgar_sd = np.array([])
        gold_sd = np.array([])
        gold_pgd = np.array([])
        MMI_pga_list = np.array([])
        MMI_pgv_list = np.array([])
        MMI_noise_pga_list = np.array([])
        MMI_noise_pgv_list = np.array([])
        
        rrups = np.array(shake_df['Rrupt(km)'])
        
        # Determine station column name in file
        try:
            stn_df['#Station'].values
            stn_header = '#Station'
        except:
            stn_header = '#station'
        
        ####################### High Frequency Intensity Measures #####################
        
        if 'hf' in stn_types:
            
            # Group channel files by station
            N = 3
            stn_files = [bb_files[n:n+N] for n in range(0, len(bb_files), N)]
            stn_noise_files = [bb_noise_files[n:n+N] for n in range(0, len(bb_noise_files), N)]
            
            # Loop over files to get the list of station names, channels, and mseed files 
            for i, group in enumerate(stn_files):
                        
                stn = group[0].split('/')[-1].split('.')[0]
                
                if stn in stn_df[stn_header].values:
                
                    components = np.array([])
                    
                    for mseed in group:
                        # Get the instrument component (E,N,Z) for this station
                        channel_code = mseed.split('/')[-1].split('bb')[-1].split('.')[-2]
                        components = np.append(components,channel_code[2])
                
                    # Get station coodinates
                    stn_id = np.where(np.array(stn_df.iloc[:,0])==stn)[0][0]
                    stlon = stn_df.iloc[stn_id,1]
                    stlat = stn_df.iloc[stn_id,2]
                    stelev = float(stn_metadata.loc[stn_id].Elevation)    
                    
                    # Compute hypocentral distance
                    hypdist = tmf.compute_rhyp(stlon,stlat,stelev,hyplon,hyplat,hypdepth)
                    hypdists = np.append(hypdists,hypdist)
                    
                    # Compute epicentral distance
                    repi = tmf.compute_repi(stlon,stlat,hyplon,hyplat)
                    epidists = np.append(epidists, repi)
                
                    # Append the earthquake and station info for this station
                    origintimes = np.append(origintimes,origintime)
                    hyplons = np.append(hyplons,hyplon)
                    hyplats = np.append(hyplats,hyplat)
                    hypdepths = np.append(hypdepths,hypdepth)
                    mws = np.append(mws,mw)
                    m0s = np.append(m0s,m0)
                    stations = np.append(stations,stn)
                    stlons = np.append(stlons,stlon)
                    stlats = np.append(stlats,stlat)
                
                
                    ############################ Read in Waveforms ############################
                
                    # Acceleration 
                    E_index = np.where(components=='E')[0][0]
                    E_record = read(group[E_index])
                
                    N_index = np.where(components=='N')[0][0]
                    N_record = read(group[N_index])
                
                    Z_index = np.where(components=='Z')[0][0]                       
                    Z_record = read(group[Z_index])
                    
                    # Acceleration w/ noise
                    E_record_noise = read(stn_noise_files[i][E_index])
                    N_record_noise = read(stn_noise_files[i][N_index])
                    Z_record_noise = read(stn_noise_files[i][Z_index])
                    
                    stsamprate = E_record[0].stats.sampling_rate
                    
                    # Convert acceleration record to velocity (HP filter)
                    E_vel = tmf.accel_to_veloc(E_record)
                    N_vel = tmf.accel_to_veloc(N_record)
                    Z_vel = tmf.accel_to_veloc(Z_record)
                    E_vel_noise = tmf.accel_to_veloc(E_record_noise)
                    N_vel_noise = tmf.accel_to_veloc(N_record_noise)
                    Z_vel_noise = tmf.accel_to_veloc(Z_record_noise)
                    
                    # High pass filter velocity 
                    E_vel = tmf.highpass(E_vel,fcorner_high,stsamprate,order,zerophase=True)
                    N_vel = tmf.highpass(N_vel,fcorner_high,stsamprate,order,zerophase=True)
                    Z_vel = tmf.highpass(Z_vel,fcorner_high,stsamprate,order,zerophase=True)
                    E_vel_noise = tmf.highpass(E_vel_noise,fcorner_high,stsamprate,order,zerophase=True)
                    N_vel_noise = tmf.highpass(N_vel_noise,fcorner_high,stsamprate,order,zerophase=True)
                    Z_vel_noise = tmf.highpass(Z_vel_noise,fcorner_high,stsamprate,order,zerophase=True)
                    
                    # Save vel as sac files
                    if not path.exists(f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal/{run}'):
                        makedirs(f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal/{run}')
                    if not path.exists(f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal_with_noise/{run}'):
                        makedirs(f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal_with_noise/{run}')
                      
                    E_vel[0].stats.channel = 'HNE'
                    E_vel_filename = f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal/{run}/{stn}.vel.{E_vel[0].stats.channel}.mseed' 
                    E_vel[0].write(E_vel_filename, format='MSEED')
                    
                    N_vel[0].stats.channel = 'HNN'
                    N_vel_filename = f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal/{run}/{stn}.vel.{N_vel[0].stats.channel}.mseed' 
                    N_vel[0].write(N_vel_filename, format='MSEED')
                    
                    Z_vel[0].stats.channel = 'HNZ'
                    Z_vel_filename = f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal/{run}/{stn}.vel.{Z_vel[0].stats.channel}.mseed' 
                    Z_vel[0].write(Z_vel_filename, format='MSEED')
                    
                    E_vel_noise[0].stats.channel = 'HNE'
                    E_vel_filename = f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal_with_noise/{run}/{stn}.vel.{E_vel[0].stats.channel}.mseed' 
                    E_vel_noise[0].write(E_vel_filename, format='MSEED')
                    
                    N_vel_noise[0].stats.channel = 'HNN'
                    N_vel_filename = f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal_with_noise/{run}/{stn}.vel.{N_vel[0].stats.channel}.mseed' 
                    N_vel_noise[0].write(N_vel_filename, format='MSEED')
                    
                    Z_vel_noise[0].stats.channel = 'HNZ'
                    Z_vel_filename = f'{home_dir}/simulations/{batch}/vel_waveforms/{stn_type}/signal_with_noise/{run}/{stn}.vel.{Z_vel[0].stats.channel}.mseed' 
                    Z_vel_noise[0].write(Z_vel_filename, format='MSEED')
                    
                
                    ####################### Calculate Intensity Measures ######################
                    
                    # Get rotd50 of acceleration components 
                    acc_avg = compute_rotd50(E_record[0].data,N_record[0].data)
                    acc_avg_noise = compute_rotd50(E_record_noise[0].data,N_record_noise[0].data)
                   
                    # Calculate PGA
                    pga_m_s2 = np.max(np.abs(acc_avg))
                    pga_g = pga_m_s2/9.81
                    pga_list = np.append(pga_list,pga_g)
                    
                    pga_noise_m_s2 = np.max(np.abs(acc_avg_noise))
                    pga_noise_g = pga_noise_m_s2/9.81
                    pga_noise_list = np.append(pga_noise_list,pga_noise_g)
                    
                    # Get rotd50 of velocity components 
                    vel_avg = compute_rotd50(E_vel[0].data,N_vel[0].data)
                    vel_avg_noise = compute_rotd50(E_vel_noise[0].data,N_vel_noise[0].data)
                    
                    # Calculate PGV
                    pgv_cm_s = np.max(np.abs(vel_avg))*100
                    pgv_list = np.append(pgv_list,pgv_cm_s)
                    
                    pgv_noise_cm_s = np.max(np.abs(vel_avg_noise))*100
                    pgv_noise_list = np.append(pgv_noise_list,pgv_noise_cm_s)
                
                    # Calculate MMI
                    MMI_pga = valid.WGRW12(np.array([pga_m_s2*100]),0)
                    MMI_pga_list = np.append(MMI_pga_list,MMI_pga)
                    
                    MMI_pgv = valid.WGRW12(np.array([pgv_cm_s]),1)
                    MMI_pgv_list = np.append(MMI_pgv_list,MMI_pgv)
                    
                    MMI_noise_pga = valid.WGRW12(np.array([pga_noise_m_s2*100]),0)
                    MMI_noise_pga_list = np.append(MMI_noise_pga_list,MMI_noise_pga)
                    
                    MMI_noise_pgv = valid.WGRW12(np.array([pgv_noise_cm_s*100]),1)
                    MMI_noise_pgv_list = np.append(MMI_noise_pgv_list,MMI_noise_pgv)
                   
                    
        ####################### Low Freqency Intensity Measures #####################
        
        if 'lf' in stn_types:
       
            # Group all files by station
            N = 3
            stn_files = [lf_files[n:n+N] for n in range(0, len(lf_files), N)]
            stn_noise_files = [lf_noise_files[n:n+N] for n in range(0, len(lf_noise_files), N)]
    
            # Loop over files to get the list of station names, channels, and mseed files 
            for i, group in enumerate(stn_files):
                
                stn = group[0].split('/')[-1].split('.')[0]
                if stn in stn_df[stn_header].values:
                    components = np.array([])
                    
                    for mseed in group:
                        # Get the instrument component (E,N,Z) for this station
                        if 'W1' in mseed:
                            channel_code = mseed.split('/')[-1].split('.')[2]
                        else:
                            channel_code = mseed.split('/')[-1].split('.')[1]
                        components = np.append(components,channel_code[2])
                
                    # Get station coodinates
                    stn_id = np.where(np.array(stn_df.iloc[:,0])==stn)[0][0]
                    stlon = stn_df.iloc[stn_id,1]
                    stlat = stn_df.iloc[stn_id,2]
                    stelev = float(stn_metadata.loc[stn_id].Elevation)
                
                    ##################### Start computations ######################        
                    
                    # Compute hypocentral distance
                    hypdist = tmf.compute_rhyp(stlon,stlat,stelev,hyplon,hyplat,hypdepth)
                    
                    if 'hf' not in stn_types:
                        hypdists = np.append(hypdists,hypdist)
                        
                        # Compute epicentral distance
                        repi = tmf.compute_repi(stlon,stlat,hyplon,hyplat)
                        epidists = np.append(epidists, repi)
                    
                        # Append the earthquake and station info for this station
                        origintimes = np.append(origintimes,origintime)
                        hyplons = np.append(hyplons,hyplon)
                        hyplats = np.append(hyplats,hyplat)
                        hypdepths = np.append(hypdepths,hypdepth)
                        mws = np.append(mws,mw)
                        m0s = np.append(m0s,m0)
                        stations = np.append(stations,stn)
                        stlons = np.append(stlons,stlon)
                        stlats = np.append(stlats,stlat)
                        
                    
                    ############################ Read in Waveforms ############################
                
                    # Displacement
                    E_index = np.where(components=='E')[0][0]
                    E_record = read(group[E_index])
                
                    N_index = np.where(components=='N')[0][0]
                    N_record = read(group[N_index])
                
                    Z_index = np.where(components=='Z')[0][0]                       
                    Z_record = read(group[Z_index])
                    
                    # Calcualte noisy GNSS data intensity measures
                    E_noise_record = read(stn_noise_files[i][E_index])
                    N_noise_record = read(stn_noise_files[i][N_index])
                    Z_noise_record = read(stn_noise_files[i][Z_index])
                
                
                    ####################### Calculate Intensity Measures ######################
                    
                    # Get euclidean norm of displacement components 
                    disp_avg = avg.get_eucl_norm_3comp(E_record[0].data, N_record[0].data, Z_record[0].data)
                    disp_noise_avg = avg.get_eucl_norm_3comp(E_noise_record[0].data, N_noise_record[0].data, Z_noise_record[0].data)
                    
                    # Calculate PGD
                    pgd = np.max(np.abs(disp_avg))
                    pgd_list = np.append(pgd_list,pgd)
                    
                    disp_noise_avg -= np.mean(disp_noise_avg[:61])
                    pgd_noise = np.max(np.abs(disp_noise_avg))
                    pgd_noise_list = np.append(pgd_noise_list,pgd_noise)     
                    
                    # Estimate PGD from Melgar et al., 2015 model
                    pgd, sd = valid.get_pgd_scaling(mw, hypdist,'MA15')
                    melgar_pgd = np.append(melgar_pgd,pgd)
                    melgar_sd = np.append(melgar_sd,sd)
                    
                    # Estimate PGD from Goldberg et al., 2021 model
                    Rp = valid.get_Rp(stlon,stlat,stelev,rupt_file)
                    pgd, sd = valid.get_pgd_scaling(mw, Rp,'GA21')
                    gold_pgd = np.append(gold_pgd,pgd)
                    gold_sd = np.append(gold_sd,sd)
        
        if 'lf' in stn_types:
            lnPGD_res_melgar = np.log(melgar_pgd) - np.log(pgd_list)
            lnPGD_res_gold = np.log(gold_pgd) - np.log(pgd_list)
            lnPGD_noise_res_melgar = np.log(melgar_pgd) - np.log(pgd_noise_list)
            lnPGD_noise_res_gold = np.log(gold_pgd) - np.log(pgd_noise_list)
        
        if 'hf' in stn_types:
            lnPGA_res_bc = np.log(BCHydro_pga) - np.log(pga_list)
            lnPGA_noise_res_bc = np.log(BCHydro_pga) - np.log(pga_noise_list)
            MMI_res_bc_pga = BCHydro_MMI_pga - MMI_pga_list
            MMI_noise_res_bc_pga = BCHydro_MMI_pga - MMI_noise_pga_list
            
            lnPGA_res_nga = np.log(NGASub_pga) - np.log(pga_list)
            lnPGA_noise_res_nga = np.log(NGASub_pga) - np.log(pga_noise_list)
            MMI_res_nga_pga = NGASub_MMI_pga - MMI_pga_list
            MMI_res_nga_pgv = NGASub_MMI_pgv - MMI_pgv_list
            MMI_noise_res_nga_pga = NGASub_MMI_pga - MMI_noise_pga_list
            MMI_noise_res_nga_pgv = NGASub_MMI_pgv - MMI_noise_pgv_list
            lnPGV_res_nga = np.log(NGASub_pgv) - np.log(pgv_list)
            lnPGV_noise_res_nga = np.log(NGASub_pgv) - np.log(pgv_noise_list)

        
        ############################### Build dataframe ###############################
        
        if 'lf' in stn_types and 'hf' in stn_types:
            dataset_dict = {'Station':stations,'stlon':stlons,'stlat':stlats,
                            'Mw':mws,'M0':m0s,'hyplon':hyplons,'hyplat':hyplats,
                            'hypdepth(km)':hypdepths,'origintime':origintimes,
                            'Rhyp(km)':hypdists,'Repi(km)':epidists,'Rrup(km)':rrups,
                            'syn_PGD_m':pgd_list,'syn_PGD_w/noise_m':pgd_noise_list,
                            'MA15_pgd_m':melgar_pgd,'MA15_sd_m':melgar_sd,
                            'GA21_pgd_m':gold_pgd,'GA21_sd_m':gold_sd,
                            'lnPGD_res_MA15':lnPGD_res_melgar,'lnPGD_noise_res_MA15':lnPGD_noise_res_melgar,
                            'lnPGD_res_GA21':lnPGD_res_gold,'lnPGD_noise_res_GA21':lnPGD_noise_res_gold,
                            'syn_PGA_g':pga_list,'syn_PGA_w/noise_g':pga_noise_list,
                            'BCHydro_pga_g':BCHydro_pga,'BCHydro_sd_g':BCHydro_sd,
                            'NGASub_pga_g':NGASub_pga,'NGA_Sub_sd_g':NGASub_pga_sd,
                            'lnPGA_res_BCHydro':lnPGA_res_bc,'lnPGA_res_w/noise_BCHydro':lnPGA_noise_res_bc,
                            'lnPGA_res_NGASub':lnPGA_res_nga,'lnPGA_res_w/noise_NGASub':lnPGA_noise_res_nga,
                            'syn_PGV_cm/s':pgv_list,'syn_PGV_w/noise_cm/s':pgv_noise_list,
                            'NGA_Sub_pgv_cm/s':NGASub_pgv,'NGA_Sub_sd_cm/s':NGASub_pgv_sd,
                            'lnPGV_res_NGASub':lnPGV_res_nga,'lnPGV_res_w/noise_NGASub':lnPGV_noise_res_nga,
                            'syn_MMI_pga':MMI_pga_list,'syn_MMI_w/noise_pga':MMI_noise_pga_list,
                            'syn_MMI_pgv':MMI_pgv_list,'syn_MMI_w/noise_pgv':MMI_noise_pgv_list,
                            'BCHydro_MMI_pga':BCHydro_MMI_pga,
                            'NGA_MMI_pga':NGASub_MMI_pga,'NGA_Sub_MMI_pgv':NGASub_MMI_pgv,
                            'MMI_res_BCHydro_pga':MMI_res_bc_pga,'MMI_res_w/noise_BCHydro_pga':MMI_noise_res_bc_pga,
                            'MMI_res_NGASub_pga':MMI_res_nga_pga,'MMI_res_w/noise_NGASub_pga':MMI_noise_res_nga_pga,
                            'MMI_res_NGASub_pgv':MMI_res_nga_pgv,'MMI_res_w/noise_NGASub_pgv':MMI_noise_res_nga_pgv}
                            
            # Make main dataframe using dictionary 
            df = pd.DataFrame(data=dataset_dict)
            
            # Save to flatfile:
            df.to_csv(flatfile_path,index=False)
        
        elif 'lf' in stn_types and 'hf' not in stn_types:
            dataset_dict = {'Station':stations,'stlon':stlons,'stlat':stlats,
                            'Mw':mws,'M0':m0s,'hyplon':hyplons,'hyplat':hyplats,
                            'hypdepth(km)':hypdepths,'origintime':origintimes,
                            'Rhyp(km)':hypdists,'Repi(km)':epidists,'Rrup(km)':rrups,
                            'syn_PGD_m':pgd_list,'syn_PGD_w/noise_m':pgd_noise_list,
                            'MA15_pgd_m':melgar_pgd,'MA15_sd_m':melgar_sd,
                            'GA21_pgd_m':gold_pgd,'GA21_sd_m':gold_sd,
                            'lnPGD_res_MA15':lnPGD_res_melgar,'lnPGD_noise_res_MA15':lnPGD_noise_res_melgar,
                            'lnPGD_res_GA21':lnPGD_res_gold,'lnPGD_noise_res_GA21':lnPGD_noise_res_gold}
        
            # Make main dataframe using dictionary 
            df = pd.DataFrame(data=dataset_dict)
            
            # Save to flatfile:
            df.to_csv(flatfile_path,index=False)
        
        elif 'hf' in stn_types and 'lf' not in stn_types:
            dataset_dict = {'Station':stations,'stlon':stlons,'stlat':stlats,
                            'Mw':mws,'M0':m0s,'hyplon':hyplons,'hyplat':hyplats,
                            'hypdepth(km)':hypdepths,'origintime':origintimes,
                            'Rhyp(km)':hypdists,'Repi(km)':epidists,'Rrup(km)':rrups,
                            'syn_PGA_g':pga_list,'syn_PGA_w/noise_g':pga_noise_list,
                            'BCHydro_pga_g':BCHydro_pga,'BCHydro_sd_g':BCHydro_sd,
                            'NGASub_pga_g':NGASub_pga,'NGA_Sub_sd_g':NGASub_pga_sd,
                            'lnPGA_res_BCHydro':lnPGA_res_bc,'lnPGA_res_w/noise_BCHydro':lnPGA_noise_res_bc,
                            'lnPGA_res_NGASub':lnPGA_res_nga,'lnPGA_res_w/noise_NGASub':lnPGA_noise_res_nga,
                            'syn_PGV_cm/s':pgv_list,'syn_PGV_w/noise_cm/s':pgv_noise_list,
                            'NGA_Sub_pgv_cm/s':NGASub_pgv,'NGA_Sub_sd_cm/s':NGASub_pgv_sd,
                            'lnPGV_res_NGASub':lnPGV_res_nga,'lnPGV_res_w/noise_NGASub':lnPGV_noise_res_nga,
                            'syn_MMI_pga':MMI_pga_list,'syn_MMI_w/noise_pga':MMI_noise_pga_list,
                            'syn_MMI_pgv':MMI_pgv_list,'syn_MMI_w/noise_pgv':MMI_noise_pgv_list,
                            'BCHydro_MMI_pga':BCHydro_MMI_pga,
                            'NGA_MMI_pga':NGASub_MMI_pga,'NGA_Sub_MMI_pgv':NGASub_MMI_pgv,
                            'MMI_res_BCHydro_pga':MMI_res_bc_pga,'MMI_res_w/noise_BCHydro_pga':MMI_noise_res_bc_pga,
                            'MMI_res_NGASub_pga':MMI_res_nga_pga,'MMI_res_w/noise_NGASub_pga':MMI_noise_res_nga_pga,
                            'MMI_res_NGASub_pgv':MMI_res_nga_pgv,'MMI_res_w/noise_NGASub_pgv':MMI_noise_res_nga_pgv}
                        
            # Make main dataframe using dictionary 
            df = pd.DataFrame(data=dataset_dict)
            
            # Save to flatfile:
            df.to_csv(flatfile_path,index=False)

