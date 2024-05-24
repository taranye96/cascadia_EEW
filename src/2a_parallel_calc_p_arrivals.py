#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 15:08:08 2021

@author: sydneydybing
modified by: tnye
"""

###############################################################################
# Script that calculates the true P-wave arrival times on the synthetics using
# the velocity model used to generate the sythetic data.
###############################################################################

# Imports
import numpy as np
import pandas as pd
from os import path,makedirs
from obspy.core import UTCDateTime
from obspy import read
from mpi4py import MPI
from numpy import genfromtxt,where,arange,ones,zeros,array,tile,argmin
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from glob import glob

# Local imports
import cascadiaEEW_main_fns as emf

###################### Set up parallelization parameters ######################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

ncpus = size

# Name of batch for simulations
batch = 'cascadia'

# Station types (used in names of waveform folders)
station_types = ['onc_onshore','onc_offshore','pnsn']

# Gather .log files
ruptures = sorted(glob(f'/home/tnye/onc/simulations/{batch}/output/ruptures/*.log'))

# Read in station metadata
metadata = pd.read_csv('/home/tnye/onc/data/station_info/station_metadata.csv')

# Build velocity model
velmod = TauPyModel(model=f'/home/tnye/onc/simulations/{batch}/structure/cascadia.npz') # velocity model


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

# Initialize file for saving P-wave arrival times
w = open('/Users/tnye/ONC/event_detection/p_waves/'+batch.split('/')[-1]+f'_arrivals_{rank}.csv', 'w')
w.write('Batch,Run,Station,Rhyp(km),Repi(km),P-wave arriv,S-wave arriv\n')

for index in recvbuf:

    file = ruptures[int(index)]
    
    # Get run number
    run = file.split('/')[-1][:-4]
    
    print(run)
    
    # Open .log file
    f = open(file,'r')
    line = f.readlines()
    
    # Get hypocenter coords
    hyp_loc_junk = line[16]
    hyplon = float(hyp_loc_junk.split(' ')[2].split('(')[1].split(')')[0].split(',')[0])
    hyplat = float(hyp_loc_junk.split(' ')[2].split('(')[1].split(')')[0].split(',')[1])
    hypdepth = float(hyp_loc_junk.split(' ')[2].split('(')[1].split(')')[0].split(',')[2])
    
    # Get hypocenter time
    hyp_time_junk = line[17]
    hyp_time = hyp_time_junk.split(' ')[2].split('Z')[0]
    
    # Epicentral time of simulations
    time_epi = UTCDateTime('2023-03-23T00:00:00.00000Z')
    
    # Loop over station types
    for stn_type in station_types:
        
        # Read in station file
        stn_df = pd.read_csv(f'/home/tnye/onc/data/station_info/{stn_type}_list.txt')
        
        # Get stations and coords
        stations = stn_df['#station']
        lon = stn_df['lon']
        lat = stn_df['lat']
        
        predictedP = 9999*ones(len(stations))
        predictedS = 9999*ones(len(stations))
    
        # Get predicted arrivals
        for ksta in range(len(stations)):
            
            # Find station coordinates
            stlon = lon[ksta]
            stlat = lat[ksta]
            stelev = metadata['Elevation'].iloc[np.where(np.array(metadata['Station'])==stations[ksta])[0][0]]
            
            # Compute Rhyp and Repi
            rhyp = emf.compute_rhyp(stlon,stlat,stelev,hyplon,hyplat,hypdepth)
            repi = emf.compute_repi(stlon,stlat,hyplon,hyplat)
            
            # Calculate P- and S-wave arrivals               
            deg = locations2degrees(stlat, stlon, hyplat, hyplon) # calculates distance between station loc and epicenter loc
            
            Ppaths=velmod.get_ray_paths(hypdepth, deg, phase_list=['P','p'])
            directP=Ppaths[0]
            predictedP[ksta] = directP.path['time'][-1]
            
            Spaths=velmod.get_ray_paths(hypdepth, deg, phase_list=['S','s'])
            directS=Spaths[0]
            predictedS[ksta] = directS.path['time'][-1]
            
            # Save arrival times to file   
            line = f'{batch},{run},{stations[ksta]},{rhyp},{repi},{time_epi + predictedP[ksta]},{time_epi + predictedS[ksta]}\n'
            w.write(line)

# Close file
w.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
