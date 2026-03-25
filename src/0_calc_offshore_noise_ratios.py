#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 09:15:29 2022

@author: tnye
"""

###############################################################################
# This script calculates the ratios of noise between one of the land sites
# and all the offshore sites based on a day of noise data. 
###############################################################################

# Imports
import numpy as np
from obspy import read
import matplotlib.pyplot as plt

stations = ['BACME','CBC27','NC89','CQS64']

st_land = read('/Users/tnye/ONC/data/station_noise/broadband/corrected_noise/NV_AL2H_2022-02-20.mseed')

db = []
for stn in stations:
    
    st_offshore = read(f'/Users/tnye/ONC/data/station_noise/offshore/uncut/NV_{stn}_2022-02-21.mseed')
    
    for tr_land in st_land:
        if tr_land.stats.channel == 'HNZ':
            break
    for tr_offshore in st_offshore:
        if tr_offshore.stats.channel == 'HNZ':
            break
    
    sind = np.where(tr_land.times()==tr_offshore.times()[0])[0][0]
    eind = np.where(tr_land.times()==tr_offshore.times()[-1])[0][0]
    land_data = tr_land.data[sind:eind+1]
    offshore_data = tr_offshore.data[::2]
    
    time = tr_land.times()[sind:eind+1]

    # Calculate difference in decibels
    land_avg = np.mean(np.abs(land_data))
    offshore_avg = np.mean(np.abs(offshore_data))
    
    db.append(20*np.log10(offshore_avg/land_avg))
    
    # Plot noise
    fig, ax = plt.subplots(1,1,figsize=(6,4))
    ax.plot(time, land_data, label='Land: AL2H')
    ax.plot(time, offshore_data, label=f'Offshore: {stn}')
    ax.set_xlabel('Time')
    ax.set_ylabel(r'Acceleration ($m/s^2$)')
    plt.legend(loc='upper right')
    plt.subplots_adjust(left=0.17, right=0.95, top=0.95, bottom=0.12)
    plt.savefig(f'/Users/tnye/ONC/plots/offshore_noise/{stn}.png',dpi=300)
    plt.close()

outfile = open('/Users/tnye/ONC/outfiles/offshore_noise_2022-02-20.txt', 'w')
file_data = np.array([stations, db],dtype=object)
file_data = file_data.T
outfile.write(f'#Station \t dB from AL2H noise\n')
np.savetxt(outfile, file_data, fmt=['%s', '%E'], delimiter='\t')
outfile.close()


