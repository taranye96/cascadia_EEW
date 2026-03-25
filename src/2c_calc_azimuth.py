#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 17:48:20 2024

@author: tnye
"""

##############################################################################
# This script calculates the azimuth and P-wave amplification for each
# simulated recording.
##############################################################################

# Imports
import numpy as np
import pandas as pd
from pyproj import Geod 
from numpy import cos,sin,deg2rad,rad2deg,arctan2

# Local imports
import cascadiaEEW_main_fns as emf

# STA/LTA dataframe
df = pd.read_csv('/Users/tnye/ONC/event_detection/sta_lta/cascadia_sta_lta_polar.csv')

# Station dataframe
stn_df = pd.read_csv('/Users/tnye/ONC/data/station_info/all_stations_metadata.csv')

# Define projection
p = Geod(ellps='WGS84')

# Path to velocity model
model = '/Users/tnye/ONC/simulations/cascadia/structure/cascadia.npz'

# Get list of runs and stations
runs = df['Run'].values
stns = df['Station'].values

# Initialize lists
az_list = []
baz_list = []
P_list = []

# Loop over rows of dataframe
for i in range(len(runs)):
    
    # Get run and station
    run = runs[i]
    stn = stns[i]

    # Read in the .log file and get event coordinates
    log_file = f'/Users/tnye/ONC/simulations/cascadia/output/ruptures/{run}.log'
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
    
    # Get station coordinates
    stn_ind = np.where(stn_df['Station']==stn)[0][0]
    stlon = stn_df['Longitude'][stn_ind]
    stlat = stn_df['Latitude'][stn_ind]

    # Calculate azimuth/back azimuth and append to lists
    az,baz,dist = p.inv(stlon,stlat,hyplon,hyplat)    
    az_list.append(az)
    baz_list.append(baz)
    
    # Read in .rupt file and get slip information
    rupt_file = np.genfromtxt(log_file.replace('.log','.rupt'))
    hyp_ind = np.where((rupt_file[:,1]==hyplon)&(rupt_file[:,2]==hyplat))[0][0]
    strike = rupt_file[hyp_ind,4]
    dip = rupt_file[hyp_ind,5]
    ss=rupt_file[hyp_ind,8]
    ds=rupt_file[hyp_ind,9]
    rake=rad2deg(arctan2(ds,ss))
    
    # Convert to radians
    rake=deg2rad(rake)                            
    strike=deg2rad(strike)
    dip=deg2rad(dip)
    azimuth=deg2rad(az)
    
    # Calculate takeoff angle
    take_off_angle=deg2rad(emf.get_takeoff_angle(hyp_ind,dist,rupt_file,model))    
    
    # Get factors for P-wave amplification                                          
    SR=sin(rake)                                                               
    CR=cos(rake)                                                               
    SD=sin(dip)                                                               
    CD=cos(dip)                                                               
    ST=sin(take_off_angle)                                                                
    CT=cos(take_off_angle)                                                                
    SS=sin(azimuth-strike)                                                            
    CS=cos(azimuth-strike)      
      
    # Calculate P-wave amplification factor                                                          
    P = CR*SD*ST**2*2*SS*CS - CR*CD*2*ST*CT*CS + SR*2*SD*CD*(CT**2-ST**2*SS**2) + SR*(CD**2-SD**2)*2*ST*CT*SS                                          
    P_list.append(P)

# Save to dataframe
data = {'Run':runs,'Station':stns,'Azimuth':az_list,'BackAzimuth':baz_list,'P':P_list}
df_out = pd.DataFrame(data)
df_out.to_csv('/Users/tnye/ONC/event_detection/azimuth/cascadia_azimuth-amplification.csv')

