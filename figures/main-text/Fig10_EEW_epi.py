#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 13:47:19 2023

@author: tnye
"""

###############################################################################
# This script makes Figure 11 in the paper, which shows the simulation
# epicenter location error from the EEW systems. 
###############################################################################

home = '/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW'


# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Geod
from obspy import UTCDateTime
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D

# Define the projection
p = Geod(ellps='WGS84')

# List of events for EEW performance test
events = ['cascadia.000015','cascadia.000018','cascadia.000051','cascadia.000065',
          'cascadia.000096','cascadia.000108']

# Read in dataframe with origintimes associated with the ONC performance tests
ONC_OT_df = pd.read_csv(f'{home}/EEW_analsysis/ONC/SimulationResults_For_Tara_23April2024/simulation_OTs.csv')

# Magnitudes of the test events
magnitudes = ['7.1','7.3','8.0','8.3','9.2','9.4']

fig, axs = plt.subplots(3,2,figsize=(3.5,5.5))
axes = [axs[0,0],axs[0,1],axs[1,0],axs[1,1],axs[2,0],axs[2,1]]

# Loop over the test events
for i, event in enumerate(events):
    
    evnum = int(event.split('.')[1])
    
    # Read in EPIC results dataframe and sort by time    
    epic_df = pd.read_csv(f'{home}/EEW_analsysis/ShakeAlert/results_logs/{event}/{event}_event_summary.txt',delimiter='\t') 
    epic_df = epic_df.set_index(epic_df['seconds since OT']).sort_index(ascending=True)
   
    # Read in ONC results and sort by time
    onc_OT = ONC_OT_df['OriginTime'][np.where(ONC_OT_df['run']==event)[0][0]]
    onc_df = pd.read_csv(f"{home}/EEW_analsysis/ONC/SimulationResults_For_Tara_23April2024/Cas.{event.split('.')[1]}_CorrelatorResults_Detections_v2.csv")
    onc_df = onc_df.set_index(onc_df.reportedTime).sort_index(ascending=True)
    
    # Get time from origin 
    t_from_OT = [] 
    for j in range(len(onc_df)):
        t_from_OT.append(UTCDateTime(onc_df['reportedTime'][j]) - UTCDateTime(onc_OT))  
        
    # Get epicenter location error
    hyplon = epic_df['event lon'].iloc[0]
    hyplat = epic_df['event lat'].iloc[0]
    _,_,onc_epi_error = p.inv([hyplon]*len(onc_df),[hyplat]*len(onc_df),onc_df['longitude'].values,onc_df['latitude'].values)
    onc_epi_error= onc_epi_error/1000
    
    axes[i].plot(epic_df['seconds since OT'],epic_df['loc err km'],c='green',lw=0.5)
    axes[i].scatter(epic_df['seconds since OT'],epic_df['loc err km'],ec='green',fc='none',marker='d',s=10,lw=0.5,label='EPIC')
    
    # Since ONC-EW only had two updates for  event 65, plot that as scatter
    axes[i].plot(t_from_OT,onc_epi_error,c='blue',lw=0.5)
    axes[i].scatter(t_from_OT,onc_epi_error,ec='blue',fc='none',marker='^',s=15,lw=0.5,label='ONC-EW')
     
    # Miscellaneous formatting
    axes[i].set_title(f'Scenario {evnum}', size=9)
    axes[i].grid(alpha=0)
    axes[i].text(0.05, 0.95, f'M {magnitudes[i]}', transform=axes[i].transAxes, ha='left', va='top', fontsize=9)
    
    axes[i].yaxis.set_major_locator(MultipleLocator(10))
    axes[i].yaxis.set_minor_locator(MultipleLocator(5))
    
    axes[i].set_ylim(-1,20)
    
    # print(event)
    # print(f"EPIC: {np.max(epic_df['loc err km'])}")
    # print(f"ONC: {np.max(onc_epi_error)}")
    # print('')
    
    if i in [1,3,5]:
        axes[i].set_yticklabels([])
    else:
        axes[i].set_ylabel(r'EEW $\delta_{epi}$ (km)', fontsize=9)
    
    if i in [4,5]:
        axes[i].set_xlabel('Time from origin (s)', fontsize=9)

axes[0].set_xlim(0,110)
axes[0].xaxis.set_major_locator(MultipleLocator(50))
axes[0].xaxis.set_minor_locator(MultipleLocator(25))

axes[1].set_xlim(0,80)
axes[1].xaxis.set_major_locator(MultipleLocator(40))
axes[1].xaxis.set_minor_locator(MultipleLocator(20))

axes[2].set_xlim(0,110)
axes[2].xaxis.set_major_locator(MultipleLocator(50))
axes[2].xaxis.set_minor_locator(MultipleLocator(25))

axes[3].set_xlim(0,60)
axes[3].xaxis.set_major_locator(MultipleLocator(30))
axes[3].xaxis.set_minor_locator(MultipleLocator(15))

axes[4].set_xlim(0,66)
axes[4].xaxis.set_major_locator(MultipleLocator(30))
axes[4].xaxis.set_minor_locator(MultipleLocator(15))

axes[5].set_xlim(0,110)
axes[5].xaxis.set_major_locator(MultipleLocator(50))
axes[5].xaxis.set_minor_locator(MultipleLocator(25))

# Set legend
epic_handle = Line2D(
    [], [],
    color='green',
    lw=0.5,
    marker='d',
    markerfacecolor='none',
    markeredgecolor='green',
    markersize=6,
    label='EPIC'
)

onc_handle = Line2D(
    [], [],
    color='blue',
    lw=0.5,
    marker='^',
    markerfacecolor='none',
    markeredgecolor='blue',
    markersize=6.5,
    label='ONC-EW'
)

axs[2,0].legend(
    handles=[epic_handle, onc_handle],
    bbox_to_anchor=(1.07, -0.55),
    loc='upper center',
    facecolor='white',
    fontsize=8,
    ncol=2
)

# Adjust suplot spacing
plt.subplots_adjust(left=0.17,right=0.95,top=0.95,hspace=0.55,wspace=0.1,bottom=0.175)

# Save figure
plt.savefig(f'{home}/manuscript/Figures/main-text/Fig10_EEW_loc-error_evolution.png',dpi=300)

