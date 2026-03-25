#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 16:38:56 2023

@author: tnye
"""

###############################################################################
# This script makes Figure 10 in the paper, which shows the simulation
# magnitude estimates from the EEW systems. 
###############################################################################

# Imports
import numpy as np
import pandas as pd
from obspy import UTCDateTime
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D

home = '/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW'

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
    
    # Plot EPIC results     
    axes[i].plot(epic_df['seconds since OT'],epic_df['mag'],c='green',lw=0.5)
    axes[i].scatter(epic_df['seconds since OT'],epic_df['mag'],ec='green',fc='none',marker='d',s=10,lw=0.5,label='EPIC')
    
    axes[i].yaxis.set_major_locator(MultipleLocator(1))
    axes[i].yaxis.set_minor_locator(MultipleLocator(0.5))
   
    # Read in ONC results and sort by time
    onc_OT = ONC_OT_df['OriginTime'][np.where(ONC_OT_df['run']==event)[0][0]]
    onc_df = pd.read_csv(f"{home}/EEW_analsysis/ONC/SimulationResults_For_Tara_23April2024/Cas.{event.split('.')[1]}_CorrelatorResults_Detections_v2.csv")
    onc_df = onc_df.set_index(onc_df.reportedTime).sort_index(ascending=True)
    
    # Get time from origin 
    t_from_OT = [] 
    for j in range(len(onc_df)):
        t_from_OT.append(UTCDateTime(onc_df['reportedTime'][j]) - UTCDateTime(onc_OT))  
    
    # Since ONC-EW only had two updates for  event 65, plot that as scatter
    if i != 3:
        axes[i].plot(t_from_OT,onc_df['magnitude'],c='blue',lw=0.5)
    else:
        axes[i].scatter(t_from_OT,onc_df['magnitude'],s=10,c='blue')
    
    axes[i].scatter(t_from_OT,onc_df['magnitude'],ec='blue',fc='none',marker='^',s=15,lw=0.5,label='ONC-EW')
        
    axes[i].set_title(f'Scenario {evnum}', size=9)
    axes[i].grid(alpha=0)
    axes[i].text(0.05, 0.95, f'M {magnitudes[i]}', transform=axes[i].transAxes, ha='left', va='top', fontsize=9)
    
    # Add line for the true magnitude
    axes[i].axhline(float(magnitudes[i]),ls='-',c='goldenrod',lw=1.2,label='True M',zorder=1)
    
    # Axes
    if i in [1,3,5]:
        axes[i].set_yticklabels([])
    else:
        axes[i].set_ylabel('EEW M Estimate', labelpad=10, fontsize=9)
    
    if i in [4,5]:
        axes[i].set_xlabel('Time from origin (s)', fontsize=9)

# Adjust maximum xlimits
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

axes[0].set_ylim(5.5,9.5)
axes[1].set_ylim(5.5,9.5)
axes[2].set_ylim(5.5,9.5)
axes[3].set_ylim(5.5,9.5)
axes[4].set_ylim(5.5,9.5)
axes[5].set_ylim(5.5,9.5)

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

trueM_handle = Line2D(
    [],[],
    color='goldenrod',
    lw=1,
    label='True M')

axs[2,0].legend(
    handles=[trueM_handle, epic_handle, onc_handle],
    bbox_to_anchor=(1.055, -0.55),
    loc='upper center',
    facecolor='white',
    fontsize=8,
    ncol=3
)

# Adjust suplot spacing
plt.subplots_adjust(left=0.17,right=0.95,top=0.95,hspace=0.55,wspace=0.1,bottom=0.175)

# Save figure
plt.savefig(f'{home}/manuscript/Figures/main-text/Fig9_EEW_mag_evolution.png',dpi=300)

