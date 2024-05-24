#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 13:47:19 2023

@author: tnye
"""

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Geod
from obspy import UTCDateTime

p = Geod(ellps='WGS84')


# List of events
events = ['cascadia.000015','cascadia.000018','cascadia.000051','cascadia.000065',
          'cascadia.000096','cascadia.000108']

ONC_OT_df = pd.read_csv('/Users/tnye/ONC/EEW_analsysis/ONC/SimulationResults_For_Tara_23April2024/simulation_OTs.csv')


magnitudes = ['7.1','7.3','8.0','8.3','9.2','9.4']


lines = ['--',':']

fig, axs = plt.subplots(3,2,figsize=(6.5,6.5))
axes = [axs[0,0],axs[0,1],axs[1,0],axs[1,1],axs[2,0],axs[2,1]]

for i, event in enumerate(events):
        
    epic_df = pd.read_csv(f'/Users/tnye/ONC/EEW_analsysis/ShakeAlert/results_logs/{event}/{event}_event_summary.txt',delimiter='\t') 
    epic_df = epic_df.set_index(epic_df['seconds since OT']).sort_index(ascending=True)
   
    onc_OT = ONC_OT_df['OriginTime'][np.where(ONC_OT_df['run']==event)[0][0]]
    onc_df = pd.read_csv(f"/Users/tnye/ONC/EEW_analsysis/ONC/SimulationResults_For_Tara_23April2024/Cas.{event.split('.')[1]}_CorrelatorResults_Detections_v2.csv")
    onc_df = onc_df.set_index(onc_df.reportedTime).sort_index(ascending=True)
    
    t_from_OT = [] 
    for j in range(len(onc_df)):
        t_from_OT.append(UTCDateTime(onc_df['reportedTime'][j]) - UTCDateTime(onc_OT))  
        
    # Get epicenter location error
    hyplon = epic_df['event lon'].iloc[0]
    hyplat = epic_df['event lat'].iloc[0]
    _,_,onc_epi_error = p.inv([hyplon]*len(onc_df),[hyplat]*len(onc_df),onc_df['longitude'].values,onc_df['latitude'].values)
    onc_epi_error= onc_epi_error/1000
    
    if i != 3:
        axes[i].plot(t_from_OT,onc_epi_error,c='indigo',lw=1,label=r'$ONC-EW$')
    else:
        axes[i].scatter(t_from_OT,onc_epi_error,s=10,c='indigo')
    axes[i].plot(epic_df['seconds since OT'],epic_df['loc err km'],c='teal',lw=1,label=r'$EPIC$')

    axes[i].axhline(0,ls='-',c='k',lw=1)
    # axes[i].text(0.02,0.8,f'M{magnitudes[i]}',transform=axes[i].transAxes,ha='left',va='bottom',fontsize=11)
    axes[i].set_title(f'{event}: M{magnitudes[i]}')
    axes[i].grid(alpha=0.25)
    
    # print(np.min(onc_epi_error))
    # print(onc_epi_error[-1])
    
    axes[i].set_ylim(ymin=-0.5)

axes[0].set_xlim(xmax=80)
axes[1].set_xlim(xmax=80)
axes[2].set_xlim(xmax=100)
axes[3].set_xlim(xmax=60)
axes[4].set_xlim(xmax=66)
fig.text(0.5, 0.1, 'Seconds from Origin Time', ha='center', va='bottom', fontsize=12)
fig.text(0.01, 0.5, 'Epicenter Location Error (km)', ha='left', va='center', rotation='vertical', fontsize=12)

handles, labels = axs[0,0].get_legend_handles_labels()
axs[2,0].legend(handles,labels,bbox_to_anchor=(1.005,-0.6),loc='upper center',facecolor='white',fontsize=10,framealpha=0.98,ncol=3)
plt.subplots_adjust(left=0.1,right=0.975,top=0.95,hspace=0.55,bottom=0.2)
plt.savefig('/Users/tnye/ONC/manuscript/Figures/EEW_loc-error_evolution.png',dpi=300)

