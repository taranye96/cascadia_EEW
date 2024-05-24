#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 16:38:56 2023

@author: tnye
"""

# Imports
import numpy as np
import pandas as pd
from obspy import UTCDateTime
import matplotlib.pyplot as plt

# List of events
events = ['cascadia.000015','cascadia.000018','cascadia.000051','cascadia.000065',
          'cascadia.000096','cascadia.000108']

ONC_OT_df = pd.read_csv('/Users/tnye/ONC/EEW_analsysis/ONC/SimulationResults_For_Tara_23April2024/simulation_OTs.csv')

magnitudes = ['7.1','7.3','8.0','8.3','9.2','9.4']

lines = ['--',':']

fig, axs = plt.subplots(3,2,figsize=(6.5,6.5))
axes = [axs[0,0],axs[0,1],axs[1,0],axs[1,1],axs[2,0],axs[2,1]]

for i, event in enumerate(events):
   
    # ONC results
    onc_OT = ONC_OT_df['OriginTime'][np.where(ONC_OT_df['run']==event)[0][0]]
    onc_df = pd.read_csv(f"/Users/tnye/ONC/EEW_analsysis/ONC/SimulationResults_For_Tara_23April2024/Cas.{event.split('.')[1]}_CorrelatorResults_Detections_v2.csv")
    onc_df = onc_df.set_index(onc_df.reportedTime).sort_index(ascending=True)
     
    t_from_OT = [] 
    for j in range(len(onc_df)):
        t_from_OT.append(UTCDateTime(onc_df['reportedTime'][j]) - UTCDateTime(onc_OT))  
    
    if i != 3:
        axes[i].plot(t_from_OT,onc_df['magnitude'],c='indigo',lw=1,label=r'M$_{ONC-EW}$')
    else:
        axes[i].scatter(t_from_OT,onc_df['magnitude'],s=10,c='indigo')
        
    # axes[i].text(0.02,0.75,f'M{magnitudes[i]}',transform=axes[i].transAxes,ha='left',va='bottom',fontsize=11)
    axes[i].set_title(f'{event}: M{magnitudes[i]}')
    axes[i].grid(alpha=0.25)
    
    try:
        
        # EPIC results
        epic_summary = pd.read_csv(f'/Users/tnye/ONC/EEW_analsysis/ShakeAlert/{event}_summary_table.txt',delimiter='\t') 
        epic_df = pd.read_csv(f'/Users/tnye/ONC/EEW_analsysis/ShakeAlert/results_logs/{event}/{event}_event_summary.txt',delimiter='\t') 
        epic_df = epic_df.set_index(epic_df['seconds since OT']).sort_index(ascending=True)
                
        axes[i].plot(epic_df['seconds since OT'],epic_df['mag'],c='teal',label=r'M$_{EPIC}$')
           
    except:
        continue
    
    axes[i].axhline(float(magnitudes[i]),ls='-',c='goldenrod',lw=1,label='True M')
   

axes[0].set_xlim(xmax=110)
axes[1].set_xlim(xmax=80)
axes[2].set_xlim(xmax=100)
axes[3].set_xlim(xmax=60)
axes[4].set_xlim(xmax=66)
axes[5].set_xlim(xmax=120)

fig.text(0.5, 0.1, 'Seconds from Origin Time', ha='center', va='bottom', fontsize=12)
fig.text(0.01, 0.5, 'EEW Magnitude Estimate', ha='left', va='center', rotation='vertical', fontsize=12)

handles, labels = axs[0,0].get_legend_handles_labels()
axs[2,0].legend(handles,labels,bbox_to_anchor=(1.005,-0.6),loc='upper center',facecolor='white',fontsize=10,framealpha=0.98,ncol=3)
plt.subplots_adjust(left=0.1,right=0.975,top=0.95,hspace=0.55,bottom=0.2)
plt.savefig('/Users/tnye/ONC/manuscript/Figures/EEW_mag_evolution.png',dpi=300)

