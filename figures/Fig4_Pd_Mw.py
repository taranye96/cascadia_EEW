#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 10:24:33 2022

@author: tnye
"""

##############################################################################
# This script makes Figure 4 in the paper, which shows the P-wave amplitude
# (Pd) – magnitude scaling for different time windows. 
##############################################################################

# Imports
from glob import glob
import pandas as pd
from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns

# Time windows
windows = [1, 2, 4, 6, 12, 20]

# Distance limits to group Pd values by
dist_lims = [200,1200,500]

# Read in Pd dataframe
Pd_df = pd.read_csv('/Users/tnye/ONC/magnitude_estimation/Pd_files/Pd_cascadia_noisy.csv')

# Gather files with Goldberg and Melgar (2020) scaling models
gold_dfs = sorted(glob('/Users/tnye/ONC/magnitude_estimation/GM20_pd-Mw_scaling/*.csv'))

# Read in file with observed Pd values (from Trugman et al., 2019)
DTAmps=genfromtxt('/Users/tnye/ONC/files/DT2019_files/DT2019_Amplitudes.txt')
obs_events = np.unique(DTAmps[:,0])
obs_Pd = np.array([])
obs_var = np.array([])
obs_M = np.array([])
for event in obs_events:
    event_Pd = []
    event_var = []
    obs_ind = np.where(DTAmps[:,0]==event)[0]
    for i in range(len(windows)):
        event_Pd.append(np.mean(DTAmps[:,i+3][obs_ind]))
        event_var.append(np.var(DTAmps[:,i+3][obs_ind]))
    obs_Pd = np.append(obs_Pd,event_Pd)
    obs_var = np.append(obs_var,event_var)
    obs_M = np.append(obs_M, DTAmps[:,1][obs_ind[0]])
obs_Pd = obs_Pd.reshape(-1, len(event_Pd))
obs_var = obs_var.reshape(-1, len(event_var))

# Make figure
plt.rcParams["font.family"] = "Helvetica"
hfont = {'fontname':'Helvetica'}
plt.rcParams['pdf.fonttype'] = 42
sns.set_style('whitegrid')

fig, axs = plt.subplots(2,3,figsize=(6.5,4.5))

# Loop over time windows
for i in range(len(windows)):
    
    TW = windows[i]
    
    if i < 3:
        ax = axs[0,i]
        ax.set_xticklabels([])
    else:
        ax = axs[1,i-3]
    
    # Plot all data points
    ax.scatter(DTAmps[:,1],DTAmps[:,i+3],marker='.',c='darkgray',s=1, alpha=0.3,label='individual recordings')
    ax.scatter(Pd_df['Magnitude'],Pd_df[f'{TW}s_Pd_logcm'],marker='.',c='darkgray',s=1, alpha=0.3)
    ax.scatter(obs_M,obs_Pd[:,i],marker='d',c='indigo',ec='indigo',lw=0.1,s=15,alpha=0.7, label=r'$\mu_{E,obs}$ (Trugman et al., 2019)')
    
    colors = ['mediumblue','green','goldenrod']
    sim_var = []
    
    # Loop over distance limits
    for j, dist_lim in enumerate(dist_lims):
    
        # Trim dataframe to only include data within distance limit
        Pd_df_trim = Pd_df[Pd_df['Repi'] < dist_lim]
        
        # Organize dataframe to get mean and variance for each event
        avg_Pd_df = Pd_df_trim.groupby('Event').agg({'Magnitude': 'first','Repi': 'first',f'{TW}s_Pd_logcm': 'mean'}).reset_index()
        var_Pd_df = Pd_df_trim.groupby('Event').agg({'Magnitude': 'first','Repi': 'first',f'{TW}s_Pd_logcm': 'var'}).reset_index()
        Pd = np.array(avg_Pd_df[f'{TW}s_Pd_logcm'])
        mw = np.array(avg_Pd_df['Magnitude'])
    
        # Read in Pd-M scaling model for the time window
        gold_df = pd.read_csv(gold_dfs[i],header=None)
    
        # Plot mean simulated event Pd
        ax.scatter(mw,Pd,marker='+',c=colors[j],ec=colors[j],lw=1,s=15,alpha=0.7,label=r'$\mu_{E,sim}$ '+f'< {dist_lim}km')
        
        if dist_lim == 200:
            sim_var.append(var_Pd_df[f'{TW}s_Pd_logcm'])
            
    ax.set_xlim([4.3,9.6])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.plot(gold_df[0],gold_df[1],c='black',label='GM20 model')
    ax.tick_params(bottom=True,left=True,labelbottom=True,labelleft=True)
    ax.text(0.975, 0.025, f'TW = {TW}s',color='k',ha='right',fontsize=10, transform=ax.transAxes)
    ax.grid(alpha=0.5)
    ax.set_xlim(4.5,9.5)
    print(round(np.mean(obs_var[:,i])/np.mean(sim_var),2))

handles, labels = axs[1,1].get_legend_handles_labels()
for lh in handles: 
    lh.set_alpha(1)
axs[1,1].legend(handles,labels,bbox_to_anchor=(0.5,-0.25),loc='upper center',facecolor='white',fontsize=10,framealpha=0.98,ncol=3,markerscale=1.5)
axs[1,1].set_xlabel('Magnitude',**hfont,fontsize=10)
fig.supylabel('log$\mathregular{_{10}}$(Pd)$\mathregular{_{10km}}$ (cm)',**hfont,fontsize=10)
plt.subplots_adjust(left=0.1,right=0.975,top=0.975,bottom=0.23,wspace=0.25,hspace=0.15)  
plt.savefig('/Users/tnye/ONC/manuscript/figures/Fig4_Pd-MW_dist-split.png',dpi=300)
