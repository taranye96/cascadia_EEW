#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 10:24:33 2022

@author: tnye
"""

# Imports
from glob import glob
import pandas as pd
from numpy import genfromtxt,where,ones,arange,polyfit,tile,zeros,linspace,mean, nanmean,array
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import seaborn as sns

# Set parameters

windows = [1, 2, 4, 6, 12, 20]

dist_lim = 1000
dist_lims = [200,1200,500]

var_folder = '100q_1r'

# Pd_df = pd.read_csv('/Users/tnye/ONC/magnitude_estimation/Pd_cascadia_longer_wfs_withnoise.csv')
# Pd_df = pd.read_csv('/Users/tnye/ONC/magnitude_estimation/Pd_files/Pd_cascadia_longer_wfs_noisy_BB.csv')
Pd_df = pd.read_csv('/Users/tnye/ONC/magnitude_estimation/Pd_files/noisy_BB_notrunc.csv')
# kalman_df = pd.read_csv(f'/Users/tnye/ONC/magnitude_estimation/Pd_cascadia_longer_wfs_kalman_withnoise_{var_folder}_rot.csv')
gold_dfs = sorted(glob(f'/Users/tnye/ONC/magnitude_estimation/pd-Mw_scaling/*.csv'))

DTAmps=genfromtxt('/Users/tnye/ONC/files/DT2019_Amplitudes.txt')
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

# Begin plotting
plt.rcParams["font.family"] = "Helvetica"
hfont = {'fontname':'Helvetica'}
plt.rcParams['pdf.fonttype'] = 42
# sns.set_style('dark')
sns.set_style('whitegrid')

fig, axs = plt.subplots(2,3,figsize=(6.5,4.5))

for i in range(len(windows)):
    
    TW = windows[i]
    
    # Choose figure axis
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
    for j, dist_lim in enumerate(dist_lims):
        
    
        # Get Pd and Mw from dataframes
        Pd_df_trim = Pd_df[Pd_df['Repi'] < dist_lim]
        avg_Pd_df = Pd_df_trim.groupby('Event').agg({'Magnitude': 'first','Repi': 'first',f'{TW}s_Pd_logcm': 'mean'}).reset_index()
        var_Pd_df = Pd_df_trim.groupby('Event').agg({'Magnitude': 'first','Repi': 'first',f'{TW}s_Pd_logcm': 'var'}).reset_index()
        Pd = np.array(avg_Pd_df[f'{TW}s_Pd_logcm'])
        mw = np.array(avg_Pd_df['Magnitude'])
        # kal_rhyp_ind = np.where(kalman_df['Rhyp']<=200)[0]
        # kalman_pd = np.array(kalman_df[f'{TW}s_Pd_logcm'])[kal_rhyp_ind]
        # kalman_mw = np.array(kalman_df['Magnitude'])[kal_rhyp_ind]
    
        gold_df = pd.read_csv(gold_dfs[i],header=None)
    
        # Plot mean simualted event Pd
        ax.scatter(mw,Pd,marker='+',c=colors[j],ec=colors[j],lw=1,s=15,alpha=0.7,label=r'$\mu_{E,sim}$ '+f'< {dist_lim}km')
        
        if dist_lim == 200:
            sim_var.append(var_Pd_df[f'{TW}s_Pd_logcm'])
            
    ax.set_xlim([4.3,9.6])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.plot(gold_df[0],gold_df[1],c='black',label='GM20 model')
    
    # ax.scatter(kalman_mw, kalman_pd, marker='.', c='goldenrod', s=15, alpha=0.3, label='Kalman-Filtered ONC')
    ax.tick_params(bottom=True,left=True,labelbottom=True,labelleft=True)
    ax.text(0.975, 0.025, f'TW = {TW}s',color='k',ha='right',fontsize=10, transform=ax.transAxes)
    # ax.text(0.025, 0.025, f"Obs={round(np.mean(obs_var[:,i]),2)}",color='k',ha='left',fontsize=10, transform=ax.transAxes)
    # ax.text(0.025, 0.075, f"Sim={round(np.mean(sim_var),2)}",color='k',ha='left',fontsize=10, transform=ax.transAxes)
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
plt.savefig(f'/Users/tnye/ONC/manuscript/figures/Pd-MW_dist-split.png',dpi=300)
