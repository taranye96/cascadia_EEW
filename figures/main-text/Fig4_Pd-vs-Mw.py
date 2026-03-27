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
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.cm as cm

home = '/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects'

# Time windows
windows = [4, 6, 12, 20]
satM = [6.91, 7.2, 7.4, 7.55]

# Distance limits to group Pd values by
# dist_lims = [200,500,1200]
dist_lims = [200]

# Read in Pd dataframe
Pd_df = pd.read_csv(f'{home}/ONC-EEW/magnitude_estimation/Pd_files/Pd_cascadia_noisy.csv')

# fc1_df = pd.read_csv('/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW/magnitude_estimation/Pd_files/matched-filter_test/Pd_fc1.0.csv')
# fc0pt1_df = pd.read_csv('/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW/magnitude_estimation/Pd_files/matched-filter_test/Pd_fc0.1.csv')
# Pd_df = pd.read_csv('/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW/magnitude_estimation/Pd_files/Pd_cascadia_fc0.3.csv')


# Gather files with Goldberg and Melgar (2020) scaling models
gold_dfs = sorted(glob(f'{home}/ONC-EEW/magnitude_estimation/GM20_pd-Mw_scaling/Pd*.csv'))

# Read in file with observed Pd values (from Trugman et al., 2019)
DTAmps = pd.read_csv(f'{home}/ONC-EEW/files/DT2019_files/DT2019_Amplitudes.txt', delimiter='\t')
obs_events = np.unique(DTAmps['EventID'])
obs_Pd = np.array([])
obs_var = np.array([])
obs_M = np.array([])
for event in obs_events:
    event_Pd = []
    event_var = []
    obs_ind = np.where(DTAmps['EventID']==event)[0]
    for i in range(len(windows)):
        W = windows[i]
        event_Pd.append(np.mean(DTAmps.loc[:,f'{W}s'][obs_ind]))
        event_var.append(np.var(DTAmps.loc[:,f'{W}s'][obs_ind]))
    obs_Pd = np.append(obs_Pd,event_Pd)
    obs_var = np.append(obs_var,event_var)
    obs_M = np.append(obs_M, DTAmps['Magnitude'][obs_ind[0]])
obs_Pd = obs_Pd.reshape(-1, len(event_Pd))
obs_var = obs_var.reshape(-1, len(event_var))

# Make figure
plt.rcParams["font.family"] = "Helvetica"
hfont = {'fontname':'Helvetica'}
plt.rcParams['pdf.fonttype'] = 42
sns.set_style('whitegrid')

fig, axs = plt.subplots(2,2,figsize=(3.5,4.25))

norm = mcolors.Normalize(vmin=0, vmax=1200)

# Loop over time windows
wi = 0
for i in range(2):
    for j in range(2):
    
        TW = windows[wi]
        
        # Plot all data points
        axs[i,j].axvline(satM[wi], c='k', lw=1, label='saturation magnitude\n(Trugman et al., 2019)', zorder=1)
        axs[i,j].scatter(DTAmps['Magnitude'].values,DTAmps[f'{TW}s'].values,marker='.',fc='lightsteelblue',ec='none',s=3, alpha=1,label='individual recordings')
        cvals = Pd_df['Repi'].to_numpy(dtype=float)
        axs[i,j].scatter(Pd_df['Magnitude'], Pd_df[f'{TW}s_Pd_logcm'], marker='.', c=cvals, cmap='magma_r', ec='none', vmin=0, vmax=1200, s=3, alpha=1, zorder=1)

        axs[i,j].scatter(obs_M,obs_Pd[:,wi], marker='d', fc='mediumblue', ec='k', lw=0.3, s=10, alpha=1, label=r'$\mu_{E,obs}$')
        
        # Loop over distance limits
        for k, dist_lim in enumerate(dist_lims):
        
            # Trim dataframe to only include data within distance limit
            if k != 0:
                Pd_df_trim = Pd_df[(Pd_df['Repi'] <= dist_lim) & (Pd_df['Repi'] > dist_lims[k-1])]
            else:
                Pd_df_trim = Pd_df[(Pd_df['Repi'] <= dist_lim)]
                                   
            # Organize dataframe to get mean and variance for each event
            avg_Pd_df = Pd_df_trim.groupby('Event').agg({'Magnitude': 'first','Repi': 'first',f'{TW}s_Pd_logcm': 'mean'}).reset_index()
            var_Pd_df = Pd_df_trim.groupby('Event').agg({'Magnitude': 'first','Repi': 'first',f'{TW}s_Pd_logcm': 'var'}).reset_index()
            Pd = np.array(avg_Pd_df[f'{TW}s_Pd_logcm'])
            mw = np.array(avg_Pd_df['Magnitude'])

            axs[i,j].scatter(mw, Pd, marker='d', fc=cm.get_cmap('magma_r')(norm(200)), ec='k', lw=0.4, s=12, alpha=1, label=r'$\mu_{E,sim}$ '+f'< {dist_lim}km', zorder=4)
        
        ## Read in Pd-M scaling model for the time window
        # gold_df = pd.read_csv(f'{home}/ONC-EEW/magnitude_estimation/GM20_pd-Mw_scaling/Pd-{TW}s.csv',header=None)
        # axs[i,j].plot(gold_df[0],gold_df[1],c='black',label='GM20 model',zorder=4)
        
        axs[i,j].xaxis.set_major_locator(MultipleLocator(1))
        axs[i,j].yaxis.set_major_locator(MultipleLocator(1))
        axs[i,j].tick_params(bottom=True,left=True,labelbottom=True,labelleft=True)
        axs[i,j].text(0.025, 0.975, f'TW = {TW}s',color='k',ha='left',va='top',fontsize=8, transform=axs[i,j].transAxes)
        axs[i,j].grid(alpha=0)
        axs[i,j].set_xlim(4.25,9.5)
        axs[i,j].set_ylim(-3.2,2.2)
        axs[i,j].tick_params(which='major', length=3, width=1)
        # print(round(np.mean(obs_var[:,i])/np.mean(sim_var),2))
        
        if j == 1:
            axs[i,j].set_yticklabels([])
        if i == 0:
            axs[i,j].set_xticklabels([])
            
        wi+=1

handles, labels = axs[0,0].get_legend_handles_labels()
opaque_handles = []
opaque_labels  = []

sat_handle = None
sat_label  = 'saturation magnitude\n(Trugman et al., 2019)'

for h, label in zip(handles, labels):
    
    if label == 'individual recordings':
        proxy = Line2D([0], [0],
                       marker='o',
                       linestyle='None',
                       markerfacecolor='darkgray',
                       markeredgecolor='darkgray',
                       markersize=1)
        out_label = label
        
    elif label.startswith(r'$\mu_{E'):
        marker = 'd'
        color = 'mediumblue' if 'obs' in label else (
            cm.get_cmap('magma_r')(norm(200))
        )
        proxy = Line2D([0], [0],
                       marker=marker,
                       linestyle='None',
                       markerfacecolor=color,
                       markeredgecolor='k',
                       markeredgewidth=0.6,
                       markersize=3)
        out_label = label
        
    elif label == 'GM20 model':
        proxy = Line2D([0], [0],
                       linestyle='-',
                       color='black')
        out_label = label
        
    elif 'saturation magnitude' in label:
        proxy = Line2D([0], [0],
                       linestyle='-',
                       color='k',
                       linewidth=1)
        sat_handle = proxy
        continue
        
    else:
        proxy = h
        out_label = label
        
    opaque_handles.append(proxy)
    opaque_labels.append(out_label)

# Append saturation entry last
if sat_handle is not None:
    opaque_handles.append(sat_handle)
    opaque_labels.append(sat_label)

# Add legend
axs[1,0].legend(opaque_handles,opaque_labels,bbox_to_anchor=(0.5,-0.375),loc='upper center',facecolor='white',fontsize=8,framealpha=0.98,ncol=1,markerscale=1.5)
fig.supxlabel('Magnitude', **hfont, fontsize=10, x=0.55, y=0.25)
fig.supylabel('log$\mathregular{_{10}}$(Pd)$\mathregular{_{10km}}$ (cm)', **hfont, fontsize=10, y=0.64)
plt.subplots_adjust(left=0.15,right=0.975,top=0.975,bottom=0.35,wspace=0.075,hspace=0.075)  

# Add colorbar
norm = mcolors.Normalize(vmin=0, vmax=1200)
cmap = plt.cm.magma_r
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cax = fig.add_axes([0.59, 0.15, 0.375, 0.025])  # left, bottom, width, height
cbar = fig.colorbar(sm, cax,
                    shrink=1,   
                    pad=0.15,
                    anchor=(0.0, 1),
                    orientation='horizontal')   
cbar.ax.xaxis.set_major_locator(MultipleLocator(250))
# cbar.ax.xaxis.set_minor_locator(MultipleLocator(250))
cbar.ax.tick_params(which='major', bottom=True, top=False, length=4, labelsize=8)
# cbar.ax.tick_params(which='minor', bottom=True, top=False, length=2)
cbar.set_label(r'$R_{epi}$ (km)', fontsize=8)


plt.savefig(f'{home}/ONC-EEW/manuscript/figures/main-text/Fig4_Pd-MW.png',dpi=500)


#%%

# Make figure
plt.rcParams["font.family"] = "Helvetica"
hfont = {'fontname':'Helvetica'}
plt.rcParams['pdf.fonttype'] = 42
sns.set_style('whitegrid')

fig, ax = plt.subplots(1,1,figsize=(3.5,3.75))
    
TW = 4
i = 2

# Plot all data points
ax.scatter(DTAmps[:,1],DTAmps[:,i+3],marker='.',c='darkgray',s=1, alpha=0.3,label='individual recordings')
ax.scatter(Pd_df['Magnitude'], Pd_df[f'{TW}s_Pd_logcm'], marker='.', c='darkgray', s=5, alpha=0.3)
ax.scatter(obs_M,obs_Pd[:,i],marker='d',fc='none',ec='indigo',lw=0.4,s=25,alpha=0.7, label=r'$\mu_{E,obs}$ (Trugman et al., 2019)')

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
    ax.scatter(mw, Pd, marker='o', fc='none', ec=colors[j], lw=0.7, s=25, alpha=1, label=r'$\mu_{E,sim}$ '+f'< {dist_lim}km')
    
    if dist_lim == 200:
        sim_var.append(var_Pd_df[f'{TW}s_Pd_logcm'])
        
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.plot(gold_df[0],gold_df[1],c='black',label='GM20 model')
ax.tick_params(bottom=True,left=True,labelbottom=True,labelleft=True)
ax.text(0.975, 0.025, f'TW = {TW}s',color='k',ha='right',fontsize=10, transform=ax.transAxes)
ax.grid(alpha=0)
ax.set_xlim(4.5,9.5)
ax.set_ylim(-4,2)
ax.set_xlabel('Magnitude')
print(round(np.mean(obs_var[:,i])/np.mean(sim_var),2))

handles, labels = axs[1,1].get_legend_handles_labels()
opaque_handles = []
for h, label in zip(handles, labels):
    if label == 'individual recordings':
        # Filled gray dot
        proxy = Line2D([0], [0],
                       marker='o',
                       linestyle='None',
                       markerfacecolor='darkgray',
                       markeredgecolor='darkgray',
                       markersize=3,
                       label=label)
    elif label.startswith(r'$\mu_{E'):  # For open diamond Trugman or open circle sim
        marker = 'd' if 'obs' in label else 'o'
        edgecolor = 'indigo' if 'obs' in label else (
            'mediumblue' if '< 200' in label else
            'green' if '< 1200' in label else
            'goldenrod'
        )
        proxy = Line2D([0], [0],
                       marker=marker,
                       linestyle='None',
                       markerfacecolor='none',
                       markeredgecolor=edgecolor,
                       markeredgewidth=0.7,
                       markersize=5,
                       label=label)
    elif label == 'GM20 model':
        proxy = Line2D([0], [0],
                       linestyle='-',
                       color='black',
                       label=label)
    else:
        proxy = h  # fallback
    opaque_handles.append(proxy)
    
labels = np.array(labels)[np.array([0,1,5,2,4,3])]
opaque_handles = np.array(opaque_handles)[np.array([0,1,5,2,4,3])]
    
ax.legend(opaque_handles,labels,bbox_to_anchor=(0.475,-0.2),loc='upper center',facecolor='white',fontsize=8,framealpha=0.98,ncol=2,markerscale=1.25)
ax.set_xlabel('Magnitude',**hfont,fontsize=10)
ax.set_ylabel('log$\mathregular{_{10}}$(Pd)$\mathregular{_{10km}}$ (cm)', **hfont, fontsize=10)

plt.subplots_adjust(left=0.15,right=0.95,top=0.96,bottom=0.325,wspace=0.1,hspace=0.15)  

plt.savefig(f'{home}/ONC-EEW/manuscript/figures/main-text/Fig4_Pd-MW_dist-split_4s.png',dpi=500)


#%%

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
    else:
        ax = axs[1,i-3]
    
    # Plot all data points
    DTAmps[:,i+3] = DTAmps[:,i+3] - np.mean(DTAmps[:,i+3])
    ax.hist(DTAmps[:,i+3],bins=np.arange(-2.5,2.6,0.25), label='observed')
    ax.set_title(f'{TW} s')
    # ax.scatter(Pd_df['Magnitude'], Pd_df[f'{TW}s_Pd_logcm'], marker='.', c='darkgray', s=5, alpha=0.3)
    # ax.scatter(obs_M,obs_Pd[:,i],marker='d',fc='none',ec='indigo',lw=0.4,s=20,alpha=0.7, label=r'$\mu_{E,obs}$ (Trugman et al., 2019)')
    
    # Loop over distance limits
    dist_lim = 1200
    
    # Trim dataframe to only include data within distance limit
    Pd_df_trim = Pd_df[Pd_df['Repi'] < dist_lim]
    Pd_df_trim.loc[:,f'{TW}s_Pd_logcm'] = Pd_df_trim.loc[:,f'{TW}s_Pd_logcm'] - np.mean(Pd_df_trim.loc[:,f'{TW}s_Pd_logcm'])
    
    ax.hist(Pd_df_trim.loc[:,f'{TW}s_Pd_logcm'],bins=np.arange(-2.5,2.6,0.25), label='simulated')
    
    # Organize dataframe to get mean and variance for each event
    avg_Pd_df = Pd_df_trim.groupby('Event').agg({'Magnitude': 'first','Repi': 'first',f'{TW}s_Pd_logcm': 'mean'}).reset_index()
    var_Pd_df = Pd_df_trim.groupby('Event').agg({'Magnitude': 'first','Repi': 'first',f'{TW}s_Pd_logcm': 'var'}).reset_index()
    Pd = np.array(avg_Pd_df[f'{TW}s_Pd_logcm'])

    sim_var = var_Pd_df[f'{TW}s_Pd_logcm']
    
    
    ax.set_xlim(-2.5,2.5)
        
    print(round(np.mean(obs_var[:,i])/np.mean(sim_var),2))

axs[1,1].set_xlabel('Demeaned log$\mathregular{_{10}}$(Pd)$\mathregular{_{10km}}$ (cm)')
axs[0,0].set_ylabel('Counts')
axs[1,0].set_ylabel('Counts')
axs[1,1].legend(loc='upper center', bbox_to_anchor=(0.5,-0.35), ncol=2)

plt.subplots_adjust(left=0.115,right=0.975,top=0.925,bottom=0.2,wspace=0.5,hspace=0.4)  
plt.savefig(f'{home}/ONC-EEW/manuscript/figures/response-to-reviews/Pd_hist.png', dpi=300)







