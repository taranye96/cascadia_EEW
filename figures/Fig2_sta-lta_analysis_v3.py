#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 16:53:28 2023

@author: tnye
"""

##############################################################################
# This script makes Figure 3 in the paper, which shows the results of the 
# STA/LTA event detection analysis. 
##############################################################################

# Imports
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.lines as mlines

home = '/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW'

# Name of batch for simulations
batch = 'cascadia'

# Are we wanting the  'obspy' or 'elarmS' picker for the STA/LTA results using
    # EPIC/ShakeAlert parameters?
shakealert_dtype = 'obspy'

# Get azimuth and P-wave amplification factor
az_df = pd.read_csv(f'{home}/event_detection/azimuth/cascadia_azimuth-amplification.csv')
az = az_df['Azimuth']
P = az_df['P']

# Read in STA/LTA file
pick_file = f'{home}/event_detection/sta_lta/{batch}_sta_lta_polar.csv'
df = pd.read_csv(pick_file)

# Get STA/LTA triggers (Yes or No) 
trig_OncEW = df['1/5/4.9_triggered']
trig_ShakeAlert = df[f'0.05/5/20_{shakealert_dtype}_triggered']

# Get successful and unsuccessful trigger indices
OncEW_nan_ind = np.where(trig_OncEW=='No')[0]
OncEW_trig_ind = np.where(trig_OncEW=='Yes')[0]
OncEW_only_nan_ind = np.where((trig_OncEW=='No')&(trig_ShakeAlert=='Yes'))[0]
OncEW_only_trig_ind = np.where((trig_OncEW=='Yes')&(trig_ShakeAlert=='No'))[0]
ShakeAlert_nan_ind = np.where(trig_ShakeAlert=='No')[0]
ShakeAlert_trig_ind = np.where(trig_ShakeAlert=='Yes')[0]
ShakeAlert_only_nan_ind = np.where((trig_ShakeAlert=='No')&(trig_OncEW=='Yes'))[0]
ShakeAlert_only_trig_ind = np.where((trig_ShakeAlert=='Yes')&(trig_OncEW=='No'))[0]
both_nan_ind = np.where((trig_ShakeAlert=='No')&(trig_OncEW=='No'))[0]
both_trig_ind = np.where((trig_ShakeAlert=='Yes')&(trig_OncEW=='Yes'))[0]

# Get P-wave pick time residuals
ppick_res_OncEW = df['1/5/4.9_P-arrival Residual (s)']
ppick_res_ShakeAlert = df[f'0.05/5/20_{shakealert_dtype}_P-arrival Residual (s)']

# Get absolute value of pick time residuals
OncEW_res_abs = np.abs(ppick_res_OncEW)
ShakeAlert_res_abs = np.abs(ppick_res_ShakeAlert)

# Get magnitude and distance metrics
repi = df['Repi(km)'].values
rrup = df['Rrup(km)'].values
rhyp = df['Rhyp(km)'].values
mag = df['Magnitude'].values

# What % of records were missed
onc_miss_perc = round(100*len(OncEW_nan_ind)/len(trig_OncEW),2)
sa_miss_perc = round(100*len(ShakeAlert_nan_ind)/len(trig_ShakeAlert),2)

# What % fo records were detected
onc_detec_perc = round(100*len(OncEW_trig_ind)/len(trig_OncEW),2)
sa_detec_perc = round(100*len(ShakeAlert_trig_ind)/len(trig_ShakeAlert),2)

# What % of records have a picktime residual < 1 s?
onc_res1_perc = round(100*len(np.where(OncEW_res_abs<1)[0])/len(np.where(np.isnan(OncEW_res_abs)==False)[0]),2)
sa_res1_perc = round(100*len(np.where(ShakeAlert_res_abs<1)[0])/len(np.where(np.isnan(ShakeAlert_res_abs)==False)[0]),2)


#%%
########## Separate subplots for STA/LTA params ##########

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
labels=[r'(a)',r'(b)',r'(c)',r'(d)',r'(e)',r'(f)']

cNorm  = colors.Normalize(vmin=-10, vmax=10)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('seismic') )

fig, axs = plt.subplots(3,2,figsize=(4,6))

# Missed events
plt.text(0.56, 0.965, r'$\bf{Missed}$ $\bf{Events}$', transform=fig.transFigure, horizontalalignment='center')

axs[0,0].grid(alpha=0.25)
axs[0,0].scatter(rhyp[ShakeAlert_nan_ind],mag[ShakeAlert_nan_ind],alpha=0.65,edgecolors='k',facecolors='none',marker='o',s=15,lw=0.3)
axs[0,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[0,0].xaxis.set_major_locator(MultipleLocator(500))
axs[0,0].xaxis.set_minor_locator(MultipleLocator(250))
axs[0,0].set_xlim(0,1300)
axs[0,0].set_ylim(6,9.7)
axs[0,0].text(-0.2,1.0,labels[0],transform=axs[0,0].transAxes,fontsize=10,va='top',ha='right')
axs[0,0].text(0.02,0.02,f'Missed: {sa_miss_perc}%',transform=axs[0,0].transAxes,fontsize=9,va='bottom',ha='left')
axs[0,0].set_ylabel('Magnitude', fontsize=10, labelpad=21)
axs[0,0].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)

axs[0,1].grid(alpha=0.25)
axs[0,1].scatter(rhyp[OncEW_nan_ind],mag[OncEW_nan_ind],alpha=0.95,edgecolors='k',facecolors='none',marker='o',s=15,lw=0.3)
axs[0,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[0,1].xaxis.set_major_locator(MultipleLocator(500))
axs[0,1].xaxis.set_minor_locator(MultipleLocator(250))
axs[0,1].set_xlim(0,1300)
axs[0,1].set_ylim(6,9.7)
axs[0,1].text(-0.075,1.0,labels[1],transform=axs[0,1].transAxes,fontsize=10,va='top',ha='right')
axs[0,1].text(0.02,0.02,f'Missed: {onc_miss_perc}%',transform=axs[0,1].transAxes,fontsize=9,va='bottom',ha='left')
axs[0,1].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[0,1].set_yticks([])


# Radiation scale pattern
plt.text(0.56, 0.65, r'$\bf{P-wave}$ $\bf{Radiation}$ $\bf{Pattern}$', transform=fig.transFigure, horizontalalignment='center')

axs[1,0].scatter(az[ShakeAlert_trig_ind], P[ShakeAlert_trig_ind], lw=0.3, s=15, marker='o', fc='none', ec='mediumblue', alpha=0.65, label='Triggered')
axs[1,0].scatter(az[ShakeAlert_nan_ind], P[ShakeAlert_nan_ind], lw=0.3, s=15, marker='o', fc='none', ec='k', alpha=0.65, label='Missed')
axs[1,0].grid(alpha=0.25)
axs[1,0].set_ylim(-1.15,1.15)
axs[1,0].set_xlabel('Azimuth')
axs[1,0].text(-0.2,1.0,labels[2],transform=axs[1,0].transAxes,fontsize=10,va='top',ha='right')
axs[1,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[1,0].set_ylabel('Amplification Factor',labelpad=15)

axs[1,1].scatter(az[OncEW_trig_ind], P[OncEW_trig_ind], lw=0.3, s=15, marker='o', fc='none', ec='mediumblue', alpha=0.65, label='Triggered')
axs[1,1].scatter(az[OncEW_nan_ind], P[OncEW_nan_ind], lw=0.3, s=15, marker='o', fc='none', ec='k', alpha=1, label='Missed')
axs[1,1].grid(alpha=0.25)
axs[1,1].set_ylim(-1.15,1.15)
axs[1,1].set_xlabel('Azimuth')
axs[1,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[1,1].text(-0.075,1.0,labels[3],transform=axs[1,1].transAxes,fontsize=10,va='top',ha='right')
axs[1,1].set_yticks([])

# Detected events
plt.text(0.56, 0.34, r'$\bf{Detected}$ $\bf{Events}$', transform=fig.transFigure, horizontalalignment='center')

axs[2,0].scatter(rhyp[ShakeAlert_nan_ind],mag[ShakeAlert_nan_ind],lw=0.3,s=15,marker='o',fc='none',ec=scalarMap.to_rgba(ppick_res_ShakeAlert),alpha=1)
axs[2,0].xaxis.set_major_locator(MultipleLocator(500))
axs[2,0].xaxis.set_minor_locator(MultipleLocator(250))
axs[2,0].tick_params(axis='x', which='minor')
axs[2,0].tick_params(axis='y', which='minor')
axs[2,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[2,0].grid(linestyle='-',alpha=0.25)
# axs[2,0].set_xlim(0,1300)
axs[2,0].set_ylim(ymin=5.9)
# axs[2,0].set_yscale('log')
axs[2,0].text(-0.2,1.0,labels[4],transform=axs[2,0].transAxes,fontsize=10,va='top',ha='right')
axs[2,0].text(0.02,0.02,f'Detected: {sa_detec_perc}%',transform=axs[2,0].transAxes,fontsize=9,va='bottom',ha='left')
axs[2,0].set_ylabel(r'Magnitude', fontsize=10)
axs[2,0].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)

axs[2,1].scatter(rhyp[OncEW_nan_ind],mag[OncEW_nan_ind],lw=0.3,s=15,marker='o',fc='none',ec=scalarMap.to_rgba(ppick_res_OncEW),alpha=1)
axs[2,1].xaxis.set_major_locator(MultipleLocator(500))
axs[2,1].xaxis.set_minor_locator(MultipleLocator(250))
axs[2,1].tick_params(axis='x', which='minor')
axs[2,1].tick_params(axis='y', which='minor')
axs[2,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[2,1].grid(linestyle='-',alpha=0.25)
# axs[2,1].set_xlim(0,1300)
axs[2,1].set_ylim(ymin=5.9)
axs[2,1].text(-0.075,1.0,labels[5],transform=axs[2,1].transAxes,fontsize=10,va='top',ha='right')
axs[2,1].text(0.02,0.02,f'Detected: {onc_detec_perc}%',transform=axs[2,1].transAxes,fontsize=9,va='bottom',ha='left')
axs[2,1].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[2,1].set_yticks([])

cax = fig.add_axes([0.825, 0.14, 0.02, 0.19])
cbar = plt.colorbar(scalarMap, label=r'$\delta_{arr}$ (s)', cax=cax, orientation="vertical")
cbar.ax.tick_params(rotation=0,labelsize=9)

handles = [mlines.Line2D([],[],marker='o',c='k',markerfacecolor='none',markersize=5,lw=0.3,linestyle='None'),
           mlines.Line2D([],[],marker='o',c='mediumblue',markerfacecolor='none',markersize=5,lw=0.3,linestyle='None')]
labels = ['Missed', 'Detected']
           
axs[2,0].legend(handles, labels,
           loc='upper center', bbox_to_anchor=(1.1675, -0.475), ncol=2)  # Adjust position and other parameters

plt.subplots_adjust(left=0.19,right=0.8,top=0.95,bottom=0.15,hspace=0.8,wspace=0.3)
plt.savefig(f'{home}/manuscript/figures/unannotated/sta-lta_results_{shakealert_dtype}_v3.png',dpi=500)
