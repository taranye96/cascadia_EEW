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

# Name of batch for simulations
batch = 'cascadia'

# Are we wanting the  'obspy' or 'elarmS' picker for the STA/LTA results using
    # EPIC/ShakeAlert parameters?
shakealert_dtype = 'obspy'

# Get azimuth and P-wave amplification factor
az_df = pd.read_csv('/Users/tnye/ONC/event_detection/azimuth/cascadia_azimuth-amplification.csv')
az = az_df['Azimuth']
P = az_df['P']

# Read in STA/LTA file
pick_file = f'/Users/tnye/ONC/event_detection/sta_lta/{batch}_sta_lta_polar.csv'
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
labels=[r'$\bf{(a)}$',r'$\bf{(b)}$',r'$\bf{(c)}$',r'$\bf{(d)}$',r'$\bf{(e)}$',r'$\bf{(f)}$']

cNorm  = colors.Normalize(vmin=6.5, vmax=9.5)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('viridis') )
cNorm2  = colors.Normalize(vmin=0, vmax=1200)
scalarMap2 = cmx.ScalarMappable(norm=cNorm2, cmap=plt.get_cmap('plasma') )

fig, axs = plt.subplots(3,2,figsize=(6.5,7))

# Missed events
plt.text(0.5, 0.975, r'$\bf{Missed}$ $\bf{Events}$', transform=fig.transFigure, horizontalalignment='center')

axs[0,0].grid(alpha=0.25)
axs[0,0].scatter(rhyp[ShakeAlert_nan_ind],mag[ShakeAlert_nan_ind],alpha=0.5,edgecolors='dimgray',facecolors='none',marker='o',s=10)
axs[0,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[0,0].xaxis.set_major_locator(MultipleLocator(500))
axs[0,0].xaxis.set_minor_locator(MultipleLocator(100))
axs[0,0].set_xlim(0,1300)
axs[0,0].set_ylim(6,9.7)
axs[0,0].text(-0.25,1.0,labels[0],transform=axs[0,0].transAxes,fontsize=10,va='top',ha='right')
axs[0,0].text(0.02,0.02,f'Missed: {sa_miss_perc}%',transform=axs[0,0].transAxes,fontsize=9,va='bottom',ha='left')
axs[0,0].set_ylabel('Magnitude',fontsize=10)
axs[0,0].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)

axs[0,1].grid(alpha=0.25)
axs[0,1].scatter(rhyp[OncEW_nan_ind],mag[OncEW_nan_ind],alpha=0.5,edgecolors='dimgray',facecolors='none',marker='o',s=10)
axs[0,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[0,1].xaxis.set_major_locator(MultipleLocator(500))
axs[0,1].xaxis.set_minor_locator(MultipleLocator(100))
axs[0,1].set_xlim(0,1300)
axs[0,1].set_ylim(6,9.7)
axs[0,1].text(-0.25,1.0,labels[1],transform=axs[0,1].transAxes,fontsize=10,va='top',ha='right')
axs[0,1].text(0.02,0.02,f'Missed: {onc_miss_perc}%',transform=axs[0,1].transAxes,fontsize=9,va='bottom',ha='left')
axs[0,1].set_ylabel('Magnitude',fontsize=10)
axs[0,1].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)


# Radiation scale pattern
plt.text(0.5, 0.67, r'$\bf{P-wave}$ $\bf{Radiation}$ $\bf{Pattern}$', transform=fig.transFigure, horizontalalignment='center')

axs[1,0].scatter(az[ShakeAlert_trig_ind], P[ShakeAlert_trig_ind], lw=0.8, s=20, marker='+', c=scalarMap2.to_rgba(rhyp[ShakeAlert_trig_ind]), alpha=0.5, label='Triggered')
axs[1,0].scatter(az[ShakeAlert_nan_ind], P[ShakeAlert_nan_ind], lw=0.8, s=20, marker='o', edgecolors=scalarMap2.to_rgba(rhyp[ShakeAlert_nan_ind]), facecolors='none', alpha=0.5, label='Missed')
axs[1,0].grid(alpha=0.25)
axs[1,0].set_ylim(-1.15,1.15)
axs[1,0].set_xlabel('Azimuth')
axs[1,0].text(-0.25,1.0,labels[2],transform=axs[1,0].transAxes,fontsize=10,va='top',ha='right')
axs[1,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[1,0].set_ylabel('Amplification Factor')

axs[1,1].scatter(az[OncEW_trig_ind], P[OncEW_trig_ind], lw=0.8, s=20, marker='+', c=scalarMap2.to_rgba(rhyp[OncEW_trig_ind]), alpha=0.5, label='Triggered')
axs[1,1].scatter(az[OncEW_nan_ind], P[OncEW_nan_ind], lw=0.8, s=20, marker='o', edgecolors=scalarMap2.to_rgba(rhyp[OncEW_nan_ind]), facecolors='none', alpha=0.5, label='Missed')
axs[1,1].grid(alpha=0.25)
axs[1,1].set_ylim(-1.15,1.15)
axs[1,1].set_xlabel('Azimuth')
axs[1,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[1,1].text(-0.25,1.0,labels[3],transform=axs[1,1].transAxes,fontsize=10,va='top',ha='right')
axs[1,1].set_ylabel('Amplification Factor')

cax = fig.add_axes([0.875, 0.455, 0.02, 0.191])
cbar = plt.colorbar(scalarMap2, label='$R_{rhyp}$ (km)', cax=cax, ticks=[0,400,800,1200], orientation="vertical")
cbar.ax.tick_params(rotation=0,labelsize=9)


# Detected events
plt.text(0.5, 0.36, r'$\bf{Detected}$ $\bf{Events}$', transform=fig.transFigure, horizontalalignment='center')

axs[2,0].scatter(rhyp,ppick_res_ShakeAlert,c=scalarMap.to_rgba(mag),marker='+',lw=0.5,alpha=0.4,s=25)
axs[2,0].axhline(0,0,ls='--',lw=0.8,color='k')
axs[2,0].xaxis.set_major_locator(MultipleLocator(500))
axs[2,0].xaxis.set_minor_locator(MultipleLocator(100))
axs[2,0].tick_params(axis='x', which='minor')
axs[2,0].tick_params(axis='y', which='minor')
axs[2,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[2,0].grid(linestyle='-',alpha=0.25)
axs[2,0].set_xlim(0,1300)
axs[2,0].text(-0.25,1.0,labels[4],transform=axs[2,0].transAxes,fontsize=10,va='top',ha='right')
axs[2,0].text(0.02,0.02,f'Detected: {sa_detec_perc}%',transform=axs[2,0].transAxes,fontsize=9,va='bottom',ha='left')
axs[2,0].set_ylabel(r'$\delta_{arr}$ (s)', fontsize=10)
axs[2,0].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)

axs[2,1].scatter(rhyp,ppick_res_OncEW,c=scalarMap.to_rgba(mag),marker='+',lw=0.5,alpha=0.4,s=25)
axs[2,1].axhline(0,0,ls='--',lw=0.8,color='k')
axs[2,1].xaxis.set_major_locator(MultipleLocator(500))
axs[2,1].xaxis.set_minor_locator(MultipleLocator(100))
axs[2,1].tick_params(axis='x', which='minor')
axs[2,1].tick_params(axis='y', which='minor')
axs[2,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[2,1].grid(linestyle='-',alpha=0.25)
axs[2,1].set_xlim(0,1300)
axs[2,1].text(-0.25,1.0,labels[5],transform=axs[2,1].transAxes,fontsize=10,va='top',ha='right')
axs[2,1].text(0.02,0.02,f'Detected: {onc_detec_perc}%',transform=axs[2,1].transAxes,fontsize=9,va='bottom',ha='left')
axs[2,1].set_ylabel(r'$\delta_{arr}$ (s)', fontsize=10)
axs[2,1].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)

cax = fig.add_axes([0.875, 0.15, 0.02, 0.19])
cbar = plt.colorbar(scalarMap, label='Magnitude', cax=cax, ticks=[6.5,7.5,8.5,9.5], orientation="vertical")
cbar.ax.tick_params(rotation=0,labelsize=9)

handles = [mlines.Line2D([],[],marker='o',c='k',markerfacecolor='none',markersize=5,linestyle='None'),
           mlines.Line2D([],[],marker='+',lw=0.5,markersize=5,c='k',linestyle='None')]
labels = ['Missed', 'Detected']
           
axs[2,0].legend(handles, labels,
           loc='upper center', bbox_to_anchor=(1.175, -0.35), ncol=2)  # Adjust position and other parameters

plt.subplots_adjust(left=0.12,right=0.86,top=0.95,bottom=0.15,hspace=0.6,wspace=0.425)
plt.savefig(f'/Users/tnye/ONC/manuscript/figures/unannotated/sta-lta_results_{shakealert_dtype}.png',dpi=400)
