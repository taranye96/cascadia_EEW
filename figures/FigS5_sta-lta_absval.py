#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 12:25:43 2024

@author: tnye
"""

##############################################################################
# This script makes Figure S5 in supporting information, which plots the
# absolute value of the STA/LTA pick time residuals.
##############################################################################

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import matplotlib.cm as cmx

# Name of batch for simulations
batch = 'cascadia'

# Are we wanting the  'obspy' or 'elarmS' picker for the STA/LTA results using
    # EPIC/ShakeAlert parameters?
shakealert_dtype = 'obspy'

# Get azimuth and P-wave amplification factor
az_df = pd.read_csv('/Users/tnye/ONC/event_detection/cascadia_azimuth-amplification.csv')
az = az_df['Azimuth']
P = az_df['P']

# Read in STA/LTA file
pick_file = f'/Users/tnye/ONC/event_detection/{batch}_sta_lta_polar.csv'
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

cNorm  = colors.Normalize(vmin=6.5, vmax=9.5)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('viridis') )
cNorm2  = colors.Normalize(vmin=0, vmax=1200)
scalarMap2 = cmx.ScalarMappable(norm=cNorm2, cmap=plt.get_cmap('plasma') )

fig, axs = plt.subplots(1,2,figsize=(6.5,3))

axs[0].scatter(rhyp,OncEW_res_abs,c=scalarMap.to_rgba(mag),marker='+',lw=0.5,alpha=0.4,s=25)
axs[0].set_yscale('log')
axs[0].set_xscale('log')
axs[0].set_xlim(8,1.5*10**3)
axs[0].set_ylim(10**-6,10**3)
axs[0].xaxis.set_minor_formatter(ticker.NullFormatter())
axs[0].tick_params(axis='x', which='minor', bottom=True)
axs[0].tick_params(axis='y', which='minor', left=True)
axs[0].tick_params(axis='x', which='major', bottom=True, top=False)
axs[0].tick_params(axis='y', which='major', left=True, right=False)
axs[0].tick_params(direction="out",labelright=False,labelsize=10)
axs[0].grid(linestyle='-',alpha=0.25)
axs[0].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[0].set_ylabel(r'| $\delta_{arr}$ | (s)', fontsize=10)
axs[0].text(0.02,0.02,r'|$\delta_{arr}$| < 1 s: '+ f'{onc_res1_perc}%',transform=axs[0].transAxes,fontsize=9,va='bottom',ha='left')

axs[1].scatter(rhyp,ShakeAlert_res_abs,c=scalarMap.to_rgba(mag),marker='+',lw=0.5,alpha=0.4,s=25)
axs[1].set_yscale('log')
axs[1].set_xscale('log')
axs[1].set_xlim(8,1.5*10**3)
axs[1].set_ylim(10**-6,10**3)
axs[1].xaxis.set_minor_formatter(ticker.NullFormatter())
axs[1].tick_params(axis='x', which='minor', bottom=True)
axs[1].tick_params(axis='y', which='minor', left=True)
axs[1].tick_params(axis='x', which='major', bottom=True, top=False)
axs[1].tick_params(axis='y', which='major', left=True, right=False)
axs[1].tick_params(direction="out",labelright=False,labelsize=10)
axs[1].grid(linestyle='-',alpha=0.25)
axs[1].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[1].set_ylabel(r'| $\delta_{arr}$ | (s)', fontsize=10)
axs[1].text(0.02,0.02,r'|$\delta_{arr}$| < 1 s: '+ f'{sa_res1_perc}%',transform=axs[1].transAxes,fontsize=9,va='bottom',ha='left')

cax = fig.add_axes([0.125, 0.15, 0.825, 0.07])
cbar = fig.colorbar(scalarMap, label='Magnitude', cax=cax, ticks=[6.5,7.5,8.5,9.5], orientation="horizontal")

plt.subplots_adjust(left=0.12,right=0.95,top=0.95,bottom=0.45,hspace=0.6,wspace=0.425)
plt.savefig('/Users/tnye/ONC/manuscript/figures/unannotated/sta-lta_absvalue.png',dpi=400)