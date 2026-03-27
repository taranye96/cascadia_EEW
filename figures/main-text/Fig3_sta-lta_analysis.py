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
shakealert_dtype = 'epic'

# Read in dataframes
az_df = pd.read_csv(f'{home}/event_detection/azimuth/cascadia_azimuth-amplification.csv')
pick_df = pd.read_csv(f'{home}/event_detection/sta_lta/{batch}_sta_lta.csv')
arriv_df = pd.read_csv(f'{home}/event_detection/p_waves/{batch}_arrivals.csv')

arriv_df['Run'] = arriv_df['Run'].str.replace('cascadia-', 'cascadia.', regex=False)
pick_df['Run'] = pick_df['Run'].str.replace('cascadia-', 'cascadia.', regex=False)

# Merge dataframes
df = pd.merge(arriv_df, pick_df, on=['Batch','Run','Station'], how='outer')
df = pd.merge(df, az_df, on=['Run', 'Station'], how='outer')

# Get azimuth and P-wave amplification factor
az = df['Azimuth']
P = df['P']

# Get STA/LTA triggers (Yes or No) 
trig_OncEW = df['1/5/4.9_triggered']
trig_ShakeAlert = df[f'0.05/5/20_triggered']

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
ppick_res_ShakeAlert = df[f'0.05/5/20_P-arrival Residual (s)']

# Get absolute value of pick time residuals
OncEW_res_abs = np.abs(ppick_res_OncEW)
ShakeAlert_res_abs = np.abs(ppick_res_ShakeAlert)

# Get magnitude and distance metrics
# repi = arriv_df['Repi(km)'].values
# rrup = arriv_df['Rrup(km)'].values
rhyp = df['Rhyp(km)'].values
mag = df['Magnitude'].values

# What % of records were missed
onc_miss_perc = round(100*len(OncEW_nan_ind)/len(trig_OncEW))
sa_miss_perc = round(100*len(ShakeAlert_nan_ind)/len(trig_ShakeAlert))

# What % fo records were detected
onc_detec_perc = round(100*len(OncEW_trig_ind)/len(trig_OncEW))
sa_detec_perc = round(100*len(ShakeAlert_trig_ind)/len(trig_ShakeAlert))

# What % of records have a picktime residual < 1 s?
onc_res1_perc = round(100*len(np.where(OncEW_res_abs<1)[0])/len(np.where(np.isnan(OncEW_res_abs)==False)[0]))
sa_res1_perc = round(100*len(np.where(ShakeAlert_res_abs<1)[0])/len(np.where(np.isnan(ShakeAlert_res_abs)==False)[0]))

onc_res_0pt5perc = round(100*len(np.where(ppick_res_OncEW>-0.5)[0])/len(np.where(np.isnan(ppick_res_OncEW)==False)[0]))
sa_res_0pt5perc = round(100*len(np.where(ppick_res_ShakeAlert>-0.5)[0])/len(np.where(np.isnan(ppick_res_ShakeAlert)==False)[0]))


#%%
########## Separate subplots for STA/LTA params ##########

from scipy.stats import gaussian_kde
import matplotlib.colors as mcolors

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
labels=['(a)','(b)','(c)']

fig, axs = plt.subplots(3,2,figsize=(3.5,6.5))

# Missed events
fig.text(0.56, 0.96, 'Missed Detections', horizontalalignment='center', fontfamily='DejaVu Sans', fontweight='bold')
xy_miss = np.vstack([rhyp[ShakeAlert_nan_ind],mag[ShakeAlert_nan_ind]])
z_miss = gaussian_kde(xy_miss)(xy_miss)
z_miss_norm = (z_miss - z_miss.min()) / (z_miss.max() - z_miss.min())
norm_miss = mcolors.Normalize(vmin=z_miss_norm.min(), vmax=z_miss_norm.max())
cmap_miss = plt.cm.Reds
fc_miss = cmap_miss(norm_miss(z_miss_norm))
axs[0,0].minorticks_on()
axs[0,0].tick_params(which='minor', bottom=True, left=True)
axs[0,0].scatter(rhyp[ShakeAlert_nan_ind],mag[ShakeAlert_nan_ind],alpha=1,edgecolors='k',facecolors=fc_miss,marker='o',s=15,lw=0.2)
axs[0,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[0,0].xaxis.set_major_locator(MultipleLocator(500))
axs[0,0].xaxis.set_minor_locator(MultipleLocator(250))
axs[0,0].yaxis.set_major_locator(MultipleLocator(1))
axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.5))
axs[0,0].grid(alpha=0)
axs[0,0].set_xscale('log')
axs[0,0].set_xlim(5,1300)
axs[0,0].set_ylim(5.9,9.99)
axs[0,0].text(-0.15,1.1,labels[0],transform=axs[0,0].transAxes,fontsize=10,va='top',ha='right')
axs[0,0].set_ylabel('Magnitude', fontsize=10, labelpad=17)
axs[0,0].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[0,0].text(0.03,0.02,f'Missed: {sa_miss_perc}%',transform=axs[0,0].transAxes,fontsize=9,va='bottom',ha='left')
axs[0,0].text(0.02,0.96,'EPIC',transform=axs[0,0].transAxes,fontsize=9,va='top',ha='left')

xy_miss = np.vstack([rhyp[OncEW_nan_ind],mag[OncEW_nan_ind]])
z_miss = gaussian_kde(xy_miss)(xy_miss)
z_miss_norm = (z_miss - z_miss.min()) / (z_miss.max() - z_miss.min())
norm_miss = mcolors.Normalize(vmin=z_miss_norm.min(), vmax=z_miss_norm.max())
cmap_miss = plt.cm.Reds
fc_miss = cmap_miss(norm_miss(z_miss_norm))
axs[0,1].minorticks_on()
axs[0,1].tick_params(which='minor', bottom=True, left=True)
axs[0,1].scatter(rhyp[OncEW_nan_ind],mag[OncEW_nan_ind],alpha=1,edgecolors='k',facecolors=fc_miss,marker='o',s=15,lw=0.2)
axs[0,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[0,1].xaxis.set_major_locator(MultipleLocator(500))
axs[0,1].xaxis.set_minor_locator(MultipleLocator(250))
axs[0,1].yaxis.set_major_locator(MultipleLocator(1))
axs[0,1].yaxis.set_minor_locator(MultipleLocator(0.5))
axs[0,1].grid(alpha=0)
axs[0,1].set_xscale('log')
axs[0,1].set_xlim(5,1300)
axs[0,1].set_ylim(5.9,9.99)
axs[0,1].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[0,1].set_yticklabels([])
axs[0,1].text(0.03,0.02,f'Missed: {onc_miss_perc}%',transform=axs[0,1].transAxes,fontsize=9,va='bottom',ha='left')
axs[0,1].text(0.02,0.96,'ONC-EW',transform=axs[0,1].transAxes,fontsize=9,va='top',ha='left')

# Radiation scale pattern
fig.text(0.56, 0.67, 'P-wave Radiation Pattern', horizontalalignment='center', fontfamily='DejaVu Sans', fontweight='bold')
xy_miss = np.vstack([rhyp[ShakeAlert_nan_ind],P[ShakeAlert_nan_ind]])
z_miss = gaussian_kde(xy_miss)(xy_miss)
z_miss_norm = (z_miss - z_miss.min()) / (z_miss.max() - z_miss.min())
norm_miss = mcolors.Normalize(vmin=z_miss_norm.min(), vmax=z_miss_norm.max())
cmap_miss = plt.cm.Reds
fc_miss = cmap_miss(norm_miss(z_miss_norm))
xy_detec = np.vstack([rhyp[ShakeAlert_trig_ind],P[ShakeAlert_trig_ind]])
z_detec = gaussian_kde(xy_detec)(xy_detec)
z_detec_norm = (z_detec - z_detec.min()) / (z_detec.max() - z_detec.min())
norm_detec = mcolors.Normalize(vmin=z_detec_norm.min(), vmax=z_detec_norm.max())
cmap_detec = plt.cm.Blues
fc_detec = cmap_detec(norm_detec(z_detec_norm))
axs[1,0].minorticks_on()
axs[1,0].tick_params(which='minor', bottom=True, left=True)
axs[1,0].scatter(az[ShakeAlert_trig_ind], P[ShakeAlert_trig_ind], lw=0.15, s=17, marker='^', fc=fc_detec, ec='k', alpha=1, label='Triggered')
axs[1,0].scatter(az[ShakeAlert_nan_ind], P[ShakeAlert_nan_ind], lw=0.2, s=15, marker='o', fc=fc_miss, ec='k', alpha=1, label='Missed')
axs[1,0].set_ylim(-1.25,1.45)
axs[1,0].set_xlabel('Azimuth (deg)')
axs[1,0].text(-0.15,1.1,labels[1],transform=axs[1,0].transAxes,fontsize=10,va='top',ha='right')
axs[1,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[1,0].set_ylabel('Radiation factor',labelpad=10)
axs[1,0].text(0.02,0.96,'EPIC',transform=axs[1,0].transAxes,fontsize=9,va='top',ha='left')
axs[1,0].yaxis.set_major_locator(MultipleLocator(1))
axs[1,0].yaxis.set_minor_locator(MultipleLocator(0.5))
axs[1,0].xaxis.set_major_locator(MultipleLocator(100))
axs[1,0].xaxis.set_minor_locator(MultipleLocator(50))
axs[1,0].grid(alpha=0)

xy_miss = np.vstack([rhyp[OncEW_nan_ind],P[OncEW_nan_ind]])
z_miss = gaussian_kde(xy_miss)(xy_miss)
z_miss_norm = (z_miss - z_miss.min()) / (z_miss.max() - z_miss.min())
norm_miss = mcolors.Normalize(vmin=z_miss_norm.min(), vmax=z_miss_norm.max())
cmap_miss = plt.cm.Reds
fc_miss = cmap_miss(norm_miss(z_miss_norm))
xy_detec = np.vstack([rhyp[OncEW_trig_ind],P[OncEW_trig_ind]])
z_detec = gaussian_kde(xy_detec)(xy_detec)
z_detec_norm = (z_detec - z_detec.min()) / (z_detec.max() - z_detec.min())
norm_detec = mcolors.Normalize(vmin=z_detec_norm.min(), vmax=z_detec_norm.max())
cmap_detec = plt.cm.Blues
fc_detec = cmap_detec(norm_detec(z_detec_norm))
axs[1,1].minorticks_on()
axs[1,1].tick_params(which='minor', bottom=True, left=True)
axs[1,1].scatter(az[OncEW_trig_ind], P[OncEW_trig_ind], lw=0.15, s=17, marker='^', fc=fc_detec, ec='k', alpha=1, label='Triggered')
axs[1,1].scatter(az[OncEW_nan_ind], P[OncEW_nan_ind], lw=0.2, s=15, marker='o', fc=fc_miss, ec='k', alpha=1, label='Missed')
axs[1,1].set_ylim(-1.25,1.45)
axs[1,1].set_xlabel('Azimuth (deg)')
axs[1,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[1,1].set_yticklabels([])
axs[1,1].text(0.02,0.96,'ONC-EW',transform=axs[1,1].transAxes,fontsize=9,va='top',ha='left')
axs[1,1].yaxis.set_major_locator(MultipleLocator(1))
axs[1,1].yaxis.set_minor_locator(MultipleLocator(0.5))
axs[1,1].xaxis.set_major_locator(MultipleLocator(100))
axs[1,1].xaxis.set_minor_locator(MultipleLocator(50))
axs[1,1].grid(alpha=0)

# Detected events
fig.text(0.56, 0.38, 'Successful Detections', horizontalalignment='center', fontfamily='DejaVu Sans', fontweight='bold')
ind = np.where(np.isnan(ppick_res_ShakeAlert)==False)[0]
xy = np.vstack([rhyp[ind],ppick_res_ShakeAlert[ind]])
z = gaussian_kde(xy)(xy)
z_norm = (z - z.min()) / (z.max() - z.min())
norm = mcolors.Normalize(vmin=z_norm.min(), vmax=z_norm.max())
cmap = plt.cm.Blues
facecolors = cmap(norm(z_norm))
axs[2,0].minorticks_on()
axs[2,0].tick_params(which='minor', bottom=True, left=True)
axs[2,0].scatter(rhyp[ind],ppick_res_ShakeAlert[ind],lw=0.15,s=17,marker='^',fc=facecolors,ec='k',alpha=1)
axs[2,0].axhline(0,0,ls='--',lw=1,color='k')
axs[2,0].xaxis.set_major_locator(MultipleLocator(500))
axs[2,0].xaxis.set_minor_locator(MultipleLocator(250))
axs[2,0].yaxis.set_major_locator(MultipleLocator(2))
axs[2,0].yaxis.set_minor_locator(MultipleLocator(1))
axs[2,0].grid(alpha=0)
axs[2,0].tick_params(axis='x', which='minor')
axs[2,0].tick_params(axis='y', which='minor')
axs[2,0].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[2,0].set_xlim(5,1300)
axs[2,0].set_ylim(-5,1)
axs[2,0].text(-0.15,1.1,labels[2],transform=axs[2,0].transAxes,fontsize=10,va='top',ha='right')
axs[2,0].set_ylabel(r'$\delta_{arr}$ (s)', fontsize=10, labelpad=10)
axs[2,0].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[2,0].text(0.02,0.96,'EPIC',transform=axs[2,0].transAxes,fontsize=9,va='top',ha='left')
axs[2,0].text(0.03,0.02,r'$\delta_{arr}$'+f' < 0.5 s: {sa_res_0pt5perc}%',
              transform=axs[2,0].transAxes,fontsize=9,va='bottom',ha='left')
axs[2,0].set_xscale('log')

ind = np.where(np.isnan(ppick_res_OncEW)==False)[0]
xy = np.vstack([rhyp[ind],ppick_res_OncEW[ind]])
z = gaussian_kde(xy)(xy)
z_norm = (z - z.min()) / (z.max() - z.min())
norm = mcolors.Normalize(vmin=z_norm.min(), vmax=z_norm.max())
cmap = plt.cm.Blues
facecolors = cmap(norm(z_norm))
axs[2,1].minorticks_on()
axs[2,1].tick_params(which='minor', bottom=True, left=True)
axs[2,1].scatter(rhyp[ind],ppick_res_OncEW[ind],lw=0.15,s=17,marker='^',fc=facecolors,ec='k',alpha=1)
axs[2,1].axhline(0,0,ls='--',lw=1,color='k')
axs[2,1].xaxis.set_major_locator(MultipleLocator(500))
axs[2,1].xaxis.set_minor_locator(MultipleLocator(250))
axs[2,1].yaxis.set_major_locator(MultipleLocator(2))
axs[2,1].yaxis.set_minor_locator(MultipleLocator(1))
axs[2,1].grid(alpha=0)
axs[2,1].tick_params(axis='x', which='minor')
axs[2,1].tick_params(axis='y', which='minor')
axs[2,1].tick_params(direction="out",labelright=False,bottom=True,left=True,labelsize=10)
axs[2,1].set_xlim(5,1300)
axs[2,1].set_ylim(-5, 1)
axs[2,1].set_xlabel(r'$R_{hyp}$ (km)', fontsize=10)
axs[2,1].set_yticklabels([])
axs[2,1].text(0.02,0.96,'ONC-EW',transform=axs[2,1].transAxes,fontsize=9,va='top',ha='left')
axs[2,1].text(0.03,0.02,r'$\delta_{arr}$'+f' < 0.5 s: {onc_res_0pt5perc}%',
              transform=axs[2,1].transAxes,fontsize=9,va='bottom',ha='left')
axs[2,1].set_xscale('log')

# Colorbars
sm1 = plt.cm.ScalarMappable(norm=norm_miss, cmap=cmap_miss)
sm1.set_array([])
cax1 = fig.add_axes([0.625, 0.1, 0.325, 0.015])  # left, bottom, width, height
cbar1 = fig.colorbar(sm1, cax1,
                    shrink=1,   
                    pad=0.15,
                    anchor=(0.0, 1),
                    orientation='horizontal')   
cbar1.set_ticklabels([])

sm2 = plt.cm.ScalarMappable(norm=norm_detec, cmap=cmap_detec)
sm2.set_array([])
cax2 = fig.add_axes([0.625, 0.085, 0.325, 0.015])  # left, bottom, width, height
cbar2 = fig.colorbar(sm2, cax2,
                    shrink=1,   
                    pad=0.15,
                    anchor=(0.0, 1),
                    orientation='horizontal')   
cbar2.set_label('Relative point density', labelpad=5)   
cbar2.set_ticks([0, 1])
cbar2.set_ticklabels(['Low', 'High'])


handles = [mlines.Line2D([],[],marker='o',c='k',markerfacecolor='orangered',markersize=5.5,markeredgewidth=0.65,linestyle='none'),
           mlines.Line2D([],[],marker='^',c='k',markerfacecolor='dodgerblue',markersize=5.5,markeredgewidth=0.65,linestyle='none')]
labels = ['Missed', 'Detected']
           
axs[2,0].legend(handles, labels,
           loc='upper center', bbox_to_anchor=(0.5, -0.45), ncol=1)  # Adjust position and other parameters

plt.subplots_adjust(left=0.175,right=0.975,top=0.95,bottom=0.2,hspace=0.7,wspace=0.1)
plt.savefig(f'{home}/manuscript/figures/main-text/Fig3_sta-lta_results.png',dpi=500)
