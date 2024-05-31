#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 13:58:12 2022

@author: tnye
"""

###############################################################################
# This script makes Figures 5–7 in the paper, which show the residuals
# between simulated intensity measures and those estiamted from GMMs.
###############################################################################

# Imports
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Name of batch for simulations
batch = 'cascadia'

# Description of station types used for naming of IM files
station_types = ['onc-onshore', 'onc-offshore','pnsn']

# Read in list of ruptures
ruptures = np.genfromtxt(f'/Users/tnye/ONC/simulations/{batch}/data/ruptures.list',dtype=str)

# Initialize lists
pga_res_nga = np.array([])
pgv_res_nga = np.array([])
MMI_res_pga_nga = np.array([])
MMI_res_pgv_nga = np.array([])
pga_res_bchydro = np.array([])
MMI_res_pga_bchydro = np.array([])
mag_list = np.array([])
rrup_list = np.array([])

pgd_mag_list = np.array([])
pgd_rhyp_list = np.array([])
pgd_res_melgar_nonoise = np.array([])
pgd_res_gold_obs_nonoise = np.array([])
pgd_res_gold_joint_nonoise = np.array([])
pgd_res_melgar_10p = np.array([])
pgd_res_gold_obs_10p = np.array([])
pgd_res_gold_joint_10p = np.array([])
pgd_res_melgar_50p = np.array([])
pgd_res_gold_obs_50p = np.array([])
pgd_res_gold_joint_50p = np.array([])
pgd_res_melgar_90p = np.array([])
pgd_res_gold_obs_90p = np.array([])
pgd_res_gold_joint_90p = np.array([])

rupt_list = np.array([])
stn_list = np.array([])

nga_pga_sd_list = np.array([])
nga_pgv_sd_list = np.array([])
nga_mmi_pga_sd_list = np.array([])
nga_mmi_pgv_sd_list = np.array([])
bchydro_pga_sd_list = np.array([])
bchydro_mmi_pga_sd_list = np.array([])
melgar_sd_list = np.array([])
gold_sd_obs_list = np.array([])
gold_sd_joint_list = np.array([])

# Loop over ruptures
for rupture in ruptures:
    
    # Read in .log file and get magnitude
    run = rupture.replace('.rupt','')
    logfile = f"/Users/tnye/ONC/simulations/{batch}/output/ruptures/{rupture.replace('.rupt','.log')}"
    f = open(logfile, 'r')
    lines = f.readlines()
    mag = float(lines[15].split(' ')[-1].split('\n')[0])
    
    # Loop over station types    
    for stations in station_types:
        
        # Read in intensity measure file
        file = f'/Users/tnye/ONC/groundmotion_analysis/IMs/{batch}/{run}_{stations}.csv'
        df = pd.read_csv(file)
        
        # Get distance and magnitude
        rrup_list = np.append(rrup_list,np.array(df['Rrup(km)']))
        mag_list = np.append(mag_list,np.array(df['Mw']))
        
        # Get NGA-Sub residuals
        pga_res_nga = np.append(pga_res_nga, np.array(df['lnPGA_res_w/noise_NGASub']))
        pgv_res_nga = np.append(pgv_res_nga, np.array(df['lnPGV_res_w/noise_NGASub']))
        MMI_res_pga_nga = np.append(MMI_res_pga_nga, np.array(df['MMI_res_w/noise_NGASub_pga']))
        MMI_res_pgv_nga = np.append(MMI_res_pgv_nga, np.array(df['MMI_res_w/noise_NGASub_pgv']))
        nga_pga_sd_list = np.append(nga_pga_sd_list, np.array(df['NGA_Sub_sd_g']))
        nga_pgv_sd_list = np.append(nga_pga_sd_list, np.array(df['NGA_Sub_sd_cm/s']))
        
        # Get BCHydro residuals
        pga_res_bchydro = np.append(pga_res_bchydro, np.array(df['lnPGA_res_w/noise_BCHydro']))
        MMI_res_pga_bchydro = np.append(MMI_res_pga_bchydro, np.array(df['MMI_res_w/noise_BCHydro_pga']))
        bchydro_pga_sd_list = np.append(bchydro_pga_sd_list, np.array(df['BCHydro_sd_g']))
        
        # Get station names and ruptures
        stn_list = np.append(stn_list,df['Station'].values)
        rupt_list = np.append(rupt_list,[rupture]*len(df))
        
        # Fore onshore GNSS stations, get IMs from different noise levels
        if stations == 'onc-onshore':
            gnss_df_10p = pd.read_csv(file.replace(stations,f'{stations}_gnss-10p'))
            gnss_df_50p = pd.read_csv(file.replace(stations,f'{stations}_gnss-50p'))
            gnss_df_90p = pd.read_csv(file.replace(stations,f'{stations}_gnss-90p'))
            pgd_res_gold_obs_10p = np.append(pgd_res_gold_obs_10p, np.array(gnss_df_10p['lnPGD_noise_res_GA21-obs']))
            pgd_res_gold_joint_10p = np.append(pgd_res_gold_joint_10p, np.array(gnss_df_10p['lnPGD_noise_res_GA21-joint']))
            pgd_res_gold_obs_50p = np.append(pgd_res_gold_obs_50p, np.array(gnss_df_50p['lnPGD_noise_res_GA21-obs']))
            pgd_res_gold_joint_50p = np.append(pgd_res_gold_joint_50p, np.array(gnss_df_50p['lnPGD_noise_res_GA21-joint']))
            pgd_res_gold_obs_90p = np.append(pgd_res_gold_obs_90p, np.array(gnss_df_90p['lnPGD_noise_res_GA21-obs']))
            pgd_res_gold_joint_90p = np.append(pgd_res_gold_joint_90p, np.array(gnss_df_90p['lnPGD_noise_res_GA21-joint']))
            pgd_res_gold_obs_nonoise = np.append(pgd_res_gold_obs_nonoise, np.array(gnss_df_50p['lnPGD_res_GA21-obs']))
            pgd_res_gold_joint_nonoise = np.append(pgd_res_gold_joint_nonoise, np.array(gnss_df_50p['lnPGD_res_GA21-joint']))
            pgd_rhyp_list = np.append(pgd_rhyp_list, np.array(gnss_df_50p['Rhyp(km)']))
            pgd_mag_list = np.append(pgd_mag_list, np.array(gnss_df_50p['Mw']))
            gold_sd_obs_list = np.append(gold_sd_obs_list, np.array(gnss_df_50p['GA21-obs_sd_m']))
            gold_sd_joint_list = np.append(gold_sd_joint_list, np.array(gnss_df_50p['GA21-joint_sd_m']))

# Get maximum standard deviations from the GMMs
pga_sd = np.max(np.concatenate((nga_pga_sd_list,bchydro_pga_sd_list)))
pgv_sd = np.max(nga_pgv_sd_list)           
pgd_sd = np.max(np.concatenate((gold_sd_obs_list,gold_sd_joint_list)))  


#%% PGA and PGV Residuals

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
fig, axs = plt.subplots(2,2,figsize=(6.5,4))

# PGA vs Magnitude
axs[0,0].scatter(mag_list, pga_res_bchydro, label='BCHydro', s=10, marker='^', facecolors='none', edgecolors='goldenrod', lw=0.3, alpha=0.6)
axs[0,0].scatter(mag_list, pga_res_nga, label='NGA-SUB', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[0,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o')
axs[0,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[0,0].axhline(np.log(pga_sd),lw=0.7,ls='--',c='k',label=r'GMM $1\sigma$')
axs[0,0].axhline(-np.log(pga_sd),lw=0.7,ls='--',c='k')
axs[0,0].axhline(0,lw=0.7,ls='-',c='k')
axs[0,0].grid('-',lw=0.7,alpha=0.5)
axs[0,0].set_ylim(-4.5,3)
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[0,0].text(-0.15, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs[0,0].transAxes)
axs[0,0].set_ylabel(r'$\delta_{PGA}$ $ln\left(\frac{pred}{sim}\right)$',fontsize=10)

# PGA vs Distance
axs[0,1].scatter(rrup_list, pga_res_bchydro, label='NGA-Sub', s=10, marker='^', facecolors='none', edgecolors='goldenrod', lw=0.3, alpha=0.6)
axs[0,1].scatter(rrup_list, pga_res_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((rrup_list,rrup_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((rrup_list,rrup_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[0,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[0,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[0,1].set_xscale('log')
axs[0,1].grid('-',lw=0.7,alpha=0.5)
axs[0,1].axhline(0,lw=0.7,ls='-',c='k')
axs[0,1].set_ylim(-4.5,3)
axs[0,1].axhline(np.log(pga_sd),lw=0.7,ls='--',c='k')
axs[0,1].axhline(-np.log(pga_sd),lw=0.7,ls='--',c='k')
axs[0,1].text(-0.1, 1, r'$\bf{(b)}$', ha='right', va='top', transform=axs[0,1].transAxes)

# PGV vs Magnitude
axs[1,0].scatter(mag_list, pgv_res_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(mag_list,pgv_res_nga,statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[1,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[1,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[1,0].grid('-',lw=0.7,alpha=0.5)
axs[1,0].set_ylim(-4.5,3)
axs[1,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[1,0].axhline(0,lw=0.7,ls='-',c='k')
axs[1,0].axhline(np.log(pgv_sd),lw=0.7,ls='--',c='k')
axs[1,0].axhline(-np.log(pgv_sd),lw=0.7,ls='--',c='k')
axs[1,0].text(-0.15, 1, r'$\bf{(c)}$', ha='right', va='top', transform=axs[1,0].transAxes)
axs[1,0].set_xlabel('Magnitude')
axs[1,0].set_ylabel(r'$\delta_{PGV}$ $ln\left(\frac{pred}{sim}\right)$',fontsize=10)

# PGV vs Distance
axs[1,1].scatter(rrup_list, pgv_res_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(rrup_list, pgv_res_nga, statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((rrup_list,rrup_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[1,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[1,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[1,1].set_xscale('log')
axs[1,1].grid('-',lw=0.7,alpha=0.5)
axs[1,1].set_ylim(-4.5,3)
axs[1,1].axhline(0,lw=0.7,ls='-',c='k')
axs[1,1].axhline(np.log(pgv_sd),lw=0.7,ls='--',c='k')
axs[1,1].axhline(-np.log(pgv_sd),lw=0.7,ls='--',c='k')
axs[1,1].text(-0.1, 1, r'$\bf{(d)}$', ha='right', va='top', transform=axs[1,1].transAxes)
axs[1,1].set_xlabel(r'$R_{hyp}$ $(km)$')

# Make legend symbols opaque
handles, labels = axs[0,0].get_legend_handles_labels()
for lh in handles: 
    try:
        lh.set_alpha(1)
    except:
        continue

# Set legend
axs[1,0].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.04,-0.405),facecolor='white',
                frameon=True,fontsize=10,title_fontsize=10,ncol=4,markerscale=2)

# Adjust spacing of subplots
plt.subplots_adjust(right=0.98,left=0.1,top=0.95,bottom=0.25,hspace=0.3,wspace=0.25)

# Save Figure
plt.savefig('/Users/tnye/ONC/manuscript/figures/Fig5_PGA-PGV_GMM_residuals.png',dpi=300)


#%% MMI Residuals

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
fig, axs = plt.subplots(2,2,figsize=(6.5,4))

# MMI_pga vs Magnitude
axs[0,0].scatter(mag_list, MMI_res_pga_bchydro, label='BCHydro', s=10, marker='^', facecolors='none', edgecolors='goldenrod', lw=0.3, alpha=0.6)
axs[0,0].scatter(mag_list, MMI_res_pga_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((MMI_res_pga_nga,MMI_res_pga_bchydro)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[0,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o')
axs[0,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[0,0].axhline(0.73,lw=0.7,ls='--',c='k',label=r'GMM 1$\sigma$')
axs[0,0].axhline(-0.73,lw=0.7,ls='--',c='k')
axs[0,0].grid('-',lw=0.7,alpha=0.5)
axs[0,0].set_ylim(-4.5,3)
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[0,0].axhline(0,lw=0.7,ls='-',c='k')
axs[0,0].text(0.02, 0.98, 'PGA', ha='left', va='top', transform=axs[0,0].transAxes)
axs[0,0].text(-0.1, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs[0,0].transAxes)

# MMI_pga vs Distance
axs[0,1].scatter(rrup_list, MMI_res_pga_bchydro, label='BCHydro', s=10, marker='^', facecolors='none', edgecolors='goldenrod', lw=0.3, alpha=0.6)
axs[0,1].scatter(rrup_list, MMI_res_pga_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((rrup_list,rrup_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((rrup_list,rrup_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[0,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[0,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[0,1].axhline(0.73,lw=0.7,ls='--',c='k')
axs[0,1].axhline(-0.73,lw=0.7,ls='--',c='k')
axs[0,1].set_xscale('log')
axs[0,1].grid('-',lw=0.7,alpha=0.5)
axs[0,1].set_ylim(-4.5,3)
axs[0,1].axhline(0,lw=0.7,ls='-',c='k')
axs[0,1].text(0.02, 0.98, 'PGA', ha='left', va='top', transform=axs[0,1].transAxes)
axs[0,1].text(-0.1, 1, r'$\bf{(b)}$', ha='right', va='top', transform=axs[0,1].transAxes)

# MMI_pgv vs Magnitude
axs[1,0].scatter(mag_list, MMI_res_pgv_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(mag_list,MMI_res_pgv_nga,statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[1,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[1,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[1,0].axhline(0.65,lw=0.7,ls='--',c='k')
axs[1,0].axhline(-0.65,lw=0.7,ls='--',c='k')
axs[1,0].set_xlabel('Magnitude')
axs[1,0].grid('-',lw=0.7,alpha=0.5)
axs[1,0].set_ylim(-4.5,3)
axs[1,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[1,0].axhline(0,lw=0.7,ls='-',c='k')
axs[1,0].text(0.02, 0.98, 'PGV', ha='left', va='top', transform=axs[1,0].transAxes)
axs[1,0].text(-0.1, 1, r'$\bf{(c)}$', ha='right', va='top', transform=axs[1,0].transAxes)

# MMI_pgv vs Distance
axs[1,1].scatter(rrup_list, MMI_res_pgv_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(rrup_list, MMI_res_pgv_nga, statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((rrup_list,rrup_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[1,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[1,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[1,1].axhline(0.65,lw=0.7,ls='--',c='k')
axs[1,1].axhline(-0.65,lw=0.7,ls='--',c='k')
axs[1,1].set_xlabel(r'$R_{rup}$ $(km)$')
axs[1,1].set_xscale('log')
axs[1,1].grid('-',lw=0.7,alpha=0.5)
axs[1,1].set_ylim(-4.5,3)
axs[1,1].axhline(0,lw=0.7,ls='-',c='k')
axs[1,1].set_xlabel(r'$R_{hyp}$ $(km)$')
axs[1,1].text(0.02, 0.98, 'PGV', ha='left', va='top', transform=axs[1,1].transAxes)
axs[1,1].text(-0.1, 1, r'$\bf{(d)}$', ha='right', va='top', transform=axs[1,1].transAxes)

# Set y-label
fig.text(0.01,0.55,r'$\delta_{MMI}$ ${\left(pred - sim\right)}$',va='center',rotation='vertical',fontsize=11)

# Make legend symbols opaque
handles, labels = axs[0,0].get_legend_handles_labels()
for lh in handles: 
    try:
        lh.set_alpha(1)
    except:
        continue

# Set legend
axs[1,0].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.04,-0.405),facecolor='white',
                frameon=True,fontsize=10,title_fontsize=10,ncol=4,markerscale=2)

# Adjust suplot spacing
plt.subplots_adjust(right=0.98,left=0.115,top=0.95,bottom=0.25,hspace=0.3,wspace=0.25)

# Save figure
plt.savefig('/Users/tnye/ONC/manuscript/figures/Fig6_MMI_residuals.png',dpi=300)


#%% PGD Residuals

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
fig, axs = plt.subplots(4,2,figsize=(6.5,6.5))

# PGD vs Magnitude, no noise
axs[0,0].scatter(pgd_mag_list, pgd_res_gold_obs_nonoise, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[0,0].scatter(pgd_mag_list, pgd_res_gold_joint_nonoise, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_nonoise)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_nonoise)),statistic='std', bins=custom_bin_edges)
axs[0,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o')
axs[0,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[0,0].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k',label=r'GMM 1$\sigma$')
axs[0,0].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[0,0].grid('-',lw=0.7,alpha=0.5)
axs[0,0].axhline(0,lw=0.7,ls='-',c='k')
axs[0,0].set_ylim(-3,3)
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[0,0].text(-0.1, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs[0,0].transAxes)
axs[0,0].text(0.02, 0.98, 'No noise', ha='left', va='top', transform=axs[0,0].transAxes)

# PGD vs Distance, no noise
axs[0,1].scatter(pgd_rhyp_list, pgd_res_gold_obs_nonoise, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[0,1].scatter(pgd_rhyp_list, pgd_res_gold_joint_nonoise, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_nonoise)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_nonoise)),statistic='std', bins=custom_bin_edges)
axs[0,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[0,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[0,1].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[0,1].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[0,1].set_ylim(-3,3)
axs[0,1].set_xscale('log')
axs[0,1].grid('-',lw=0.7,alpha=0.5)
axs[0,1].axhline(0,lw=0.7,ls='-',c='k')
axs[0,1].text(-0.1, 1, r'$\bf{(b)}$', ha='right', va='top', transform=axs[0,1].transAxes)
axs[0,1].text(0.02, 0.98, 'No noise', ha='left', va='top', transform=axs[0,1].transAxes)

# PGD vs Magnitude, 10th perc noise
axs[1,0].scatter(pgd_mag_list, pgd_res_gold_obs_10p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[1,0].scatter(pgd_mag_list, pgd_res_gold_joint_10p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_10p,pgd_res_gold_joint_10p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_10p,pgd_res_gold_joint_10p)),statistic='std', bins=custom_bin_edges)
axs[1,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o')
axs[1,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[1,0].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[1,0].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[1,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[1,0].grid('-',lw=0.7,alpha=0.5)
axs[1,0].axhline(0,lw=0.7,ls='-',c='k')
axs[1,0].text(-0.1, 1, r'$\bf{(c)}$', ha='right', va='top', transform=axs[1,0].transAxes)
axs[1,0].text(0.02, 0.98, '10th perc', ha='left', va='top', transform=axs[1,0].transAxes)

# PGD vs Distance, 10th perc noise
axs[1,1].scatter(pgd_rhyp_list, pgd_res_gold_obs_10p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[1,1].scatter(pgd_rhyp_list, pgd_res_gold_joint_10p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_10p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_10p)),statistic='std', bins=custom_bin_edges)
axs[1,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[1,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[1,1].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[1,1].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[1,1].set_xscale('log')
axs[1,1].grid('-',lw=0.7,alpha=0.5)
axs[1,1].axhline(0,lw=0.7,ls='-',c='k')
axs[1,1].text(-0.1, 1, r'$\bf{(d)}$', ha='right', va='top', transform=axs[1,1].transAxes)
axs[1,1].text(0.02, 0.98, '10th perc', ha='left', va='top', transform=axs[1,1].transAxes)

# PGD vs Magnitude, 50th perc noise
axs[2,0].scatter(pgd_mag_list, pgd_res_gold_obs_50p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[2,0].scatter(pgd_mag_list, pgd_res_gold_joint_50p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_50p,pgd_res_gold_joint_50p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_50p,pgd_res_gold_joint_50p)),statistic='std', bins=custom_bin_edges)
axs[2,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o')
axs[2,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[2,0].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[2,0].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[2,0].yaxis.set_label_coords(-.225, .5)
axs[2,0].grid('-',lw=0.7,alpha=0.5)
axs[2,0].axhline(0,lw=0.7,ls='-',c='k')
axs[2,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[2,0].text(-0.1, 1, r'$\bf{(e)}$', ha='right', va='top', transform=axs[2,0].transAxes)
axs[2,0].text(0.02, 0.98, '50th perc', ha='left', va='top', transform=axs[2,0].transAxes)

# PGD vs Distance, 50th perc noise
axs[2,1].scatter(pgd_rhyp_list, pgd_res_gold_obs_50p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[2,1].scatter(pgd_rhyp_list, pgd_res_gold_joint_50p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_50p,pgd_res_gold_joint_50p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_50p,pgd_res_gold_joint_50p)),statistic='std', bins=custom_bin_edges)
axs[2,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[2,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[2,1].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[2,1].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[2,1].set_xscale('log')
axs[2,1].grid('-',lw=0.7,alpha=0.5)
axs[2,1].axhline(0,lw=0.7,ls='-',c='k')
axs[2,1].text(-0.1, 1, r'$\bf{(f)}$', ha='right', va='top', transform=axs[2,1].transAxes)
axs[2,1].text(0.02, 0.98, '50th perc', ha='left', va='top', transform=axs[2,1].transAxes)

# PGD vs Magnitude, 90th perc noise
axs[3,0].scatter(pgd_mag_list, pgd_res_gold_obs_90p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[3,0].scatter(pgd_mag_list, pgd_res_gold_joint_90p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_90p,pgd_res_gold_joint_90p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_90p,pgd_res_gold_joint_90p)),statistic='std', bins=custom_bin_edges)
axs[3,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[3,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[3,0].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[3,0].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[3,0].set_xlabel('Magnitude',fontsize=11)
axs[3,0].grid('-',lw=0.7,alpha=0.5)
axs[3,0].axhline(0,lw=0.7,ls='-',c='k')
axs[3,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[3,0].text(-0.1, 1, r'$\bf{(g)}$', ha='right', va='top', transform=axs[3,0].transAxes)
axs[3,0].text(0.02, 0.98, '90th perc', ha='left', va='top', transform=axs[3,0].transAxes)

# PGD vs Distance, 90th perc noise
axs[3,1].scatter(pgd_rhyp_list, pgd_res_gold_obs_90p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[3,1].scatter(pgd_rhyp_list, pgd_res_gold_joint_90p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_90p,pgd_res_gold_joint_90p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_90p,pgd_res_gold_joint_90p)),statistic='std', bins=custom_bin_edges)
axs[3,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[3,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[3,1].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[3,1].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[3,1].set_xlabel(r'$R_{hyp}$ $(km)$',fontsize=11)
axs[3,1].set_xscale('log')
axs[3,1].grid('-',lw=0.7,alpha=0.5)
axs[3,1].axhline(0,lw=0.7,ls='-',c='k')
axs[3,1].text(-0.1, 1, r'$\bf{(h)}$', ha='right', va='top', transform=axs[3,1].transAxes)
axs[3,1].text(0.02, 0.98, '90th perc', ha='left', va='top', transform=axs[3,1].transAxes)

# Set y label
fig.text(0.01,0.57,r'$\delta_{PGD}$ $ln\left(\frac{pred}{sim}\right)$',va='center',rotation='vertical',fontsize=11)

# Make legend symbols opaque
handles, labels = axs[0,0].get_legend_handles_labels()
for lh in handles: 
    try:
        lh.set_alpha(1)
    except:
        continue

# Set legend
axs[3,0].legend(handles, labels, loc='upper left', bbox_to_anchor=(0,-0.405),
                facecolor='white',frameon=True,fontsize=10,title_fontsize=10,ncol=4,markerscale=2)

# Adjust spacing of subplots
plt.subplots_adjust(right=0.98,left=0.125,top=0.985,bottom=0.14,hspace=0.3,wspace=0.25)

# Save figure
plt.savefig('/Users/tnye/ONC/manuscript/figures/Fig7_PGD_GMM_residuals.png',dpi=300)

