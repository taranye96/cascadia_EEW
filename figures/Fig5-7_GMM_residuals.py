#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 13:58:12 2022

@author: tnye
"""

###############################################################################
# Script used to make the IM residual figures. 
###############################################################################

# Imports
from os import path, makedirs
import numpy as np
import pandas as pd
from os import path, makedirs
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

batch = 'cascadia_longer_wfs'
station_types = ['onc-onshore', 'onc-offshore','pnsn']

ruptures = np.genfromtxt(f'/Users/tnye/ONC/simulations/{batch}/data/ruptures.list',dtype=str)

pga_res_nga = np.array([])
pgv_res_nga = np.array([])
MMI_res_pga_nga = np.array([])
MMI_res_pgv_nga = np.array([])
pga_res_bchydro = np.array([])
MMI_res_pga_bchydro = np.array([])
mag_list = np.array([])
rrup_list = np.array([])
pga_list = np.array([])

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

for rupture in ruptures:
    
    run = rupture.replace('.rupt','')
    logfile = f"/Users/tnye/ONC/simulations/{batch}/output/ruptures/{rupture.replace('.rupt','.log')}"
    
    f = open(logfile, 'r')
    lines = f.readlines()
    mag = float(lines[15].split(' ')[-1].split('\n')[0])
        
    for stations in station_types:
        
        file = f'/Users/tnye/ONC/flatfiles/IMs/{batch}/{run}_{stations}.csv'
    
        df = pd.read_csv(file)
        
        rrup_list = np.append(rrup_list,np.array(df['Rrup(km)']))
        mag_list = np.append(mag_list,np.array(df['Mw']))
        
        # NGA-Sub
        pga_res_nga = np.append(pga_res_nga, np.array(df['lnPGA_res_w/noise_NGASub']))
        pgv_res_nga = np.append(pgv_res_nga, np.array(df['lnPGV_res_w/noise_NGASub']))
        MMI_res_pga_nga = np.append(MMI_res_pga_nga, np.array(df['MMI_res_w/noise_NGASub_pga']))
        MMI_res_pgv_nga = np.append(MMI_res_pgv_nga, np.array(df['MMI_res_w/noise_NGASub_pgv']))
        nga_pga_sd_list = np.append(nga_pga_sd_list, np.array(df['NGA_Sub_sd_g']))
        nga_pgv_sd_list = np.append(nga_pga_sd_list, np.array(df['NGA_Sub_sd_cm/s']))
        
        pga_list = np.append(pga_list, np.array(df['syn_PGA_w/noise_g']))
        
        # BCHydro
        pga_res_bchydro = np.append(pga_res_bchydro, np.array(df['lnPGA_res_w/noise_BCHydro']))
        MMI_res_pga_bchydro = np.append(MMI_res_pga_bchydro, np.array(df['MMI_res_w/noise_BCHydro_pga']))
        bchydro_pga_sd_list = np.append(bchydro_pga_sd_list, np.array(df['BCHydro_sd_g']))
        
        stn_list = np.append(stn_list,df['Station'].values)
        rupt_list = np.append(rupt_list,[rupture]*len(df))
        
        # if stations != 'onc-offshore':
            # pgd_res_melgar_noise = np.append(pgd_res_melgar_noise, np.array(df['lnPGD_noise_res_MA15']))
            # pgd_res_gold_noise = np.append(pgd_res_gold_noise, np.array(df['lnPGD_noise_res_GA21']))
            # pgd_res_melgar_nonoise = np.append(pgd_res_melgar_nonoise, np.array(df['lnPGD_res_MA15']))
            # pgd_res_gold_nonoise = np.append(pgd_res_gold_nonoise, np.array(df['lnPGD_res_GA21']))
            # pgd_rhyp_list = np.append(pgd_rhyp_list, np.array(df['Rhyp(km)']))
            # pgd_mag_list = np.append(pgd_mag_list, np.array(df['Mw']))
            # melgar_sd_list = np.append(nga_pga_sd_list, np.array(df['MA15_sd_m']))
            # gold_sd_list = np.append(gold_sd_list, np.array(df['GA21_sd_m']))
        
        if stations == 'onc-onshore':
            gnss_df_10p = pd.read_csv(file.replace(stations,f'{stations}_gnss-10p'))
            gnss_df_50p = pd.read_csv(file.replace(stations,f'{stations}_gnss-50p'))
            gnss_df_90p = pd.read_csv(file.replace(stations,f'{stations}_gnss-90p'))
            # pgd_res_melgar_10p = np.append(pgd_res_melgar_10p, np.array(gnss_df_10p['lnPGD_noise_res_MA15']))
            pgd_res_gold_obs_10p = np.append(pgd_res_gold_obs_10p, np.array(gnss_df_10p['lnPGD_noise_res_GA21-obs']))
            pgd_res_gold_joint_10p = np.append(pgd_res_gold_joint_10p, np.array(gnss_df_10p['lnPGD_noise_res_GA21-joint']))
            # pgd_res_melgar_50p = np.append(pgd_res_melgar_50p, np.array(gnss_df_50p['lnPGD_noise_res_MA15']))
            pgd_res_gold_obs_50p = np.append(pgd_res_gold_obs_50p, np.array(gnss_df_50p['lnPGD_noise_res_GA21-obs']))
            pgd_res_gold_joint_50p = np.append(pgd_res_gold_joint_50p, np.array(gnss_df_50p['lnPGD_noise_res_GA21-joint']))
            # pgd_res_melgar_90p = np.append(pgd_res_melgar_90p, np.array(gnss_df_90p['lnPGD_noise_res_MA15']))
            pgd_res_gold_obs_90p = np.append(pgd_res_gold_obs_90p, np.array(gnss_df_90p['lnPGD_noise_res_GA21-obs']))
            pgd_res_gold_joint_90p = np.append(pgd_res_gold_joint_90p, np.array(gnss_df_90p['lnPGD_noise_res_GA21-joint']))
            # pgd_res_melgar_nonoise = np.append(pgd_res_melgar_nonoise, np.array(gnss_df_50p['lnPGD_res_MA15']))
            pgd_res_gold_obs_nonoise = np.append(pgd_res_gold_obs_nonoise, np.array(gnss_df_50p['lnPGD_res_GA21-obs']))
            pgd_res_gold_joint_nonoise = np.append(pgd_res_gold_joint_nonoise, np.array(gnss_df_50p['lnPGD_res_GA21-joint']))
            pgd_rhyp_list = np.append(pgd_rhyp_list, np.array(gnss_df_50p['Rhyp(km)']))
            pgd_mag_list = np.append(pgd_mag_list, np.array(gnss_df_50p['Mw']))
            # melgar_sd_list = np.append(nga_pga_sd_list, np.array(gnss_df_50p['MA15_sd_m']))
            gold_sd_obs_list = np.append(gold_sd_obs_list, np.array(gnss_df_50p['GA21-obs_sd_m']))
            gold_sd_joint_list = np.append(gold_sd_joint_list, np.array(gnss_df_50p['GA21-joint_sd_m']))

pga_sd = np.max(np.concatenate((nga_pga_sd_list,bchydro_pga_sd_list)))
pgv_sd = np.max(nga_pgv_sd_list)           
# pgd_sd = np.max(np.concatenate((melgar_sd_list,gold_sd_list)))  
pgd_sd = np.max(np.concatenate((gold_sd_obs_list,gold_sd_joint_list)))  


#%% PGD Residuals

figdir = f'/Users/tnye/ONC/manuscript/figures/PGD_GMM_residuals_GA21_updated.png'


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
# fig, axs = plt.subplots(3,2,figsize=(6.5,4.86))
fig, axs = plt.subplots(4,2,figsize=(6.5,6.5))

axs[0,0].scatter(pgd_mag_list, pgd_res_gold_obs_nonoise, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[0,0].scatter(pgd_mag_list, pgd_res_gold_joint_nonoise, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_nonoise)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_nonoise)),statistic='std', bins=custom_bin_edges)
# bin_means, _, _ = stats.binned_statistic(pgd_mag_list,pgd_res_gold_nonoise,statistic='mean', bins=custom_bin_edges)
# bin_std, _, _ = stats.binned_statistic(pgd_mag_list,pgd_res_gold_nonoise,statistic='std', bins=custom_bin_edges)
axs[0,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o')
axs[0,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[0,0].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k',label=r'GMM 1$\sigma$')
axs[0,0].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
# axs[0,0].xaxis.set_ticklabels([])
axs[0,0].grid('-',lw=0.7,alpha=0.5)
axs[0,0].axhline(0,lw=0.7,ls='-',c='k')
axs[0,0].set_ylim(-3,3)
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[0,0].text(-0.1, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs[0,0].transAxes)
axs[0,0].text(0.02, 0.98, 'No noise', ha='left', va='top', transform=axs[0,0].transAxes)

axs[0,1].scatter(pgd_rhyp_list, pgd_res_gold_obs_nonoise, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[0,1].scatter(pgd_rhyp_list, pgd_res_gold_joint_nonoise, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_nonoise)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_nonoise)),statistic='std', bins=custom_bin_edges)
# bin_means, _, _ = stats.binned_statistic(pgd_rhyp_list,pgd_res_gold_nonoise,statistic='mean', bins=custom_bin_edges)
# bin_std, _, _ = stats.binned_statistic(pgd_rhyp_list,pgd_res_gold_nonoise,statistic='std', bins=custom_bin_edges)
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
# axs[0,1].xaxis.set_ticklabels([])

axs[1,0].scatter(pgd_mag_list, pgd_res_gold_obs_10p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[1,0].scatter(pgd_mag_list, pgd_res_gold_joint_10p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_10p,pgd_res_gold_joint_10p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_10p,pgd_res_gold_joint_10p)),statistic='std', bins=custom_bin_edges)
# bin_means, _, _ = stats.binned_statistic(pgd_mag_list,pgd_res_gold_10p,statistic='mean', bins=custom_bin_edges)
# bin_std, _, _ = stats.binned_statistic(pgd_mag_list,pgd_res_gold_10p,statistic='std', bins=custom_bin_edges)
axs[1,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o')
axs[1,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[1,0].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[1,0].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
# axs[1,0].xaxis.set_ticklabels([])
axs[1,0].xaxis.set_minor_locator(MultipleLocator(0.5))
# axs[1,0].set_ylabel('ln(pred / sim)')
axs[1,0].grid('-',lw=0.7,alpha=0.5)
axs[1,0].axhline(0,lw=0.7,ls='-',c='k')
axs[1,0].text(-0.1, 1, r'$\bf{(c)}$', ha='right', va='top', transform=axs[1,0].transAxes)
axs[1,0].text(0.02, 0.98, '10th perc', ha='left', va='top', transform=axs[1,0].transAxes)

axs[1,1].scatter(pgd_rhyp_list, pgd_res_gold_obs_10p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[1,1].scatter(pgd_rhyp_list, pgd_res_gold_joint_10p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_10p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_nonoise,pgd_res_gold_joint_10p)),statistic='std', bins=custom_bin_edges)
# bin_means, _, _ = stats.binned_statistic(pgd_rhyp_list,pgd_res_gold_10p,statistic='mean', bins=custom_bin_edges)
# bin_std, _, _ = stats.binned_statistic(pgd_rhyp_list,pgd_res_gold_10p,statistic='std', bins=custom_bin_edges)
axs[1,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[1,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[1,1].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[1,1].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[1,1].set_xscale('log')
# axs[1,1].xaxis.set_ticklabels([])
axs[1,1].grid('-',lw=0.7,alpha=0.5)
axs[1,1].axhline(0,lw=0.7,ls='-',c='k')
axs[1,1].text(-0.1, 1, r'$\bf{(d)}$', ha='right', va='top', transform=axs[1,1].transAxes)
axs[1,1].text(0.02, 0.98, '10th perc', ha='left', va='top', transform=axs[1,1].transAxes)

axs[2,0].scatter(pgd_mag_list, pgd_res_gold_obs_50p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[2,0].scatter(pgd_mag_list, pgd_res_gold_joint_50p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_50p,pgd_res_gold_joint_50p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_50p,pgd_res_gold_joint_50p)),statistic='std', bins=custom_bin_edges)
# bin_means, _, _ = stats.binned_statistic(pgd_mag_list,pgd_res_gold_50p,statistic='mean', bins=custom_bin_edges)
# bin_std, _, _ = stats.binned_statistic(pgd_mag_list,pgd_res_gold_50p,statistic='std', bins=custom_bin_edges)
axs[2,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o')
axs[2,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[2,0].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[2,0].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[2,0].yaxis.set_label_coords(-.225, .5)
axs[2,0].grid('-',lw=0.7,alpha=0.5)
axs[2,0].axhline(0,lw=0.7,ls='-',c='k')
# axs[2,0].xaxis.set_ticklabels([])
axs[2,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[2,0].text(-0.1, 1, r'$\bf{(e)}$', ha='right', va='top', transform=axs[2,0].transAxes)
axs[2,0].text(0.02, 0.98, '50th perc', ha='left', va='top', transform=axs[2,0].transAxes)
# axs[2,0].set_ylabel(r'$\delta_{PGD}$ ln(pred / sim)',fontsize=12)
# axs[2,0].set_ylabel(r'$\delta_{PGD}$ $ln\left(\frac{pred}{sim}\right)$',fontsize=12)

axs[2,1].scatter(pgd_rhyp_list, pgd_res_gold_obs_50p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[2,1].scatter(pgd_rhyp_list, pgd_res_gold_joint_50p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_50p,pgd_res_gold_joint_50p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_50p,pgd_res_gold_joint_50p)),statistic='std', bins=custom_bin_edges)
# bin_means, _, _ = stats.binned_statistic(pgd_rhyp_list,pgd_res_gold_50p,statistic='mean', bins=custom_bin_edges)
# bin_std, _, _ = stats.binned_statistic(pgd_rhyp_list,pgd_res_gold_50p,statistic='std', bins=custom_bin_edges)
axs[2,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[2,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[2,1].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[2,1].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[2,1].set_xscale('log')
axs[2,1].grid('-',lw=0.7,alpha=0.5)
axs[2,1].axhline(0,lw=0.7,ls='-',c='k')
# axs[2,1].xaxis.set_ticklabels([])
axs[2,1].text(-0.1, 1, r'$\bf{(f)}$', ha='right', va='top', transform=axs[2,1].transAxes)
axs[2,1].text(0.02, 0.98, '50th perc', ha='left', va='top', transform=axs[2,1].transAxes)

axs[3,0].scatter(pgd_mag_list, pgd_res_gold_obs_90p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[3,0].scatter(pgd_mag_list, pgd_res_gold_joint_90p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_90p,pgd_res_gold_joint_90p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_mag_list,pgd_mag_list)),np.concatenate((pgd_res_gold_obs_90p,pgd_res_gold_joint_90p)),statistic='std', bins=custom_bin_edges)
# bin_means, _, _ = stats.binned_statistic(pgd_mag_list,pgd_res_gold_90p,statistic='mean', bins=custom_bin_edges)
# bin_std, _, _ = stats.binned_statistic(pgd_mag_list,pgd_res_gold_90p,statistic='std', bins=custom_bin_edges)
axs[3,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[3,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o", label=r'Bin mean and $1\sigma$')
axs[3,0].axhline(np.log(pgd_sd),lw=0.7,ls='--',c='k')
axs[3,0].axhline(-np.log(pgd_sd),lw=0.7,ls='--',c='k')
# axs[3,0].set_ylabel('ln(pred / sim)')
axs[3,0].set_xlabel('Magnitude',fontsize=11)
axs[3,0].grid('-',lw=0.7,alpha=0.5)
axs[3,0].axhline(0,lw=0.7,ls='-',c='k')
axs[3,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[3,0].text(-0.1, 1, r'$\bf{(g)}$', ha='right', va='top', transform=axs[3,0].transAxes)
axs[3,0].text(0.02, 0.98, '90th perc', ha='left', va='top', transform=axs[3,0].transAxes)

axs[3,1].scatter(pgd_rhyp_list, pgd_res_gold_obs_90p, label=r'GA21$_{obs}$', marker='1', s=10, lw=0.3, alpha=0.6)
axs[3,1].scatter(pgd_rhyp_list, pgd_res_gold_joint_90p, label=r'GA21$_{joint}$', marker='x', s=10, lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_90p,pgd_res_gold_joint_90p)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((pgd_rhyp_list,pgd_rhyp_list)),np.concatenate((pgd_res_gold_obs_90p,pgd_res_gold_joint_90p)),statistic='std', bins=custom_bin_edges)
# bin_means, _, _ = stats.binned_statistic(pgd_rhyp_list,pgd_res_gold_90p,statistic='mean', bins=custom_bin_edges)
# bin_std, _, _ = stats.binned_statistic(pgd_rhyp_list,pgd_res_gold_90p,statistic='std', bins=custom_bin_edges)
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

fig.text(0.01,0.57,r'$\delta_{PGD}$ $ln\left(\frac{pred}{sim}\right)$',va='center',rotation='vertical',fontsize=11)

handles, labels = axs[0,0].get_legend_handles_labels()
for lh in handles: 
    try:
        lh.set_alpha(1)
    except:
        continue

# axs[2,0].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.325,-0.4),facecolor='white',
#                 frameon=True,fontsize=10,title_fontsize=10,ncol=3,markerscale=2)

# axs[3,0].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.05,-0.405),facecolor='white',
#                 frameon=True,fontsize=10,title_fontsize=10,ncol=4,markerscale=2)

axs[3,0].legend(handles, labels, loc='upper left', bbox_to_anchor=(0,-0.405),
                facecolor='white',frameon=True,fontsize=10,title_fontsize=10,ncol=4,markerscale=2)

# plt.subplots_adjust(right=0.98,left=0.125,top=0.985,bottom=0.175,hspace=0.15,wspace=0.25)
plt.subplots_adjust(right=0.98,left=0.125,top=0.985,bottom=0.14,hspace=0.3,wspace=0.25)
plt.savefig(figdir,dpi=300)


#%% High Frequency IM Residuals


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
# fig, axs = plt.subplots(3,2,figsize=(6.5,7.75))
fig, axs = plt.subplots(4,2,figsize=(6.5,6.25))

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
axs[0,0].xaxis.set_ticklabels([])
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
# axs[0,0].tick_params(axis='both', direction='out', length=5, width=2, colors='black')
axs[0,0].text(0.02, 0.98, r'$\delta_{PGA}$', ha='left', va='top', transform=axs[0,0].transAxes)
axs[0,0].text(-0.15, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs[0,0].transAxes)

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
axs[0,1].xaxis.set_ticklabels([])
axs[0,1].set_ylim(-4.5,3)
axs[0,1].axhline(np.log(pga_sd),lw=0.7,ls='--',c='k')
axs[0,1].axhline(-np.log(pga_sd),lw=0.7,ls='--',c='k')
axs[0,1].text(0.02, 0.98, r'$\delta_{PGA}$', ha='left', va='top', transform=axs[0,1].transAxes)
axs[0,1].text(-0.1, 1, r'$\bf{(b)}$', ha='right', va='top', transform=axs[0,1].transAxes)

axs[1,0].scatter(mag_list, pgv_res_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(mag_list,pgv_res_nga,statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[1,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[1,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[1,0].set_ylabel('ln(pred / sim)')
axs[1,0].grid('-',lw=0.7,alpha=0.5)
axs[1,0].set_ylim(-4.5,3)
axs[1,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[1,0].axhline(0,lw=0.7,ls='-',c='k')
axs[1,0].xaxis.set_ticklabels([])
axs[1,0].axhline(np.log(pgv_sd),lw=0.7,ls='--',c='k')
axs[1,0].axhline(-np.log(pgv_sd),lw=0.7,ls='--',c='k')
axs[1,0].text(0.02, 0.98, r'$\delta_{PGV}$', ha='left', va='top', transform=axs[1,0].transAxes)
axs[1,0].text(-0.15, 1, r'$\bf{(c)}$', ha='right', va='top', transform=axs[1,0].transAxes)

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
axs[1,1].xaxis.set_ticklabels([])
axs[1,1].axhline(np.log(pgv_sd),lw=0.7,ls='--',c='k')
axs[1,1].axhline(-np.log(pgv_sd),lw=0.7,ls='--',c='k')
axs[1,1].text(0.02, 0.98, r'$\delta_{PGV}$', ha='left', va='top', transform=axs[1,1].transAxes)
axs[1,1].text(-0.1, 1, r'$\bf{(d)}$', ha='right', va='top', transform=axs[1,1].transAxes)

axs[2,0].scatter(mag_list, MMI_res_pga_bchydro, label='BCHydro', s=10, marker='^', facecolors='none', edgecolors='goldenrod', lw=0.3, alpha=0.6)
axs[2,0].scatter(mag_list, MMI_res_pga_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((MMI_res_pga_nga,MMI_res_pga_bchydro)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[2,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[2,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[2,0].axhline(0.73,lw=0.7,ls='--',c='k')
axs[2,0].axhline(-0.73,lw=0.7,ls='--',c='k')
axs[2,0].set_ylabel('pred - sim')
axs[2,0].grid('-',lw=0.7,alpha=0.5)
axs[2,0].set_ylim(-4.5,3)
axs[2,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[2,0].axhline(0,lw=0.7,ls='-',c='k')
axs[2,0].xaxis.set_ticklabels([])
axs[2,0].text(0.02, 0.98, r'$\delta_{MMI, PGA}$', ha='left', va='top', transform=axs[2,0].transAxes)
axs[2,0].text(-0.15, 1, r'$\bf{(e)}$', ha='right', va='top', transform=axs[2,0].transAxes)

axs[2,1].scatter(rrup_list, MMI_res_pga_bchydro, label='BCHydro', s=10, marker='^', facecolors='none', edgecolors='goldenrod', lw=0.3, alpha=0.6)
axs[2,1].scatter(rrup_list, MMI_res_pga_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(np.concatenate((rrup_list,rrup_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((rrup_list,rrup_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[2,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[2,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[2,1].axhline(0.73,lw=0.7,ls='--',c='k')
axs[2,1].axhline(-0.73,lw=0.7,ls='--',c='k')
axs[2,1].set_xscale('log')
axs[2,1].grid('-',lw=0.7,alpha=0.5)
axs[2,1].set_ylim(-4.5,3)
axs[2,1].axhline(0,lw=0.7,ls='-',c='k')
axs[2,1].xaxis.set_ticklabels([])
axs[2,1].text(0.02, 0.98, r'$\delta_{MMI, PGA}$', ha='left', va='top', transform=axs[2,1].transAxes)
axs[2,1].text(-0.1, 1, r'$\bf{(f)}$', ha='right', va='top', transform=axs[2,1].transAxes)

axs[3,0].scatter(mag_list, MMI_res_pgv_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.arange(6.5,9.75,0.25)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(mag_list,MMI_res_pgv_nga,statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((mag_list,mag_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[3,0].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[3,0].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[3,0].axhline(0.65,lw=0.7,ls='--',c='k')
axs[3,0].axhline(-0.65,lw=0.7,ls='--',c='k')
axs[3,0].set_ylabel('pred - sim')
axs[3,0].set_xlabel('Magnitude')
axs[3,0].grid('-',lw=0.7,alpha=0.5)
axs[3,0].set_ylim(-4.5,3)
axs[3,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[3,0].axhline(0,lw=0.7,ls='-',c='k')
axs[3,0].text(0.02, 0.98, r'$\delta_{MMI, PGV}$', ha='left', va='top', transform=axs[3,0].transAxes)
axs[3,0].text(-0.15, 1, r'$\bf{(g)}$', ha='right', va='top', transform=axs[3,0].transAxes)

axs[3,1].scatter(rrup_list, MMI_res_pgv_nga, label='NGA-Sub', s=10, facecolors='none', edgecolors='mediumpurple', lw=0.3, alpha=0.6)
custom_bin_edges = np.logspace(np.log10(1), np.log10(1300), 13)
bin_centers = (custom_bin_edges[1:] + custom_bin_edges[:-1]) / 2
bin_means, _, _ = stats.binned_statistic(rrup_list, MMI_res_pgv_nga, statistic='mean', bins=custom_bin_edges)
bin_std, _, _ = stats.binned_statistic(np.concatenate((rrup_list,rrup_list)),np.concatenate((pga_res_nga,pga_res_bchydro)),statistic='std', bins=custom_bin_edges)
axs[3,1].scatter(bin_centers, bin_means, s=5, color='k', marker='o', label='Bin mean')
axs[3,1].errorbar(bin_centers, bin_means, yerr=bin_std, markersize=1, elinewidth=0.75, color='k', fmt="o")
axs[3,1].axhline(0.65,lw=0.7,ls='--',c='k')
axs[3,1].axhline(-0.65,lw=0.7,ls='--',c='k')
axs[3,1].set_xlabel(r'$R_{rup}$ $(km)$')
axs[3,1].set_xscale('log')
axs[3,1].grid('-',lw=0.7,alpha=0.5)
axs[3,1].set_ylim(-4.5,3)
axs[3,1].axhline(0,lw=0.7,ls='-',c='k')
axs[3,1].set_xlabel(r'$R_{hyp}$ $(km)$')
axs[3,1].text(0.02, 0.98, r'$\delta_{MMI, PGV}$', ha='left', va='top', transform=axs[3,1].transAxes)
axs[3,1].text(-0.1, 1, r'$\bf{(h)}$', ha='right', va='top', transform=axs[3,1].transAxes)

handles, labels = axs[0,0].get_legend_handles_labels()
for lh in handles: 
    try:
        lh.set_alpha(1)
    except:
        continue

axs[3,0].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.05,-0.405),facecolor='white',
                frameon=True,fontsize=10,title_fontsize=10,ncol=4,markerscale=2)


plt.subplots_adjust(right=0.98,left=0.1,top=0.985,bottom=0.14,hspace=0.15,wspace=0.25)
plt.savefig('/Users/tnye/ONC/manuscript/figures/HF_GMM_residuals.png',dpi=300)


#%% PGA and PGV Residuals

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
# fig, axs = plt.subplots(3,2,figsize=(6.5,7.75))
fig, axs = plt.subplots(2,2,figsize=(6.5,4))

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
# axs[0,0].xaxis.set_ticklabels([])
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
# axs[0,0].tick_params(axis='both', direction='out', length=5, width=2, colors='black')
# axs[0,0].text(0.02, 0.98, r'$\delta_{PGA}$', ha='left', va='top', transform=axs[0,0].transAxes)
axs[0,0].text(-0.15, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs[0,0].transAxes)
# axs[0,0].set_ylabel(r'$\delta_{PGA}$ ln(pred / sim)')
axs[0,0].set_ylabel(r'$\delta_{PGA}$ $ln\left(\frac{pred}{sim}\right)$',fontsize=10)

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
# axs[0,1].xaxis.set_ticklabels([])
axs[0,1].set_ylim(-4.5,3)
axs[0,1].axhline(np.log(pga_sd),lw=0.7,ls='--',c='k')
axs[0,1].axhline(-np.log(pga_sd),lw=0.7,ls='--',c='k')
# axs[0,1].text(0.02, 0.98, r'$\delta_{PGA}$', ha='left', va='top', transform=axs[0,1].transAxes)
axs[0,1].text(-0.1, 1, r'$\bf{(b)}$', ha='right', va='top', transform=axs[0,1].transAxes)

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
# axs[1,0].text(0.02, 0.98, r'$\delta_{PGV}$', ha='left', va='top', transform=axs[1,0].transAxes)
axs[1,0].text(-0.15, 1, r'$\bf{(c)}$', ha='right', va='top', transform=axs[1,0].transAxes)
axs[1,0].set_xlabel('Magnitude')
# axs[1,0].set_ylabel(r'$\delta_{PGV}$ ln(pred / sim)')
axs[1,0].set_ylabel(r'$\delta_{PGV}$ $ln\left(\frac{pred}{sim}\right)$',fontsize=10)

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
# axs[1,1].text(0.02, 0.98, r'$\delta_{PGV}$', ha='left', va='top', transform=axs[1,1].transAxes)
axs[1,1].text(-0.1, 1, r'$\bf{(d)}$', ha='right', va='top', transform=axs[1,1].transAxes)
axs[1,1].set_xlabel(r'$R_{hyp}$ $(km)$')

handles, labels = axs[0,0].get_legend_handles_labels()
for lh in handles: 
    try:
        lh.set_alpha(1)
    except:
        continue

axs[1,0].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.04,-0.405),facecolor='white',
                frameon=True,fontsize=10,title_fontsize=10,ncol=4,markerscale=2)


plt.subplots_adjust(right=0.98,left=0.1,top=0.95,bottom=0.25,hspace=0.3,wspace=0.25)
plt.savefig('/Users/tnye/ONC/manuscript/figures/PGA-PGV_GMM_residuals.png',dpi=300)


#%% MMI Residuals

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
fig, axs = plt.subplots(2,2,figsize=(6.5,4))

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
# axs[0,0].xaxis.set_ticklabels([])
axs[0,0].text(0.02, 0.98, 'PGA', ha='left', va='top', transform=axs[0,0].transAxes)
axs[0,0].text(-0.1, 1, r'$\bf{(a)}$', ha='right', va='top', transform=axs[0,0].transAxes)
# axs[0,0].set_ylabel(r'$\delta_{MMI}$ ${\left(pred - sim\right)}$',fontsize=10)

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
# axs[0,1].xaxis.set_ticklabels([])
axs[0,1].text(0.02, 0.98, 'PGA', ha='left', va='top', transform=axs[0,1].transAxes)
axs[0,1].text(-0.1, 1, r'$\bf{(b)}$', ha='right', va='top', transform=axs[0,1].transAxes)

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
# axs[1,0].set_ylabel(r'$\delta_{MMI}$ ${\left(pred - sim\right)}$',fontsize=10)

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

fig.text(0.01,0.55,r'$\delta_{MMI}$ ${\left(pred - sim\right)}$',va='center',rotation='vertical',fontsize=11)

handles, labels = axs[0,0].get_legend_handles_labels()
for lh in handles: 
    try:
        lh.set_alpha(1)
    except:
        continue

axs[1,0].legend(handles, labels, loc='upper left', bbox_to_anchor=(0.04,-0.405),facecolor='white',
                frameon=True,fontsize=10,title_fontsize=10,ncol=4,markerscale=2)


plt.subplots_adjust(right=0.98,left=0.115,top=0.95,bottom=0.25,hspace=0.3,wspace=0.25)
plt.savefig('/Users/tnye/ONC/manuscript/figures/MMI_residuals.png',dpi=300)


