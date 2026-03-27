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
from matplotlib.ticker import MultipleLocator, LogLocator, NullFormatter
import matplotlib.colors as mcolors
from matplotlib import cm

home = '/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW'

# Name of batch for simulations
batch = 'cascadia'

# Description of station types used for naming of IM files
station_types = ['onc-onshore', 'onc-offshore','pnsn']

# Read in list of ruptures
ruptures = np.genfromtxt(f'{home}/simulations/{batch}/data/ruptures.list',dtype=str)

# Initialize lists
pga_res_nga = np.array([])
pgv_res_nga = np.array([])
MMI_res_pga_nga = np.array([])
MMI_res_pgv_nga = np.array([])
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
melgar_sd_list = np.array([])
gold_sd_obs_list = np.array([])
gold_sd_joint_list = np.array([])

# Loop over ruptures
for rupture in ruptures:
    
    # Read in .log file and get magnitude
    run = rupture.replace('.rupt','')
    # logfile = f"{home}/simulations/{batch}/output/ruptures/{rupture.replace('.rupt','.log')}"
    # f = open(logfile, 'r')
    # lines = f.readlines()
    # mag = float(lines[15].split(' ')[-1].split('\n')[0])
    
    # Loop over station types    
    for stations in station_types:
        
        # Read in intensity measure file
        # file = f'{home}/groundmotion_analysis_NGASub-global/IMs/{batch}/{run}_{stations}.csv'
        file = f'{home}/groundmotion_analysis_NGASub-Cascadia/IMs/{batch}/{run}_{stations}.csv'
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
        nga_pgv_sd_list = np.append(nga_pgv_sd_list, np.array(df['NGA_Sub_sd_cm/s']))
        
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
max_pga_sd = np.max(nga_pga_sd_list)
max_pgv_sd = np.max(nga_pgv_sd_list)        
max_pgd_sd = np.max(gold_sd_obs_list)


#%%

## Open quake stuff:
from openquake.hazardlib import imt, const
from openquake.hazardlib.contexts import RuptureContext
from openquake.hazardlib.contexts import DistancesContext
from openquake.hazardlib.contexts import SitesContext
from openquake.hazardlib.gsim.parker_2020 import ParkerEtAl2020SInter
import pandas as pd
import numpy as np

pga_sd_list = []
pgv_sd_list = []

rrup_model = np.arange(0,1201,10)

## Set GMPEs:
NGASub = ParkerEtAl2020SInter(region='Cascadia')

## Set the empty arrays:
median_NGASub = np.array([])
sd_NGASub = np.array([])


## Make contexts:
rctx = RuptureContext()
dctx = DistancesContext()
sctx = SitesContext()

## Magnitude in rupture context is always the same, so grab first one:
rctx.mag = 8.0

## Site context - seems to now need to be from a "site collection", which seems to be a pandas dataframe.
##   Boore 2020 needs vs30, and nothing else - so set them for Nans. These assume the function is used on a single
###    value to predict - if you're doing it in array form (many records for a single earthquake), then make these arrays.
sitecol_dict = {'sids':[1],'vs30':[760],'vs30measured':[np.nan],'z1pt0':[np.nan],'z2pt5':[np.nan]}
sitecollection = pd.DataFrame(sitecol_dict)

## Then put into a sites context:
sctx = SitesContext(sitecol=sitecollection)   

for rrup in rrup_model: 

    ## Use  dataframe distance and vs30 in dist and site context as arrays:
    dctx.rrup = [rrup]
    
    
    _,pga_sd = NGASub.get_mean_and_stddevs(sctx, rctx, dctx, imt.PGA(), [const.StdDev.TOTAL])
    _,pgv_sd = NGASub.get_mean_and_stddevs(sctx, rctx, dctx, imt.PGV(), [const.StdDev.TOTAL])
    
    pga_sd_list.append(pga_sd[0][0])
    pgv_sd_list.append(pgv_sd[0][0])



#%% PGA and PGV Residuals

bres = 13

x_list = [mag_list,rrup_list,mag_list,rrup_list]
y_list = [pga_res_nga,pga_res_nga,pgv_res_nga,pgv_res_nga]
x_edges = [np.linspace(6.5,9.5,bres), np.logspace(np.log10(8),np.log10(1250),bres),
           np.linspace(6.5,9.5,bres),np.logspace(np.log10(8),np.log10(1250),bres)]
y_edges = [np.linspace(-3.5, 2.5, bres)]*4


all_H = []  # Store 2D histogram counts
for i in range(4):
    H, xedges, yedges = np.histogram2d(x_list[i], y_list[i], bins=[x_edges[i], y_edges[i]])
    all_H.append(H)
vmax = max([H.max() for H in all_H])

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'
fig, axs = plt.subplots(2,2,figsize=(3.5,4))

# PGA vs Magnitude
H, xedges, yedges = np.histogram2d(x_list[0], y_list[0], bins=[x_edges[0], y_edges[0]])
X, Y = np.meshgrid(xedges, yedges)
pcm = axs[0,0].pcolormesh(X, Y, H.T, cmap='Blues', norm=plt.matplotlib.colors.LogNorm(vmin=1, vmax=vmax))
bin_centers = (x_edges[0][:-1] + x_edges[0][1:]) / 2
bin_means, _, _ = stats.binned_statistic(x_list[0], y_list[0], statistic='mean', bins=x_edges[0])
bin_std, _, _ = stats.binned_statistic(x_list[0], y_list[0], statistic='std', bins=x_edges[0])
axs[0,0].errorbar(bin_centers, bin_means, yerr=bin_std, fmt='o', markerfacecolor='none', c='orange', markersize=3, elinewidth=1, markeredgewidth=1, capsize=2)
axs[0,0].set_ylim(-3.5, 2.5)
axs[0,0].axhline(0, lw=1, ls='-', c='k')
axs[0,0].axhline(max_pga_sd, lw=1, ls='--', c='k')
axs[0,0].axhline(-max_pga_sd, lw=1, ls='--', c='k')
axs[0,0].text(-0.175, 1, r'(a)', ha='right', va='top', transform=axs[0,0].transAxes)
axs[0,0].set_ylabel(r'$\delta_{PGA}$ $ln\left(\frac{GMM}{sim}\right)$',fontsize=10)
axs[0,0].set_xticklabels([])
axs[0,0].xaxis.set_major_locator(MultipleLocator(1))
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[0,0].yaxis.set_major_locator(MultipleLocator(2))
axs[0,0].yaxis.set_minor_locator(MultipleLocator(1))

# PGA vs Distance
H, xedges, yedges = np.histogram2d(x_list[1], y_list[1], bins=[x_edges[1], y_edges[1]])
X, Y = np.meshgrid(xedges, yedges)
pcm = axs[0,1].pcolormesh(X, Y, H.T, cmap='Blues', norm=plt.matplotlib.colors.LogNorm(vmin=1, vmax=vmax))
bin_centers = (x_edges[1][:-1] + x_edges[1][1:]) / 2
bin_means, _, _ = stats.binned_statistic(x_list[1], y_list[1], statistic='mean', bins=x_edges[1])
bin_std, _, _ = stats.binned_statistic(x_list[1], y_list[1], statistic='std', bins=x_edges[1])
axs[0,1].errorbar(bin_centers, bin_means, yerr=bin_std, fmt='o', markerfacecolor='none', c='orange', markersize=3, elinewidth=1, markeredgewidth=1, capsize=2)
axs[0,1].set_xscale('log')
axs[0,1].set_ylim(-3.5, 2.5)
axs[0,1].axhline(0, lw=1, ls='-', c='k')
axs[0,1].plot(rrup_model, pga_sd_list, lw=1, ls='--', c='k')
axs[0,1].plot(rrup_model, -1*np.array(pga_sd_list), lw=1, ls='--', c='k')
axs[0,1].set_yticklabels([])
axs[0,1].set_xticklabels([])
axs[0,1].yaxis.set_major_locator(MultipleLocator(2))
axs[0,1].yaxis.set_minor_locator(MultipleLocator(1))

# PGV vs Magnitude
H, xedges, yedges = np.histogram2d(x_list[2], y_list[2], bins=[x_edges[2], y_edges[2]])
X, Y = np.meshgrid(xedges, yedges)
pcm = axs[1,0].pcolormesh(X, Y, H.T, cmap='Blues', norm=plt.matplotlib.colors.LogNorm(vmin=1, vmax=vmax))
bin_centers = (x_edges[2][:-1] + x_edges[2][1:]) / 2
bin_means, _, _ = stats.binned_statistic(x_list[2], y_list[2], statistic='mean', bins=x_edges[2])
bin_std, _, _ = stats.binned_statistic(x_list[2], y_list[2], statistic='std', bins=x_edges[2])
axs[1,0].errorbar(bin_centers, bin_means, yerr=bin_std, fmt='o', markerfacecolor='none', c='orange',
                  markersize=3, elinewidth=1, markeredgewidth=1, capsize=2,
                  label='Bin means and\nstd dev')
axs[1,0].set_ylim(-3.5, 2.5)
axs[1,0].axhline(0, lw=1, ls='-', c='k')
axs[1,0].axhline(max_pgv_sd, lw=1, ls='--', c='k')
axs[1,0].axhline(-max_pgv_sd, lw=1, ls='--', c='k')
axs[1,0].text(-0.175, 1, r'(b)', ha='right', va='top', transform=axs[1,0].transAxes)
axs[1,0].set_ylabel(r'$\delta_{PGV}$ $ln\left(\frac{GMM}{sim}\right)$',fontsize=10)
axs[1,0].set_xlabel('Magnitude')
axs[1,0].xaxis.set_major_locator(MultipleLocator(1))
axs[1,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[1,0].yaxis.set_major_locator(MultipleLocator(2))
axs[1,0].yaxis.set_minor_locator(MultipleLocator(1))

# PGV vs Distance
H, xedges, yedges = np.histogram2d(x_list[3], y_list[3], bins=[x_edges[3], y_edges[3]])
X, Y = np.meshgrid(xedges, yedges)
pcm = axs[1,1].pcolormesh(X, Y, H.T, cmap='Blues', norm=plt.matplotlib.colors.LogNorm(vmin=1, vmax=vmax))
bin_centers = (x_edges[3][:-1] + x_edges[3][1:]) / 2
bin_means, _, _ = stats.binned_statistic(x_list[3], y_list[3], statistic='mean', bins=x_edges[3])
bin_std, _, _ = stats.binned_statistic(x_list[3], y_list[3], statistic='std', bins=x_edges[3])
axs[1,1].errorbar(bin_centers, bin_means, yerr=bin_std, fmt='o', markerfacecolor='none',
                  c='orange', markersize=3, elinewidth=1, markeredgewidth=1, capsize=2)
axs[1,1].set_xscale('log')
axs[1,1].set_ylim(-3.5, 2.5)
axs[1,1].axhline(0, lw=1, ls='-', c='k')
axs[1,1].plot(rrup_model, pgv_sd_list, lw=1, ls='--', c='k')
axs[1,1].plot(rrup_model, -1*np.array(pgv_sd_list), lw=1, ls='--', c='k')
axs[1,1].set_xlabel(r'$R_{rup}$ (km)')
axs[1,1].set_yticklabels([])
axs[1,1].yaxis.set_major_locator(MultipleLocator(2))
axs[1,1].yaxis.set_minor_locator(MultipleLocator(1))

# Create a ScalarMappable for the colorbar
norm = plt.matplotlib.colors.LogNorm(vmin=1, vmax=vmax)
cmap = cm.get_cmap('Blues')
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.58, 0.11, 0.385, 0.03])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal', ticks=LogLocator(subs=range(10)))
cbar.ax.xaxis.set_major_locator(LogLocator(base=10))
cbar.ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(2,3,4,5,6,7,8,9)))
cbar.ax.tick_params(which='major', bottom=True, top=False, length=4, labelsize=8)
cbar.ax.tick_params(which='minor', bottom=True, top=False, length=2)
cbar.set_label(r'Record count per bin', fontsize=8)

# Set legend
axs[1,0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.45),facecolor='white',
                frameon=True,fontsize=8,ncol=1,markerscale=1)

# Adjust spacing of subplots
plt.subplots_adjust(right=0.96,left=0.17,top=0.975,bottom=0.3,hspace=0.1,wspace=0.1)

# Save Figure
plt.savefig(f'{home}/manuscript/figures/main-text/Fig6_PGA-PGV_GMM_residuals_cascadia.png',dpi=500)


#%% MMI Residuals

# Define bin resolution
bres = 13

# Set up bin edges
x_list = [mag_list, rrup_list, mag_list, rrup_list]
y_list = [MMI_res_pga_nga, MMI_res_pga_nga, MMI_res_pgv_nga, MMI_res_pgv_nga]
x_edges = [np.linspace(6.5, 9.5, bres), np.logspace(np.log10(8), np.log10(1250), bres),
           np.linspace(6.5, 9.5, bres), np.logspace(np.log10(8), np.log10(1250), bres)]
y_edges = [np.linspace(-3.5, 3, bres)] * 4

# Compute max value across histograms
all_H = [np.histogram2d(x_list[i], y_list[i], bins=[x_edges[i], y_edges[i]])[0] for i in range(4)]
vmax = max([H.max() for H in all_H])

# Initialize figure
fig, axs = plt.subplots(2, 2, figsize=(3.5, 4))

# Labels and annotations
panel_labels = ['(a)', '', '(b)', '']
comp_labels = ['PGA', 'PGA', 'PGV', 'PGV']
ylabels = [r'$\delta_{MMI}$ (GMM–sim)', '',
           r'$\delta_{MMI}$ (GMM–sim)', '']
xlabels = ['', '', 'Magnitude', r'$R_{rup}$ (km)']
xscales = ['linear', 'log', 'linear', 'log']
ref_lines = [0.73, 0.73, 0.65, 0.65]

# Plot
for i in range(4):
    row, col = divmod(i, 2)
    ax = axs[row, col]
        
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
        

    # Histogram
    H, xedges, yedges = np.histogram2d(x_list[i], y_list[i], bins=[x_edges[i], y_edges[i]])
    X, Y = np.meshgrid(xedges, yedges)
    pcm = ax.pcolormesh(X, Y, H.T, cmap='Blues',
                        norm=mpl.colors.LogNorm(vmin=1, vmax=vmax))

    # Bin stats
    bin_centers = (x_edges[i][:-1] + x_edges[i][1:]) / 2
    bin_means, _, _ = stats.binned_statistic(x_list[i], y_list[i], statistic='mean', bins=x_edges[i])
    bin_std, _, _ = stats.binned_statistic(x_list[i], y_list[i], statistic='std', bins=x_edges[i])
    ax.errorbar(bin_centers, bin_means, yerr=bin_std,
                fmt='o', markerfacecolor='none', c='orange',
                markersize=3, elinewidth=1, markeredgewidth=1, capsize=2,
                label='Bin means and\nstd dev' if i == 2 else None)

    # Aesthetics
    ax.set_ylim(-3.5, 3)
    ax.set_xscale(xscales[i])
    ax.axhline(0, lw=1, ls='-', c='k')
    ax.axhline(ref_lines[i], lw=1, ls='--', c='k')
    ax.axhline(-ref_lines[i], lw=1, ls='--', c='k')
    ax.text(0.02, 0.98, comp_labels[i], ha='left', va='top', transform=ax.transAxes, size=9)
    ax.text(-0.175, 1, panel_labels[i], ha='right', va='top', transform=ax.transAxes)
    if ylabels[i]:
        ax.set_ylabel(ylabels[i], labelpad=10, fontsize=10)
    if xlabels[i]:
        ax.set_xlabel(xlabels[i])
    
    if i == 0:
        ax.set_xticklabels([])
    elif i == 1:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    elif i == 3:
        ax.set_yticklabels([])

axs[0,0].xaxis.set_major_locator(MultipleLocator(1))
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[1,0].xaxis.set_major_locator(MultipleLocator(1))
axs[1,0].xaxis.set_minor_locator(MultipleLocator(0.5))

# Create a ScalarMappable for the colorbar
norm = plt.matplotlib.colors.LogNorm(vmin=1, vmax=vmax)
cmap = cm.get_cmap('Blues')
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.58, 0.11, 0.385, 0.03])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal', ticks=LogLocator(subs=range(10)))
cbar.ax.xaxis.set_major_locator(LogLocator(base=10))
cbar.ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(2,3,4,5,6,7,8,9)))
cbar.ax.tick_params(which='major', bottom=True, top=False, length=4, labelsize=8)
cbar.ax.tick_params(which='minor', bottom=True, top=False, length=2)
cbar.set_label(r'Record count per bin', fontsize=8)


# Legend
axs[1, 0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.45), facecolor='white',
                 frameon=True, fontsize=8, ncol=1, markerscale=1)

# Layout
plt.subplots_adjust(right=0.96,left=0.17,top=0.975,bottom=0.3,hspace=0.1,wspace=0.1)


# Save figure
plt.savefig(f'{home}/manuscript/figures/main-text/Fig7_MMI_residuals_cascadia.png',dpi=500)


#%% PGD Residuals

# Define bin resolution
bres = 13
# bres = 7 

# Set up bin edges
x_list = [pgd_mag_list, pgd_rhyp_list, pgd_mag_list, pgd_rhyp_list, pgd_mag_list, pgd_rhyp_list, pgd_mag_list, pgd_rhyp_list]
y_list = [pgd_res_gold_obs_nonoise, pgd_res_gold_obs_nonoise, pgd_res_gold_obs_10p, pgd_res_gold_obs_10p,
          pgd_res_gold_obs_50p, pgd_res_gold_obs_50p, pgd_res_gold_obs_90p, pgd_res_gold_obs_90p]
x_edges = [np.linspace(6.5, 9.5, bres), np.logspace(np.log10(20), np.log10(1250), bres),
           np.linspace(6.5, 9.5, bres), np.logspace(np.log10(20), np.log10(1250), bres),
           np.linspace(6.5, 9.5, bres), np.logspace(np.log10(20), np.log10(1250), bres),
           np.linspace(6.5, 9.5, bres), np.logspace(np.log10(20), np.log10(1250), bres)]
y_edges = [np.linspace(-4.5, 3, bres)] * 8

# Compute max value across histograms
all_H = [np.histogram2d(x_list[i], y_list[i], bins=[x_edges[i], y_edges[i]])[0] for i in range(8)]
vmax = max([H.max() for H in all_H])

fig, axs = plt.subplots(4,2,figsize=(3.5,6.6))

# Labels and annotations
panel_labels = ['(a)', '', '(b)', '', '(c)', '', '(d)', '']
comp_labels = ['No noise', 'No noise', '10th perc', '10th perc', '50th perc',
               '50th perc', '90th perc', '90th perc']
xlabels = ['', '', '', '', '', '', 'Magnitude', r'$R_{hyp}$ (km)']
xscales = ['linear', 'log', 'linear', 'log', 'linear', 'log', 'linear', 'log']

# Plot
for i in range(8):
    row, col = divmod(i, 2)
    ax = axs[row, col]
        
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(1))

    # Histogram
    H, xedges, yedges = np.histogram2d(x_list[i], y_list[i], bins=[x_edges[i], y_edges[i]])
    X, Y = np.meshgrid(xedges, yedges)
    pcm = ax.pcolormesh(X, Y, H.T, cmap='Blues',
                        norm=mpl.colors.LogNorm(vmin=1, vmax=vmax))

    # Bin stats
    bin_centers = (x_edges[i][:-1] + x_edges[i][1:]) / 2
    bin_means, _, _ = stats.binned_statistic(x_list[i], y_list[i], statistic='mean', bins=x_edges[i])
    bin_std, _, _ = stats.binned_statistic(x_list[i], y_list[i], statistic='std', bins=x_edges[i])
    ax.errorbar(bin_centers, bin_means, yerr=bin_std,
                fmt='o', markerfacecolor='none', c='orange',
                markersize=3, elinewidth=1, markeredgewidth=1, capsize=2,
                label='Bin means and\nstd dev' if i == 6 else None)

    # Aesthetics
    ax.set_ylim(-4.5, 3)
    ax.set_xscale(xscales[i])
    ax.axhline(0, lw=1, ls='-', c='k')
    ax.axhline(max_pgd_sd, lw=1, ls='--', c='k')
    ax.axhline(-max_pgd_sd, lw=1, ls='--', c='k')
    ax.text(0.02, 0.985, comp_labels[i], ha='left', va='top', transform=ax.transAxes, size=8)

    if xlabels[i]:
        ax.set_xlabel(xlabels[i])
        
    if col == 0:
        ax.text(-0.225, 1, panel_labels[i], ha='right', va='top', transform=ax.transAxes)
    else:
        ax.set_xlim(xmin=20)
    
    if row != 3:
        ax.set_xticklabels([])
    if col == 1:
        ax.set_yticklabels([])
    else:
        # Set y label
        ax.set_ylabel(r'$\delta_{PGD}$ $ln\left(\frac{GMM}{sim}\right)$')

axs[0,0].xaxis.set_major_locator(MultipleLocator(1))
axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[1,0].xaxis.set_major_locator(MultipleLocator(1))
axs[1,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[2,0].xaxis.set_major_locator(MultipleLocator(1))
axs[2,0].xaxis.set_minor_locator(MultipleLocator(0.5))
axs[3,0].xaxis.set_major_locator(MultipleLocator(1))
axs[3,0].xaxis.set_minor_locator(MultipleLocator(0.5))

# Create a ScalarMappable for the colorbar
norm = plt.matplotlib.colors.LogNorm(vmin=1, vmax=vmax)
cmap = cm.get_cmap('Blues')
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.6, 0.07, 0.375, 0.02])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal', ticks=LogLocator(subs=range(10)))
cbar.ax.xaxis.set_major_locator(LogLocator(base=10))
cbar.ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(2,3,4,5,6,7,8,9)))
cbar.ax.tick_params(which='major', bottom=True, top=False, length=4, labelsize=8)
cbar.ax.tick_params(which='minor', bottom=True, top=False, length=2)
cbar.set_label(r'Record count per bin', fontsize=8)

# Legend
axs[3,0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.47), facecolor='white',
                 frameon=True, fontsize=8, ncol=1, markerscale=1)

# Layout
plt.subplots_adjust(right=0.96, left=0.2, top=0.975, bottom=0.185, hspace=0.1, wspace=0.1)

# Save figure
plt.savefig(f'{home}/manuscript/figures/main-text/Fig5_PGD_GMM_residuals.png',dpi=500)

