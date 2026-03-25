#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 10:34:27 2025

@author: tnye
"""

import pygmt
import os
import io

home = "/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects"
out = "cascadia_events_stations"

# File paths
topo_grd = f"{home}/ONC-EEW/GMT/data/DEMs/cascadia30.grd"
event_coords = f"{home}/ONC-EEW/files/events_coords.txt"
pnsn_coords = f"{home}/ONC-EEW/files/fakequakes/station_info/pnsn_coords.txt"
offshore_coords = f"{home}/ONC-EEW/files/fakequakes/station_info/offshore_coords.txt"
onshore_coords = f"{home}/ONC-EEW/files/fakequakes/station_info/onshore_coords.txt"
trench = f"{home}/ONC-EEW/GMT/data/CascadiaTrench.txt"
topo_no = f"{home}/ONC-EEW/GMT/topo/topo_no.grd"
topo_grad = f"{home}/ONC-EEW/GMT/topo/topo.grad"
topo_cpt = f"{home}/ONC-EEW/GMT/cpt/topo.cpt"
bathy_cpt = f"{home}/ONC-EEW/GMT/cpt/bathy.cpt"
ruptures_path=f'{home}/ONC-EEW/GMT/plotting_files'
rupt_xy = f"{ruptures_path}/rupture_xy/cascadia/cascadia-000089.xy"
rupt_epi = f"{ruptures_path}/rupture_epi/cascadia/cascadia.000089_epi.txt"
slip_cpt = f"{home}/ONC-EEW/GMT/cpt/cascadia_000089.cpt"

inc=10
Mw=8.9

region = [-129, -121.25, 39.75, 51]
projection = "M2.25i"

datamap_grid = pygmt.datasets.load_earth_relief(resolution='01m', region=region)
rupture_grid = pygmt.datasets.load_earth_relief(resolution='01m', region=region)

#%%

# # Make cpt files
# slip_cpt = pygmt.makecpt(cmap="/Users/tnye/cpt/magma_white.cpt",series=[0, 20, 0.02],reverse=True,continuous=True,)
# topo_cpt = pygmt.makecpt(cmap="/Users/tnye/cpt/gray_lt.cpt",series=[-0.65, 0.65, 0.01],continuous=True,)
# bathy_cpt = pygmt.makecpt(cmap="/Users/tnye/cpt/gr58_hult.cpt",series=[-5000, 0, 1],continuous=True,reverse=True,)

# Set general config
fig = pygmt.Figure()
pygmt.config(
    FONT_ANNOT_PRIMARY="10p,Helvetica",
    MAP_FRAME_TYPE="plain",
    PS_MEDIA="a2",
    FORMAT_GEO_MAP="ddd",
    FONT_LABEL="16p"
)

with fig.subplot(nrows=1, ncols=1, figsize=("2.75i", "6i"), margins=["0.15i"]):
    
    with fig.set_panel(0):
        
        # Plot base
        fig.basemap(region=region, projection=projection, frame=["WesN", "2p"])
        fig.grdimage(grid=datamap_grid, cmap=bathy_cpt, shading='+d', transparency=20, region=region, projection=projection)
        
        # Coastline and bathymetry
        # fig.coast(land='white', transparency=20, region=region, projection=projection)
        fig.coast(water="197/212/253", land='white', transparency=30, region=region, projection=projection)
        fig.coast(shorelines='0.25p,black', transparency=0, region=region, projection=projection)
        fig.coast(borders=["1/0.75p,black", "2/0.25p,black"], region=region, projection=projection)
        
        # Trench
        fig.plot(data=trench, pen="0.5p,black", style="f0.75c/0.15c+l+t", fill="black", region=region, projection=projection)
        fig.plot(x=[0,0], y=[0,0], pen='0.5,black', region=region, projection=projection, label='Cascadia Trench')
        
        # Events
        fig.plot(data=event_coords, style="a0.2c", pen="0.5p,black", region=region, projection=projection, label='Simulation epicenters')
        
        # Stations
        fig.plot(data=pnsn_coords, style="i0.18c", fill="54/156/110", pen="black", region=region, projection=projection, label='PNSN seismometers')
        fig.plot(data=offshore_coords, style="d0.2c", fill="240/138/36", pen="black", region=region, projection=projection, label='ONC offshore seismometers')
        fig.plot(data=onshore_coords, style="n0.18c", fill="142/83/201", pen="black", region=region, projection=projection, label='ONC collocated GNSS stations')
        fig.plot(data=onshore_coords, style="n0.2c", fill="white", pen="white", transparency=100, region=region, projection=projection, label='and seismometers')
        
        # fig.legend(position="JBL+jBL+o0.2c", box="+gwhite+p1p")
        with pygmt.config(FONT_ANNOT_PRIMARY='8p'):
            fig.legend(position="JBC+jTC+o0c/0.25c", box="+gwhite+p1p", region=region, projection=projection)        

fig.show()


# Save figure
fig.savefig(f'{home}/ONC-EEW/manuscript/figures/subfigs/{out}.png', dpi=500)
