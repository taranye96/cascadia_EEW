#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 10:34:27 2025

@author: tnye
"""

import numpy as np
import pygmt
from glob import glob

home = "/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects"

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

region = [-129, -121.25, 39.75, 51]
projection = "M2.25i"

datamap_grid = pygmt.datasets.load_earth_relief(resolution='01m', region=region)
rupture_grid = pygmt.datasets.load_earth_relief(resolution='01m', region=region)

rupt_files = sorted(glob(f'{home}/ONC-EEW//simulations/cascadia/output/ruptures/*.rupt'))

#%%

for ruptfile in rupt_files:
    
    rupture = ruptfile.split('/')[-1].replace('.rupt','')
    rupt_xy = f"{ruptures_path}/rupture_xy/cascadia/{rupture}.xy"
    rupt_epi = f"{ruptures_path}/rupture_epi/cascadia/{rupture}_epi.txt"
    
    rupt_data = np.genfromtxt(ruptfile)
    ss_slip = rupt_data[:,8]
    ds_slip = rupt_data[:,9]
    total_slip = np.sqrt(ss_slip**2 + ds_slip**2)
    max_slip = np.max(total_slip)
    
    if max_slip < 1:
        inc = 0.25
        scale = 1
    elif max_slip >= 1 and max_slip < 2:
        inc = 0.5
        scale = 2
    elif max_slip >=1 and max_slip < 5:
        inc = 1
        scale = 5
    elif max_slip >=5 and max_slip < 10:
        inc = 2
        scale = 10
    elif max_slip >=10 and max_slip < 20:
        inc = 5
        scale = 20
    elif max_slip >=20 and max_slip < 50:
        inc = 10
        scale = 50
    elif max_slip >= 50:
        inc = 20
        scale = max_slip
    
    logfile = ruptfile.replace('.rupt','.log')
    f = open(logfile, "r")
    lines = f.readlines()
    Mw = round(float(lines[15].split(' ')[-1].strip('\n')),1)
    f.close()

    
    # Make cpt file
    slip_cpt = f"{home}/ONC-EEW/GMT/cpt/{rupture}.cpt"
    pygmt.makecpt(cmap="/Users/tnye/cpt/magma_white.cpt",series=[0, scale, 0.02],
                  reverse=True,continuous=True,output=slip_cpt)

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
            
            # Basemap
            fig.basemap(region=region, projection=projection, frame=["WesN", "2p"])
            fig.grdimage(grid=datamap_grid, cmap=bathy_cpt, shading='+d', transparency=20, region=region, projection=projection)
            
            # Coastline and bathymetry
            # fig.coast(land='white', transparency=20, region=region, projection=projection)
            fig.coast(water="197/212/253", land='white', transparency=30, region=region, projection=projection)
            fig.coast(shorelines='0.25p,black', transparency=0, region=region, projection=projection)
            fig.coast(borders=["1/0.75p,black", "2/0.25p,black"], region=region, projection=projection)
         
            # Plot rupt model
            fig.plot(data=rupt_xy, cmap=slip_cpt, region=region, projection=projection, fill='+z')
            fig.plot(data=rupt_xy, pen='0.01p,180/180/180', region=region, projection=projection)
    
            # Plot epicenter as a filled circle (Sa0.45c in GMT)
            fig.plot(data=rupt_epi, region=region, projection=projection, style="a0.45c", pen="0.5p,black")
            
            # Slip scale bar
            with pygmt.config(FONT_ANNOT_PRIMARY='16p'):
                fig.colorbar(cmap=slip_cpt, frame=[f"xa{inc}+lSlip (m)"], position="JBC+o0c/0.8c+w1.6i/0.15i", projection=projection, region=region)
            
            # Add magnitude label text (at fixed lon,lat with white font)
            fig.text(x=-128.15, y=40.12, text=f"M{Mw}", font="12p,Helvetica,black", justify="CM", region=region, projection=projection)
            
    
    fig.show()
    
    
    # Save figure
    fig.savefig(f'{home}/ONC-EEW/manuscript/figures/DATAVERSE/{rupture}.png', dpi=300)
