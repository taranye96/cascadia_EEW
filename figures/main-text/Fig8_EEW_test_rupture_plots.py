#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 09:43:28 2025

@author: tnye
"""

# Imports
import os
import numpy as np
import pygmt
import make_xy


home='/Users/tnye/Library/CloudStorage/OneDrive-DOI/UO-projects/ONC-EEW'

# Input files
slab2='/Users/tnye/slab_models/Cascadia/cas_slab2_dep_02.24.18.grd'
DEM=f'{home}/GMT/topo/cascadia30.grd'
ruptures_path=f'{home}/GMT/plotting_files'
oncew_coords=f'{home}/GMT/data/station_info/onc-ew_coords.txt'
shakealert_coords=f'{home}/GMT/data/station_info/shakealert_coords.txt'
overlap_coords_pnsn=f'{home}/GMT/data/station_info/overlap_coords_pnsn.txt'
overlap_coords_onc=f'{home}/GMT/data/station_info/overlap_coords_onc.txt'
epic_grid_kml='/Users/tnye/GoogleEarth/EPIC_detection_grid.kml'
epic_grid_txt='/Users/tnye/GoogleEarth/EPIC_detection_grid.txt'
topo_no=f'{home}/GMT/topo/topo_no.grd'
topo=f'{home}/GMT/topo/topo.grad'
slip_cpt=f'{home}/GMT/cpt/slip.cpt'
topo_cpt=f'{home}/GMT/cpt/topo.cpt'
bathy_cpt=f'{home}/GMT/cpt/bathy.cpt'
meshfile = f'{home}/files/fakequakes/model_info/cascadia_slab2_40km.mshout'

rupt_files = [f'{home}/simulations/cascadia/output/ruptures/cascadia-000015.rupt',
              f'{home}/simulations/cascadia/output/ruptures/cascadia-000018.rupt',
              f'{home}/simulations/cascadia/output/ruptures/cascadia-000051.rupt',
              f'{home}/simulations/cascadia/output/ruptures/cascadia-000065.rupt',
              f'{home}/simulations/cascadia/output/ruptures/cascadia-000096.rupt',
              f'{home}/simulations/cascadia/output/ruptures/cascadia-000108.rupt']

ruptures_path=f'{home}/GMT/plotting_files'

# ONC bounds
box_coords = [
    [-132.75, 46],
    [-122, 46],
    [-122, 52],
    [-132.75, 52],
    [-132.75, 46]  # close the polygon by repeating the first point
]


region = [-128.75, -121.5, 39.75, 51]
projection = "M1.4i"

rupture_grid = pygmt.datasets.load_earth_relief(resolution='01m', region=region)

def get_colorbar_parameters(ss_slip, ds_slip):
    
    total_slip = np.sqrt(ss_slip**2 + ds_slip**2)
    
    max_slip = np.max(total_slip)
    
    if max_slip < 1:
        inc = 0.5
        scale = 1
    elif max_slip>=1 and max_slip<5:
        inc = 2
        scale = 5
    elif max_slip>=5 and max_slip<10:
        inc = 5
        scale = 10
    elif max_slip >=10 and max_slip < 20:
        inc = 10
        scale = 20
    elif max_slip >=20 and max_slip < 50:
        inc = 20
        scale = 50
    elif max_slip >= 50:
        inc = 20
        scale = max_slip
    
    return(inc, scale)
    

# Set general config
fig = pygmt.Figure()
pygmt.config(
    FONT_ANNOT_PRIMARY="8p,Helvetica",

    MAP_FRAME_TYPE="plain",
    PS_MEDIA="a2",
    FORMAT_GEO_MAP="ddd",
    FONT_LABEL="16p"
)

with fig.subplot(nrows=2, ncols=3, subsize=("1.4i", "3.5i"), margins=["0.4c", "0.8c"]):
    
    for i, rupt_file in enumerate(rupt_files):
        
        # Read in rupture data
        rupt_data = np.genfromtxt(rupt_file)
        
        log_file = rupt_file.replace('.rupt','.log')
        f = open(log_file, "r")
        lines = f.readlines()
        
        hypline = lines[16]
        hyplon,hyplat,hypdepth = lines[16].split(':')[1].strip('\n').strip(')').strip(' (').split(',')
        hyplon = float(hyplon)
        hyplat = float(hyplat)
        run = rupt_file.split('/')[-1].strip('.rupt')
        run_num = run.split('-')[-1]
        
        mag = round(float(lines[15].split(' ')[-1].strip('\n')),1)
        
        # Make slip cpt
        rupt_cpt=f'{home}/GMT/cpt/{run}_slip.cpt'
        inc, scale = get_colorbar_parameters(rupt_data[:,8], rupt_data[:,9])
        pygmt.makecpt(cmap=slip_cpt, series=[0,scale], output=rupt_cpt)
        
        rupt_epi=f'{ruptures_path}/rupture_epi/cascadia/{run}_epi.txt'
                
        rupt_xy = f'{home}/GMT/plotting_files/rupture_xy/cascadia/{run}_active.xy'
        # make_xy.xy_active_faults(meshfile,rupt_file,rupt_xy)
        
        with fig.set_panel(i):
            
            # Plot base
            fig.basemap(region=region, projection=projection, frame=["WesN", "2p"])
            fig.grdimage(grid=rupture_grid, cmap=bathy_cpt, shading='+d', transparency=20, region=region, projection=projection)
            
            # Coastline and bathymetry
            fig.coast(water="197/212/253", land='white', transparency=30, region=region, projection=projection)
            fig.coast(shorelines='0.25p,black', transparency=0, region=region, projection=projection)
            fig.coast(borders=["1/0.25p,black", "2/0.25p,black"], region=region, projection=projection)
            
            # Rupture
            fig.plot(data=rupt_xy, fill="white", pen="white", region=region, projection=projection)
            fig.plot(data=rupt_xy, cmap=rupt_cpt, region=region, projection=projection, fill='+z')
            # fig.plot(data=rupt_xy, pen='0.01p,180/180/180', region=region, projection=projection)
            
            # Scalebar
            with pygmt.config(FONT_ANNOT_PRIMARY='16p'):
                fig.colorbar(cmap=rupt_cpt, frame=[f"xa{inc}+lSlip (m)"], position="JBC+o0c/0.5c+w1.6i/0.125i", projection=projection, region=region)
                
            # Run number
            fig.text(x=-127.8, y=40.1, text=f"M{mag}", font="10p,Helvetica,black", justify="CM", region=region, projection=projection)
            
            # Stations
            fig.plot(data=shakealert_coords, style="i0.19c", fill="53/176/82", pen="black", region=region, projection=projection)
            fig.plot(data=oncew_coords, style="i0.19c", fill="235/252/3", pen="black", region=region, projection=projection)
            fig.plot(data=overlap_coords_onc, style="i0.19c", fill="3/252/240", pen="black", region=region, projection=projection)
            fig.plot(data=overlap_coords_pnsn, style="i0.19c", fill="3/252/240", pen="black", region=region, projection=projection)
            
            # Event
            fig.plot(data=[[hyplon, hyplat]], style="a0.35c", pen="0.75p,black", region=region, projection=projection)
            
            # Detection bounds
            fig.plot(data=epic_grid_txt, pen=".8p,black,", region=region, projection=projection)
            fig.plot(data=box_coords, pen=".8p,black,--.", region=region, projection=projection)


fig.show()


# Save figure
fig.savefig(f'{home}/manuscript/figures/subfigs/Fig8_test_ruptures.png', dpi=500)

