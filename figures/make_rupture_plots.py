#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 21:18:18 2022

@author: tnye
"""

# Imports
import numpy as np
from glob import glob
from os import chdir, path, mkdir
import subprocess
from mudpy import gmttools

chdir('/Users/tnye/ONC/GMT')

batch = 'cascadia_longer_wfs'

rupt_files = sorted(glob(f'/Users/tnye/ONC/simulations/{batch}/output/ruptures/*.rupt'))

# Make output folder
if not path.exists(f'/Users/tnye/ONC/rupture_models/{batch}'):
    mkdir(f'/Users/tnye/ONC/rupture_models/{batch}')

if not path.exists(f'/Users/tnye/ONC/figures/png/rupture_models/{batch}'):
    mkdir(f'/Users/tnye/ONC/figures/png/rupture_models/{batch}')

for rupt_file in rupt_files:
    
    rupt_data = np.genfromtxt(rupt_file)
    
    log_file = rupt_file.replace('.rupt','.log')
    f = open(log_file, "r")
    lines = f.readlines()
    
    hypline = lines[16]
    lon,lat,depth = lines[16].split(':')[1].strip('\n').strip(')').strip(' (').split(',')
    lon = float(lon)
    lat = float(lat)
    run = rupt_file.split('/')[-1].strip('.rupt')
    run_num = run.split('.')[-1]
    
    f = open(f'/Users/tnye/ONC/rupture_models/{run}_epi.txt', 'w')
    f.write(f'{lon},{lat}')
    f.close()

    
    mag = round(float(lines[15].split(' ')[-1].strip('\n')),1)
    
    ss_slip = rupt_data[:,8]
    ds_slip = rupt_data[:,9]
    
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

    # Make .xy file
    meshfile = '/Users/tnye/ONC/data/fault_info/cascadia_slab2_40km.mshout'
    outfile = f'/Users/tnye/ONC/rupture_models/{batch}/{run}.xy'

    gmttools.triangular_rupt_2_gmt(meshfile,rupt_file,outfile)

    
    subprocess.run(['./plot_one_rupture.gmt',batch,f'{max_slip}',f'{mag}',f'{run}',f'{run_num}',f'{lon}',f'{lat}',f'{inc}',f'{scale}'])

