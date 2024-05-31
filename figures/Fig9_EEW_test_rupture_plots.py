#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 21:18:18 2022

@author: tnye
"""

##############################################################################
# This script makes Figure 9 in the paper, which shows slip patterns of the
# test ruptures used in the EEW performance analysis. This script calls the
# GMT script plot_test_ruptures.gmt.
##############################################################################

# Imports
import numpy as np
from os import path, makedirs
import subprocess
from mudpy import gmttools

batch = 'cascadia'

rupt_files = ['/Users/tnye/ONC/simulations/cascadia/output/ruptures/cascadia.000015.rupt',
              '/Users/tnye/ONC/simulations/cascadia/output/ruptures/cascadia.000018.rupt',
              '/Users/tnye/ONC/simulations/cascadia/output/ruptures/cascadia.000051.rupt',
              '/Users/tnye/ONC/simulations/cascadia/output/ruptures/cascadia.000065.rupt',
              '/Users/tnye/ONC/simulations/cascadia/output/ruptures/cascadia.000096.rupt',
              '/Users/tnye/ONC/simulations/cascadia/output/ruptures/cascadia.000108.rupt']

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
    
    f = open(f'/Users/tnye/ONC/GMT/plotting_files/rupture_epi/{batch}/{run}_epi.txt', 'w')
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
    
    subprocess.run(['./plot_test_ruptures.gmt',batch,f'{max_slip}',f'{mag}',f'{run}',f'{run_num}',f'{lon}',f'{lat}',f'{inc}',f'{scale}'])

