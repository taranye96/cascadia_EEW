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
import csv
import math
import numpy

chdir('/Users/tnye/ONC/GMT')

batch = 'cascadia_longer_wfs'

rupt_files = sorted(glob(f'/Users/tnye/ONC/simulations/{batch}/output/ruptures/*.rupt'))

# Make output folder

if not path.exists(f'/Users/tnye/ONC/figures/png/rupture_models/{batch}/onset'):
    mkdir(f'/Users/tnye/ONC/figures/png/rupture_models/{batch}/onset')

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
    
    mag = round(float(lines[15].split(' ')[-1].strip('\n')),1)
    
    onset = rupt_data[:,12]
    
    max_onset = np.max(onset)
    
    if max_onset < 0.5:
        inc = 0.1
    elif max_onset>=0.5 and max_onset<1:
        inc = 0.5
    elif max_onset>=1 and max_onset<5:
        inc = 1
    elif max_onset>=5 and max_onset<10:
        inc = 5
    elif max_onset>=10 and max_onset<50:
        inc = 10
    elif max_onset>=50 and max_onset<100:
        inc = 25
    elif max_onset>=100 and max_onset<200:
        inc=50
    else:
        inc=100
        
    #Read mesh file first
    meshfile = '/Users/tnye/ONC/data/fault_info/cascadia_slab2_40km.mshout'
    outfile = f'/Users/tnye/ONC/rupture_models/{batch}/onset/{run}.xy'
    
    def triangular_onset_2_gmt(meshfile,rupt_file,outfile,kinematic_out_folder=None,percentage=0,is_total_model=True):
        meshnum = list()
        meshlon1 = list()
        meshlon2 = list()
        meshlon3 = list()
        meshlat1 = list()
        meshlat2 = list()
        meshlat3 = list()
        meshdep1 = list()
        meshdep2 = list()
        meshdep3 = list()
        faultarea = list()
        with open(meshfile, 'r') as f:
            next(f)
            reader = csv.reader(f,delimiter='\t')
            for row in reader:
                meshnum.append(int(row[0]))
                meshlon1.append(float(row[4]))
                meshlat1.append(float(row[5]))
                meshdep1.append(float(row[6]))
                meshlon2.append(float(row[7]))
                meshlat2.append(float(row[8]))
                meshdep2.append(float(row[9]))
                meshlon3.append(float(row[10]))
                meshlat3.append(float(row[11]))
                meshdep3.append(float(row[12]))
                faultarea.append(float(row[14]))
        
        
        #Read inverse file
        
        invnum = list()
        ss = list()
        ds = list()
        rupttime = list()
        risetime = list()
        duration = list()
        rig = list()
        
        
        with open(rupt_file, 'r') as f:
            next(f)
            reader = csv.reader(f,delimiter='\t')
            for row in reader:
                invnum.append(int(row[0]))
                ss.append(float(row[8]))
                ds.append(float(row[9]))
                rupttime.append(float(row[12]))
                risetime.append(float(row[6]))
                duration.append(float(row[7]))
                if is_total_model:
                    rig.append(float(row[12]))
                else:
                    rig.append(float(row[13]))
        
        
        INVN = numpy.asarray(invnum)
        MESN = numpy.asarray(meshnum)
        
        SS = numpy.asarray(ss)
        DS = numpy.asarray(ds)
        
        RT = numpy.asarray(rupttime)
        DR = numpy.asarray(duration)
        
        TOTS = numpy.sqrt(numpy.power(SS,2)+numpy.power(DS,2))
        FA = numpy.asarray(faultarea)
        RIG = numpy.asarray(rig)
        
        
        #Total slip model
        moment = 0
        fso = open(outfile,'w')
        slip_threshold=(percentage/100)*TOTS.max()
        for i in range(0, numpy.amax(MESN)):
            a1 = numpy.where(MESN[i] == INVN)[0]
            totslip = numpy.sum(TOTS[a1])
            lon1 = "{0:.4f}".format(meshlon1[i])
            lon2 = "{0:.4f}".format(meshlon2[i])
            lon3 = "{0:.4f}".format(meshlon3[i])
            lat1 = "{0:.4f}".format(meshlat1[i])
            lat2 = "{0:.4f}".format(meshlat2[i])
            lat3 = "{0:.4f}".format(meshlat3[i])
            dep1 = "{0:.4f}".format(meshdep1[i])
            dep2 = "{0:.4f}".format(meshdep2[i])
            dep3 = "{0:.4f}".format(meshdep3[i])
            on = "{0:.4f}".format(onset[i])
            fso.write('> -Z'+on+'\n')
            fso.write(lon1+' '+lat1+' '+dep1+'\n')
            fso.write(lon2+' '+lat2+' '+dep2+'\n')
            fso.write(lon3+' '+lat3+' '+dep3+'\n')
            fso.write(lon1+' '+lat1+' '+dep1+'\n')
        
        fso.close()
    
    triangular_onset_2_gmt(meshfile,rupt_file,outfile)

    subprocess.run(['./plot_onset.gmt',batch,f'{max_onset}',f'{mag}',f'{run}',f'{run_num}',f'{lon}',f'{lat}',f'{inc}'])

