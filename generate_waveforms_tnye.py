#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 10:23:57 2022

@author: tnye
"""

# Imports
from mudpy import fakequakes,runslip,forward
import numpy as np
from obspy.core import UTCDateTime
from shutil import copy, copy2
import os
import time

########                            GLOBALS                             ########
#home='/home/tnye/onc/simulations/'
#project_name='batch_09'
home='/home/tnye/onc/simulations/'
project_name='cascadia_longer_wfs'
run_name='cascadia'
################################################################################


##############         What do you want to do?
init=0
make_ruptures=0
make_GFs=0
make_synthetics=0
make_waveforms=1
make_hf_waveforms=0
match_filter=0
# Things that only need to be done once
load_distances=0
G_from_file=1


# DEFINE PARAMETERS
# Runtime parameters 
ncpus=16                                     # how many CPUS you want to use for parallelization (needs ot be at least 2)
Nrealizations=4                              # Number of fake ruptures to generate per magnitude bin
hot_start=0                                  # If code quits in the middle of running, it will pick back up at this index

# File parameters
model_name='cascadia.mod'                    # Velocity model file name
#fault_base = 'cascadia_north45_25km'
fault_base='cascadia_slab2_40km'       # Fault model name
#fault_base='leech_river'
fault_name = f'{fault_base}.fault'
mean_slip_name=None                          # Set to path of .rupt file if patterning synthetic runs after a mean rupture model
rupture_list='ruptures.list'                 # Name of list of ruptures that are used to generate waveforms. 'ruptures.list' uses the full list of ruptures FakeQuakes creates. If you create file with a sublist of ruptures, use that file name.
distances_name='cascadia'                        # Name of matrix with estimated distances between subfaults i and j for every subfault pair


# Source parameters
#UTM_zone='10T' # UTM_zone for rupture region 
UTM_zone='10U'
time_epi=UTCDateTime('2023-03-22T00:00:00Z') # Origin time of event (can set to any time, as long as it's not in the future)
target_Mw=np.array([6.8])
#target_Mw=np.linspace(6.8,9.5,28)     # Desired magnitude(s), can either be one value or an array
#target_Mw=np.arange(4.5,7.0,0.5)
#target_Mw = np.array([4.5, 5, 5.5, 6])
#target_Mw = np.array([6.5,7.0,7.5,8.0,8.5,9,9.5])
#target_Mw = np.array([6,6.5,7,7.5,8,8.5,9])
#target_Mw = np.array([6])

hypocenter=None                              # Coordinates of subfault closest to desired hypocenter, or set to None for random
force_hypocenter=False                       # Set to True if hypocenter specified
rake=90                                      # Average rake for subfaults
scaling_law='T'                              # Type of rupture: T for thrust, S for strike-slip, N for normal
force_magnitude=False                       # Set to True if you want the rupture magnitude to equal the exact target magnitude
force_area=False                            # Set to True if you want the ruptures to fill the whole fault model

force_lat=False
#latmin=44.5
#latmax=50.8
latmin=None
latmax=None

# Correlation function parameters
hurst=0.4                                    # Hurst exponent form Melgar and Hayes 2019
Ldip='auto'                                  # Correlation length scaling: 'auto' uses Melgar and Hayes 2019, 'MB2002' uses Mai and Beroza 2002
Lstrike='auto'                               # Same as above
slip_standard_deviation=0.9                  # Standard deviation for slip statistics: Keep this at 0.9
lognormal=True                               # Keep this as True to solve the problem of some negative slip subfaults that are produced

# Rupture propagation parameters
rise_time = 'MH2017' 
rise_time_depths=[10,15]                     # Transition depths for rise time scaling (if slip shallower than first index, rise times are twice as long as calculated)
#rise_time_depths=[1,2]
max_slip=60                                  # Maximum sip (m) allowed in the model
max_slip_rule=False                          # If true, uses a magntidude-depence for max slip
shear_wave_fraction_shallow=0.49             # 0.8 is a standard value (Mai and Beroza 2002)
shear_wave_fraction_deep=0.8
source_time_function='dreger'                # options are 'triangle' or 'cosine' or 'dreger'
stf_falloff_rate=4                           # Only affects Dreger STF, 4-8 are reasonable values
num_modes=200                                 # Number of modes in K-L expansion
slab_name='cas_slab2_dep_02.24.18.xyz'       # Slab 2.0 Ascii file for 3D geometry, set to None for simple 2D geometry
#slab_name='leech_river_slab.xyz'
#mesh_name='cascadia_slab2_40km.mshout'       # GMSH output file for 3D geometry, set to None for simple 2D geometry
mesh_name = f'{fault_base}.mshout'

# WAVEFORM PARAMETERS
# Green's Functions parameters
GF_list='pnsn.gflist'                          # Stations file name
G_name='pnsn'                                  # Basename you want for the Green's functions matrices

# fk parameters
# used to solve wave eq in frequency domain 
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20           # Should be set to 0.1, 0, 1, 20
custom_stf=None                              # Assumes specified source time function above if set to None

# Low frequency waveform parameters
dt=0.5                                       # Sampling interval of LF data 
NFFT=1024                                     # Number of samples in LF waveforms (should be in powers of 2)
                                             # dt*NFFT = length of low-frequency dispalcement record
                                             # want this value to be close to duration (length of high-frequency record)

# High frequency waveform parameters
stress_parameter=50                          # Stress drop measured in bars (standard value is 50)
moho_depth_in_km=30.0                        # Average depth to Moho in this region 
Pwave=True                                   # Calculates P-waves as well as S-waves if set to True, else just S-Waves
kappa=None                                   # Station kappa values: Options are GF_list for station-specific kappa, a singular value for all stations, or the default 0.04s for every station if set to None
hf_dt=0.01                                   # Sampling interval of HF data
duration=510                                 # Duration (in seconds) of HF record
high_stress_depth=30                         # Doesn't do anything, but still shows up as a parameter. Set to whatever you want. 

# Intrinsic attenatuation parameters
Qexp=0.8
Qmethod='no_moho'

# Scattering attenuation parameters
scattering='on'
Qc_exp=0
baseline_Qc=150

# Match filter parameters
zero_phase=False                              # If True, filters waveforms twice to remove phase, else filters once
order=4                                      # Number of poles for filters
fcorner_low=0.998                                # Corner frequency at which to filter waveforms (needs to be between 0 and the Nyquist frequency)
fcorner_high=0.1
###############################################################################

start = time.time()

#Initalize project folders
if init==1:
    fakequakes.init(home,project_name)
copy2(f'/home/tnye/onc/files/{model_name}', home+project_name+'/structure')
copy2(f'/home/tnye/onc/files/{fault_name}', home+project_name+'/data/model_info')
copy2('/home/tnye/onc/files/' + GF_list, home+project_name+'/data/station_info')
copy2(f'/home/tnye/onc/files/{slab_name}', home+project_name+'/data/model_info')
copy2(f'/home/tnye/onc/files/{mesh_name}', home+project_name+'/data/model_info')

if load_distances==1:
    source_folder = f'/home/tnye/onc/data/distances/{fault_base}/'
    destination_folder = home+project_name+'/data/distances/'
    for file_name in os.listdir(source_folder):
        source = source_folder + file_name
        destination = destination_folder + file_name
        copy(source, destination)

#if G_from_file==1:
    #source_folder = f'/home/tnye/onc/data/GFs/{fault_base}.{G_name}/'
    #destination_folder = home+project_name+'/GFs/matrices/'
    #for file_name in os.listdir(source_folder):
     #   source = source_folder + file_name
      #  destination = destination_folder + file_name
       # copy(source, destination)


# Make a text file with the parameters
print('Writing parameter file list')
w = open(f'{home}{project_name}/parameter_file.txt', 'w')
w.write(f'model_name: {model_name}\n')
w.write(f'fault_name: {fault_name}\n')
w.write(f'mesh_name: {mesh_name}\n')
w.write(f'shear_wave_fraction_shallow: {shear_wave_fraction_shallow}\n')
w.write(f'shear_wave_fraction_deep: {shear_wave_fraction_deep}\n')
w.write(f'stations: {GF_list}\n')
w.write(f'stress_parameter: {stress_parameter}\n')
w.write(f'Qexp: {Qexp}\n')
w.write(f'fcorner_low: {fcorner_low}\n')
w.write(f'fcorner_high: {fcorner_high}')
w.close()

#Generate rupture models
if make_ruptures==1: 
    #fakequakes.generate_ruptures(home,project_name,run_name,fault_name,slab_name,mesh_name,load_distances,
    #    distances_name,UTM_zone,target_Mw,model_name,hurst,Ldip,Lstrike,num_modes,Nrealizations,rake,
    #    rise_time_depths,time_epi,max_slip,source_time_function,lognormal,slip_standard_deviation,scaling_law,
    #    ncpus,mean_slip_name=mean_slip_name,force_magnitude=force_magnitude,force_area=force_area,
    #    hypocenter=hypocenter,force_hypocenter=force_hypocenter,shear_wave_fraction_shallow=shear_wave_fraction_shallow,
    #    shear_wave_fraction_deep=shear_wave_fraction_deep,max_slip_rule=max_slip_rule)    
    
    fakequakes.generate_ruptures(home,project_name,run_name,fault_name,slab_name,mesh_name,load_distances,
                distances_name,UTM_zone,target_Mw,model_name,hurst,Ldip,Lstrike,num_modes,Nrealizations,rake,rise_time,
                rise_time_depths,time_epi,max_slip,source_time_function,lognormal,slip_standard_deviation,scaling_law,
                ncpus,mean_slip_name=mean_slip_name,force_magnitude=force_magnitude,force_area=force_area,
                hypocenter=hypocenter,force_hypocenter=force_hypocenter,shear_wave_fraction_shallow=shear_wave_fraction_shallow,
                shear_wave_fraction_deep=shear_wave_fraction_deep,max_slip_rule=max_slip_rule)

ncpus = 16

# Prepare waveforms and synthetics       
if make_GFs==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,None,fault_name,model_name,
        dt,None,NFFT,None,make_GFs,make_synthetics,dk,pmin,
        pmax,kmax,0,time_epi,hot_start,ncpus,custom_stf,impulse=True) 

#Make low frequency waveforms
if make_waveforms==1:
    forward.waveforms_fakequakes(home,project_name,fault_name,rupture_list,GF_list,
                model_name,run_name,dt,NFFT,G_from_file,G_name,source_time_function,
                stf_falloff_rate)

#Make semistochastic HF waveforms         
if make_hf_waveforms==1:
    #forward.hf_waveforms(home,project_name,fault_name,rupture_list,GF_list,
    #            model_name,run_name,dt,NFFT,G_from_file,G_name,rise_time_depths,
    #            moho_depth_in_km,ncpus,Qexp,source_time_function=source_time_function,
    #            duration=duration,stf_falloff_rate=stf_falloff_rate,hf_dt=hf_dt,
    #            Pwave=Pwave,hot_start=hot_start,stress_parameter=stress_parameter,
    #            high_stress_depth=high_stress_depth,kappa=kappa)

    forward.hf_waveforms(home,project_name,fault_name,rupture_list,GF_list,
                model_name,run_name,dt,NFFT,G_from_file,G_name,rise_time_depths,
                moho_depth_in_km,ncpus,source_time_function=source_time_function,
                duration=duration,stf_falloff_rate=stf_falloff_rate,hf_dt=hf_dt,
                Pwave=Pwave,hot_start=hot_start,stress_parameter=stress_parameter,
                high_stress_depth=high_stress_depth,kappa=kappa,Qexp=Qexp,
                Qmethod=Qmethod,scattering=scattering,Qc_exp=Qc_exp,baseline_Qc=baseline_Qc)

# Combine LF and HF waveforms with match filter                              
if match_filter==1:
    forward.match_filter(home,project_name,fault_name,rupture_list,GF_list,
            zero_phase,order,fcorner_low,fcorner_high)

end = time.time()
#print ("2 Mag 9 LF, 6 stations, 16 CPUs:", end - start)

