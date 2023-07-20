#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 14:29:41 2023

@author: tnye
"""

# Imports
from mudpy.hfsims import windowed_gaussian,apply_spectrum
from mudpy.forward import gnss_psd
import numpy as np
from obspy import read

def add_synthetic_gnss_noise(st_E, st_N, st_Z, percentile=50):
    
    dt = st_E[0].stats.delta
    
    if dt == 0.5:
        st_E[0].decimate(factor=2)
        st_N[0].decimate(factor=2)
        st_Z[0].decimate(factor=2)
        
        st_E[0].stats.delta = 1.0
        st_N[0].stats.delta = 1.0
        st_Z[0].stats.delta = 1.0
    
    dt = st_E[0].stats.delta
    duration=st_E[0].stats.npts*dt
    
    # get white noise
    std=1.0 # this is a dummy parameter give it any value
    E_noise=windowed_gaussian(duration,dt,std=std,window_type=None)
    N_noise=windowed_gaussian(duration,dt,std=std,window_type=None)
    Z_noise=windowed_gaussian(duration,dt,std=std,window_type=None)
    
    f,Epsd,Npsd,Zpsd=gnss_psd(level=percentile,return_as_frequencies=True,return_as_db=False)
    
    #Covnert PSDs to amplitude spectrum
    Epsd = Epsd**0.5
    Npsd = Npsd**0.5
    Zpsd = Zpsd**0.5
    
    #apply the spectrum
    E_noise=apply_spectrum(E_noise,Epsd,f,dt,is_gnss=True)
    N_noise=apply_spectrum(N_noise,Npsd,f,dt,is_gnss=True)
    Z_noise=apply_spectrum(Z_noise,Zpsd,f,dt,is_gnss=True)
    
    #Remove mean for good measure
    E_noise -= np.mean(E_noise)
    N_noise -= np.mean(N_noise)
    Z_noise -= np.mean(Z_noise)
    
    st_E_noisy = st_E.copy()
    st_N_noisy = st_N.copy()
    st_Z_noisy = st_Z.copy()
    
    st_E_noisy[0].data = st_E_noisy[0].data + E_noise[:-1]
    st_N_noisy[0].data = st_N_noisy[0].data + N_noise[:-1]
    st_Z_noisy[0].data = st_Z_noisy[0].data + Z_noise[:-1]

    return(st_E_noisy,st_N_noisy,st_Z_noisy,E_noise[:-1],N_noise[:-1],Z_noise[:-1])
