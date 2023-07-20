#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 15:40:24 2022

@author: tnye
"""

##############################################################################

##############################################################################

def get_Rp(stlon, stlat, stelev, rupture, exp=-2.3):
    """
    Calculates generalized mean rupture distance (Rp) for a station.
    
    Inputs:
        stlon(float): Station longitude
        stlat(float): Station latitude
        stelev(float): Station elevation (km)
        rupt_file(string): Path to .rupt file
        exp(float): Power of the mean
    
    Returns:
        Rp(float): Generalized mean rupture distance  
    """
    
    import numpy as np
    import pandas as pd
    import tsueqs_main_fns as main
    
    rupt = np.genfromtxt(rupture)
    ind = np.where((rupt[:,8]!=0) & (rupt[:,9]!=0))[0]
    
    total_slip = np.sum(np.sqrt(rupt[:,8]**2+rupt[:,9]**2))
    
    weighted_dist = []
    for i in ind:
        w_i = np.sqrt(rupt[i,8]**2 + rupt[i,9]**2)/total_slip
        R_i = main.compute_rhyp(stlon,stlat,stelev,rupt[i,1],rupt[i,2],rupt[i,3])
        weighted_dist.append(w_i*(R_i**exp))
    
    Pd = np.sum(weighted_dist)**(1/exp)
    
    
    return(Pd)

def get_pgd_scaling(Mw, R, model):
    """
    Empirically estimates PGD from hypocentral distance using the scaling
    relation from Goldberg et al. (2021).
    
    Inputs:
        MW(float): Moment magnitude
        R(float): Hypocentral distance in km. 
    
    Returns:
        PGD(float): Peak ground displacement (m) 
    """
    
    import numpy as np
    
    if model == 'GA21':
        A = -5.902
        B = 1.303
        C = -0.168
        sigma = 0.255
        logpgd = A + B*Mw + C*Mw*np.log10(R)
        pgd = 10**logpgd
        pgd = pgd/100
    
    elif model == 'MA15':
        A = -4.434
        B = 1.047
        C = -0.138
        sigma = 0.27
        logpgd = A + B*Mw + C*Mw*np.log10(R)
        pgd = 10**logpgd
        pgd = pgd/100
        
    return(pgd, sigma)


def WGRW12(y, mode):
    '''
    Compute MMI with Worden 2012 given either PGA or PGV
    Input:
        y:      Array - Ground-motion intensity measure in cm/s/s if PGA 
                    or cm/s if PGV
        mode:   Integer - 0 if PGA, 1 if PGV
    Returns:
        MMI:    Array - Intensity measure in modified mercalli intensity
    '''
    import numpy as np
    
    if (mode == 0):
        pgalen = len(y)
        MMI = np.zeros(pgalen)
        pgalarge = np.where(np.log10(y) >= 1.57)[0]
        pgasmall = np.where((np.log10(y) < 1.57) & (y > 0))[0]
        pgazero = np.where(y == 0)[0]
        
        MMI[pgalarge] = 3.70*np.log10(y[pgalarge])-1.60
        MMI[pgasmall] = 1.55*np.log10(y[pgasmall])+1.78
        MMI[pgazero] = -10
    else:
        pgvlen = len(y)
        MMI = np.zeros(pgvlen)
        pgvlarge = np.where(np.log10(y) >= 0.53)[0]
        pgvsmall = np.where(np.log10(y) < 0.53)[0]
        try:
            MMI[pgvlarge] = 3.16*np.log10(y[pgvlarge])+2.89
        except:
            pass
        MMI[pgvsmall] = 1.47*np.log10(y[pgvsmall])+3.78
        
    return(MMI)


    