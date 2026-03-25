#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 00:31:29 2023

@author: tnye
"""

#%%

import obspy
from obspy import *
import os
import pandas as pd
import numpy as np


def SUBROUTINE_TO_GROUND(tr,G=1):
    '''
    Input: raw waveform (here an obspy trace)
           Station gain
    
    Output:
        Displacement (D)
        Velocity (V)
        Acceleration (A)
        Predominant period (Tp)
        Noise Level (N)
        Signal Level (S)
    '''
    # Constants:
    
    Q = 0.994       # High pass filter onstant
    B = 2/(1+Q)
    #G = get_gain_from_ShakeAlert_channel_file(tr)
    
    dt      = tr.stats.delta
    df      = tr.stats.sampling_rate
    channel = tr.stats.channel
    
    offset = 0
    trig_data = tr.data[offset:len(tr.data)]
    
    Tsig   = 0.05        # signal timescale (sec)
    Tnoise = 30       # noise timescale (sec)
    Tk     = 3           # relative importance of values vs derivative  
    
    alphaS = (Tsig-dt)/float(Tsig)
    alphaN = (Tnoise - dt)/float(Tnoise)
    
    A_series = [0.0];  V_series = [0.0];  D_series = [0.0]
    S        =[ 0.0];  N        = [0.0]
    
    if (channel[1] == 'N') or (channel[1]=='L'):
        # Channel is Acceleration
        for i in range(1,len(trig_data)):
            A_series.append(   ((trig_data[i] - trig_data[i-1]) / (B*G))  + Q*A_series[i-1])
            V_series.append(   ((A_series[i]+A_series[i-1])*dt /  (2*B))  + Q*V_series[i-1]  )
            D_series.append(   ((V_series[i]+V_series[i-1])*dt/   (2*B))  + Q*D_series[i-1]  )
        
    elif channel[1] == 'H':
        # Channel is Velocity
        for i in range(1,len(trig_data)):
            A_series.append(   ((trig_data[i] - trig_data[i-1])     / (dt*G)) )
            V_series.append(   ((trig_data[i] - trig_data[i-1])     / (G*B))  + Q*V_series[i-1]  )
            D_series.append(   ((V_series[i] +   V_series[i-1])*dt  / (2*B))  + Q*D_series[i-1]  )
            
    
    # High-passed raw value (not output)
    R_E2= [0.0]
    for i in range(1,len(trig_data)):
        R_E2.append(Q*R_E2[i-1] + (trig_data[i] - trig_data[i-1]))   
        Cf = R_E2[i]**2 + Tk*((trig_data[i] - trig_data[i-1])**2)
        S.append(alphaS*S[i-1] + (1-alphaS)*Cf)  
        N.append(alphaN*N[i-1] + (1-alphaN)*Cf)
        
      
    return(D_series,V_series,A_series,N,S)


def SUBROUTINE_PICKER(V_series,tr):
    
    dt      = tr.stats.delta
    df      = tr.stats.sampling_rate
    offset=0
    trig_data = tr.data[offset:len(tr.data)]
    Q = 0.994       # High pass filter onstant
    trig_level = 20
    
    
    # short term average (not output)
    Tsta = 0.05     # short term average timescale (sec)
    Tlta = 5        # long term average timetimescale (sec)
    Tk   = 3        # relative importance of values vs derivative
    
    alphaSTA = (Tsta - dt)/Tsta
    alphaLTA = (Tlta - dt)/Tlta
    
    STA   = [0.0]
    LTA   = [0.0]
    R     = [0.0]
    E2snr = [0.0]
    
    for i in range(1,len(trig_data)):
        R.append(Q*R[i-1] + (V_series[i] - V_series[i-1]))   
        Cf = R[i]**2 +  Tk*(V_series[i] - V_series[i-1])**2
        
        STA.append( alphaSTA*STA[i-1] + (1-alphaSTA)*Cf )
        LTA.append( alphaLTA*LTA[i-1] + (1-alphaLTA)*Cf )
    
    
    for i in range(1,len(STA)):
        try:
            E2snr.append(float(STA[i])/float(LTA[i]))
        except Exception as e:
            print(e, i)
            E2snr.append(0.0)
    #-----------------------------------------------------------------------------#
    #                                TRIGGERING
    
    trig_yn = 0;  trig_time = 0; trig_pos=-1; t_trig = 0
    for k in range(int(Tlta*df),len(E2snr)):
        if E2snr[k]>trig_level:
            t_trig = tr.stats.starttime+(k/tr.stats.sampling_rate)
            trig_yn = 1
            trig_time = (t_trig-tr.stats.starttime+(offset/tr.stats.sampling_rate))
            trig_pos =  (t_trig-tr.stats.starttime+offset)*tr.stats.sampling_rate
            
            break  
        else:
            trig_yn = 0
            continue
    return(t_trig,trig_time,trig_pos,E2snr)
            
           
            
           
def SUBROUTINE_TRIGGER_PARAMS(D_series,V_series,A_series,N,S,trig_time,trig_pos,tr):
   
    
    offset=0
    trig_data = tr.data[offset:len(tr.data)]
    dt      = tr.stats.delta
    df      = tr.stats.sampling_rate

    Dur = 4   # duration of reporting (s)
     
    snr_log = [0.0]
    pd_log  = [0.0]
   
    pd = 0;  pv=0; pa=0;
    pd_time = 0;  pv_time = 0; pa_time = 0;
    pd_snr = 0; 
    for i in range(1,len(trig_data)):
        t_trig = tr.stats.starttime+(i/tr.stats.sampling_rate)
      
        
        
        if i < int(trig_pos + (tr.stats.sampling_rate*.2)):
            snr_log.append(S[i]/N[i])
        elif (i >= int(trig_pos + tr.stats.sampling_rate*.2))  &  (i < trig_pos + (df*Dur) ):
            snr_log.append(S[i]/N[int(trig_pos)])
            
            if np.abs(D_series[i]) > pd:
                pd = np.abs(D_series[i])
                
                
                
                pd_time = (t_trig-tr.stats.starttime+(offset/tr.stats.sampling_rate))
                pd_snr = snr_log[i]
                
            if np.abs(V_series[i]) > pv:
                pv = np.abs(V_series[i])
                pv_time = (t_trig-tr.stats.starttime+(offset/tr.stats.sampling_rate))
                
            if np.abs(A_series[i]) > pa:
                pa = np.abs(A_series[i])
                pa_time = (t_trig-tr.stats.starttime+(offset/tr.stats.sampling_rate))
                
                
        else:
            snr_log.append(S[i]/N[int(trig_pos)])
        pd_log.append(pd)
        
        
        
        
    return(snr_log,pd_log,pd,pd_time,pd_snr,pa,pa_time,pv,pv_time)


#%%
# #############################################################################

# from obspy.signal.trigger import classic_sta_lta
# import tsueqs_main_fns as tmf

# E_file = '/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-PNSN-Onshore-StrongMotion/Cas-PNSN_SignalwithNoise/Cas-PNSN-SigNoise_cascadia-000020/Cas-PNSN-SigNoise_cascadia-000020_ALSE-HNE.mseed'
# N_file = '/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-PNSN-Onshore-StrongMotion/Cas-PNSN_SignalwithNoise/Cas-PNSN-SigNoise_cascadia-000020/Cas-PNSN-SigNoise_cascadia-000020_ALSE-HNN.mseed'
# Z_file = '/Users/tnye/ONC/simulations/cascadia_longer_wfs/waveforms_data_curation/Cas-PNSN-Onshore-StrongMotion/Cas-PNSN_SignalwithNoise/Cas-PNSN-SigNoise_cascadia-000020/Cas-PNSN-SigNoise_cascadia-000020_ALSE-HNZ.mseed'

# st = Stream()
# st.append(read(N_file)[0])
# st.append(read(E_file)[0])
# st.append(read(Z_file)[0])

# st[0].stats.channel = 'HHN'
# st[1].stats.channel = 'HHE'
# st[2].stats.channel = 'HHZ'

# E_vel =  tmf.highpass(tmf.accel_to_veloc(tmf.highpass(read(E_file),1/15,100,4,zerophase=False)),1/15,100,4,zerophase=False)
# N_vel =  tmf.highpass(tmf.accel_to_veloc(tmf.highpass(read(N_file),1/15,100,4,zerophase=False)),1/15,100,4,zerophase=False)
# Z_vel =  tmf.highpass(tmf.accel_to_veloc(tmf.highpass(read(Z_file),1/15,100,4,zerophase=False)),1/15,100,4,zerophase=False)

# pick_times = []


# for ii in range(len(st)):
    
    # tr = st[ii]
    # G = 1 # YOU NEED TO HAVE A GAIN VALUE HERE FOR
    # (D_series,V_series,A_series,N,S)  =    SUBROUTINE_TO_GROUND(tr,G)
    # (trig_time,trig_pos,stalta)       =    SUBROUTINE_PICKER(V_series,tr)
    # # if trig_pos > 0:
    # #     (snr_log,pd_log,pd,pd_time,pd_snr,pa,pa_time,pv,pv_time) =    SUBROUTINE_TRIGGER_PARAMS(D_series,V_series,A_series,N,S,trig_time,trig_pos,tr)
    # pick_times.append(trig_time)


# cft_E = classic_sta_lta(E_vel[0].data, int(0.05 * 100), int(5 * 100))
# cft_N = classic_sta_lta(N_vel[0].data, int(0.05 * 100), int(5 * 100))
# cft_Z = classic_sta_lta(Z_vel[0].data, int(0.05 * 100), int(5 * 100))

# obspy_pick = []

# # Determine if algorithm was triggered
# if len(np.where(cft_E>=20)[0])>0:
#     pick = E_vel[0].times()[np.where(cft_E>=20)[0][0]]
#     obspy_pick.append(pick)
# else:
#     obspy_pick.append(0)
# if len(np.where(cft_N>=20)[0])>0:
#     pick = N_vel[0].times()[np.where(cft_N>=20)[0][0]]
#     obspy_pick.append(pick)
# else:
#     obspy_pick.append(0)
# if len(np.where(cft_Z>=20)[0])>0:
#     pick = Z_vel[0].times()[np.where(cft_Z>=20)[0][0]]
#     obspy_pick.append(pick)
# else:
#     obspy_pick.append(0)



# # plt.figure
# # plt.plot(st[0].times(),V_series,alpha=0.7,label='EPIC')
# # # plt.plot(st_N_vel_filt[0].times(),st_N_vel_filt[0].data,alpha=0.7,label='integrate')
# # plt.legend()


