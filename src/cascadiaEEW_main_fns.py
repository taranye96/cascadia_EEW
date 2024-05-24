#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 15:40:24 2022

@author: tnye
"""

##############################################################################
# This module contains all the functions called within the other scripts of
# this project.  
##############################################################################


def highpass(datastream,fcorner,fsample,order,zerophase=True):
    '''
    Make a highpass zero phase filter on stream object
    Input:
        datastream:                 Obspy stream object with data to filter
        fcorner:                    Corner frequency at which to highpass filter
        fsample:                    Sample rate in (Hz) - or 1/dt
        order:                      The numper of poles in the filter (use 2 or 4)
        zerophase:                  Boolean for whether or not the filter is zero phase
                                        (that does/doesn't advance or delay 
                                        certain frequencies). Butterworth filters 
                                        are not zero phase, they introduce a phase
                                        shift. If want zero phase, need to filter twice.
    Output:
        highpassedstream:           Obspy stream wth highpass filtered data
    '''
    from scipy.signal import butter,filtfilt,lfilter
    from numpy import array
    
    data = datastream[0].data
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    if zerophase==True:
        data_filt=filtfilt(b,a,data)
    else:
        data_filt=lfilter(b,a,data)
    
    
    ## Make a copy of the original stream object:
    highpassedstream = datastream.copy()
    
    ## Add the highpassed data to it:
    highpassedstream[0].data = data_filt
    
    return highpassedstream


def bandpass(datastream,lowcut,highcut,fsample,order,zerophase=True):
    """
    Applies a bandpass filter to an obspy stream.
    
    Inputs:
        datastream: Obspy stream
        lowcut(float): Lower corner frequency of filter
        highcut(float): Higher coner frequency of filter
        fsample(float): Sampling rate
        order(int): Number of poles
        zerophase(True/False): If True, filter run forwards and backwards to remove phase shift
        
    Return:
        bandpassedstream: Bandpass-filtered obspy stream
        
    """
    
    from scipy.signal import butter,filtfilt,lfilter
    
    data = datastream[0].data
    
    b, a = butter(order, [lowcut, highcut], fs=fsample, btype='band')
    if zerophase==True:
        data_filt=filtfilt(b,a,data)
    else:
        data_filt=lfilter(b,a,data)
    
    
    ## Make a copy of the original stream object:
    bandpassedstream = datastream.copy()
    
    ## Add the highpassed data to it:
    bandpassedstream[0].data = data_filt
    
    return bandpassedstream


def accel_to_veloc(acc_timeseries_stream):
    '''
    Integrate an acceleration time series to get velocity
    Input:
        acc_timeseries_stream:              Obspy stream object with a time series
                                                of acceleration, units in m/s/s
                                                baseline and gain corrected.
    Output:
        vel_timeseries_stream:          
    '''
    
    from scipy.integrate import cumtrapz
    
    ## Get the bsaeline corrected and gain corrected time series:
    acc_amplitude = acc_timeseries_stream[0].data
    
    ## And times:
    acc_times = acc_timeseries_stream[0].times()
    
    ## INtegrate the acceration time series to get velocity:
    vel_amplitude = cumtrapz(acc_amplitude,x=acc_times,initial=0)
    
    ## Make a copy of the old stream object:
    vel_timeseries_stream = acc_timeseries_stream.copy()
    
    ## Put the integrated velocity in there in place of accel:
    vel_timeseries_stream[0].data = vel_amplitude
    
    ## Return:
    return vel_timeseries_stream
  


def rotate_gnss_local(st_XYZ, stlon, stlat):
    """
    Transforms GNSS displacements from a global cartesian frame (X,Y,Z) into a 
    local North, East, Up frame (N,E,U).
    
    Inputs:
        st_XYZ: Obspy stream with three traces for the three channels.
        stlon(float): Station longitude
        stlat(float): Staiton latitude
    
    Return:
        st_NEU: Rotated GNSS stream
        
    """
    
    import numpy as np
    from numpy import sin, cos, deg2rad, dot
    
    # Make a stream for the rotated data
    st_NEU = st_XYZ.copy()
    
    # Set up rotation matrix
    R = np.array([[-sin(deg2rad(stlat))*cos(deg2rad(stlon)),
                  -sin(deg2rad(stlon))*sin(deg2rad(stlat)),
                  cos(deg2rad(stlat))],
                 [-sin(deg2rad(stlon)),
                  cos(deg2rad(stlon)),
                  0],
                 [cos(deg2rad(stlon))*cos(deg2rad(stlat)),
                  cos(deg2rad(stlat))*sin(deg2rad(stlon)),
                  sin(deg2rad(stlat))]])
    
    for i in range(len(st_XYZ[0].times())):
        
        dx = st_XYZ[0].data[i] - st_XYZ[0].data[0]
        dy = st_XYZ[1].data[i] - st_XYZ[1].data[0]
        dz = st_XYZ[2].data[i] - st_XYZ[2].data[0]
    
        [dn, de, du] = dot(R, np.array([dx,dy,dz]).T)
        
        st_NEU[0].data[i] = dn
        st_NEU[1].data[i] = de
        st_NEU[2].data[i] = du
        
    return(st_NEU)


def add_synthetic_gnss_noise(st_E, st_N, st_Z, percentile=50):
    
    from mudpy.forward import gnss_psd
    from mudpy.hfsims import windowed_gaussian,apply_spectrum
    from mudpy.forward import gnss_psd
    import numpy as np
    from obspy import read
    
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


def get_eucl_norm_3comp(E_record, N_record, Z_record):
    """
    Get the euclidean norm of the three components of a record.  This is
    equivalent to calculating the magnitude of a vector. 

    Inputs:
        E_record(array): East-West component trace data.
        N_record(array): North-South component trace data.
        Z_record(array): Vertical component trace data.
    
    Return:
        eucl_norm(array): Record of euclidian norm.
    """

    import numpy as np

    eucl_norm = np.sqrt(E_record**2 + N_record**2 + Z_record**2)

    return (eucl_norm)



def compute_rotd50(accel1, accel2):
    """
    Computes rotd50 of a timeseries
    
    Input:
        accel1(array): Timeseries 1 accelerations
        accel2(array): Timeseries 2 accelerations
    
    Return:
        rotd50_accel(array): 50th percentile rotated timeseries
    """
    
    import numpy as np
    
    percentile = 50
    
    accels = np.array([accel1,accel2])
    angles = np.arange(0, 180, step=0.5)
    
    # Compute rotated time series
    radians = np.radians(angles)
    coeffs = np.c_[np.cos(radians), np.sin(radians)]
    rotated_time_series = np.dot(coeffs, accels)
    
    # Sort this array based on the response
    peak_responses = np.abs(rotated_time_series).max(axis=1)
    
    # Get the peak response at the requested percentiles
    p_peak_50 = np.percentile(peak_responses, percentile)
    
    # Get index of 50th percentile
    ind = np.where(peak_responses==peak_responses.flat[np.abs(peak_responses - p_peak_50).argmin()])[0][0]
    
    # Get rotd50 timeseries
    rotd50_accel = rotated_time_series[ind]
    
    return(rotd50_accel)




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
    
    Rp = np.sum(weighted_dist)**(1/exp)
    
    
    return(Rp)


def compute_repi(stlon,stlat,hypolon,hypolat):
    '''
    Compute the hypocentral distance for a given station lon, lat and event hypo
    Input:
        stlon:          Float with the station/site longitude
        stlat:          Float with the station/site latitude
        hypolon:        Float with the hypocentral longitude
        hypolat:        Float with the hypocentral latitude
    Output:
        repi:           Float with the epicentral distance, in km
    '''
    
    from pyproj import Geod
    
    
    ## Make the projection:
    p = Geod(ellps='WGS84')
    
    ## Apply the projection to get azimuth, backazimuth, distance (in meters): 
    az,backaz,horizontal_distance = p.inv(stlon,stlat,hypolon,hypolat)

    ## Put them into kilometers:
    horizontal_distance = horizontal_distance/1000.
 
    ## Epicentral distance is horizontal distance:
    repi = horizontal_distance
    
    return repi


def compute_rhyp(stlon,stlat,stelv,hypolon,hypolat,hypodepth):
    '''
    Compute the hypocentral distance for a given station lon, lat and event hypo
    Input:
        stlon:          Float with the station/site longitude
        stlat:          Float with the station/site latitude
        stelv:          Float with the station/site elevation (m)
        hypolon:        Float with the hypocentral longitude
        hypolat:        Float with the hypocentral latitude
        hypodepth:      Float with the hypocentral depth (km)
    Output:
        rhyp:           Float with the hypocentral distance, in km
    '''
    
    import numpy as np
    from pyproj import Geod
    
    
    ## Make the projection:
    p = Geod(ellps='WGS84')
    
    ## Apply the projection to get azimuth, backazimuth, distance (in meters): 
    az,backaz,horizontal_distance = p.inv(stlon,stlat,hypolon,hypolat)

    ## Put them into kilometers:
    horizontal_distance = horizontal_distance/1000.
    stelv = stelv/1000
    ## Hypo deptha lready in km, but it's positive down. ST elevation is positive
    ##    up, so make hypo negative down (so the subtraction works out):
    hypodepth = hypodepth * -1
    
    ## Get the distance between them:
    rhyp = np.sqrt(horizontal_distance**2 + (stelv - hypodepth)**2)
    
    return rhyp

def get_pgd_scaling(Mw, R, model):
    """
    Empirically estimates PGD from hypocentral distance using the scaling
    relation from Goldberg et al. (2021).
    
    Inputs:
        MW(float): Moment magnitude
        R(float): Distance, either Rp for GA21 or Rhyp for MA15 
    
    Returns:
        PGD(float): Peak ground displacement (m) 
    """
    
    import numpy as np
    
    if model == 'GA21_joint':
        A = -5.902
        B = 1.303
        C = -0.168
        sigma = 0.255
        logpgd = A + B*Mw + C*Mw*np.log10(R)
        pgd = 10**logpgd
        pgd = pgd/100
        
    elif model == 'GA21_obs':
        A = -3.841
        B = 0.919
        C = -0.122
        sigma = 0.252
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


def polarization(st_Z, st_N, st_E, noise_Z, noise_N, noise_E):
    """
    Performs polarization analysis to isolate P- and S-waves, from Rosenberger 2019.
    
    Input:
        st_Z(Stream): Vertical component obspy stream
        st_N(Stream): N component obspy stream
        st_E(Stream): E component obspy stream
        noise_Z(Stream): Noise from vertical component obspy stream
        noise_N(Stream): Noise from N component obspy stream
        noise_E(Stream): Noise from E component obspy stream
    
    Output:
        st_P(Stream): Stream with 3 components of P-wave
        st_S(Stream): Stream with 3 components of S-wave
    """
    
    import numpy as np
    from obspy import read, Stream
    
    var = np.max([np.var(noise_E[0].data),np.var(noise_N[0].data),np.var(noise_Z[0].data)])

    M = len(st_E[0].data)

    # Initialize matrices U and S
    U = np.zeros((3, 2))
    S = np.zeros((2, 2))

    st_P = Stream()
    st_P.append(st_Z[0].copy())
    st_P.append(st_N[0].copy())
    st_P.append(st_E[0].copy())

    st_S = Stream()
    st_S.append(st_Z[0].copy())
    st_S.append(st_N[0].copy())
    st_S.append(st_E[0].copy())

    lambd = (M-1)/M

    p_array = np.zeros(M)

    for i in range(M):
       
        # Construct data vector d
        d = np.array([st_Z[0].data[i - 1], st_N[0].data[i - 1], st_E[0].data[i - 1]])

        # Simple initialization
        if i == 0:
            n2d = np.linalg.norm(d, 2)
            U[:, 0] = d / n2d  # Eq. 26
            S[0, 0] = n2d  # Eq. 27

        # Compute m and m_p
        m = np.dot(U.T, d)  # m || U, Eq. 16
        m_p = d - np.dot(U, m)  # projection _| U, Eq. 17
        
        r = np.dot(d.T, d)
        s = np.dot(m.T, m)
        t = np.dot(U, m)
        
        # Innovation
        pp = r - 2 * s + np.dot(t.T, t)
        
        if np.abs(pp) < var:
            pp = 0
        
        p = np.sqrt(pp)
        
        if p > 0:
            # Rank increase Q: 3x3, Eq. 19
            Q = np.array([[lambd * S[0, 0], m[0]],
                          [0, 0],
                          [p, 0]])
        else:
            # Rank preserving Q: 2x3, Eq. 22
            Q = np.array([[lambd * S[0, 0], m[0]],
                          [0, 0]])
        
        # Standard generic SVD for the diagonalization of Q, Eq. 24
        A, Ss, V = np.linalg.svd(Q)
        
        # Updating S and U
        S[:2, :2] = np.diag(Ss[:2])  # discard smallest singular value
        
        if p > 0:
            U = np.column_stack((U, m_p / p)).dot(A[:, :2])  # rank reduction
        else:
            U = U.dot(A[:, :2])  # rank preserving
        
        re = 1.0
        
        if S[0, 0] > 0:
            re = S[1, 1] / S[0, 0]  # recti-linearity
        
        incl = np.abs(U[0, 0])  # Eq. 29, cos(incl) would be the angle of incidence vs. Z in radians
        rn = np.dot(np.dot(U,U.T),d)
        Rn = 1 - (S[0,1]/S[0,0])
        
        pn = np.dot(incl*Rn,rn)
        st_P[0].data[i] = pn[0]
        st_P[1].data[i] = pn[1]
        st_P[2].data[i] = pn[2]
        
        sn = np.dot((1-incl)*Rn,rn)
        st_S[0].data[i] = sn[0]
        st_S[1].data[i] = sn[1]
        st_S[2].data[i] = sn[2]
    
    return(st_P, st_S)
    

def get_takeoff_angle(hyp_ind,dist,rupt_file,model):
    
    from obspy.geodetics import kilometer2degrees
    from obspy.taup import TauPyModel
    
    velmod = TauPyModel(model=model)
    
    dist_in_degs=kilometer2degrees(dist/1000.)
    zs=rupt_file[hyp_ind,3]
    
    #Get ray paths for all direct P arrivals
    try:
        Ppaths=velmod.get_ray_paths(zs,dist_in_degs,phase_list=['P','p'])
    except:
        zs = zs+0.0001
        Ppaths=velmod.get_ray_paths(zs,dist_in_degs,phase_list=['P','p'])
    
    directP=Ppaths[0]
    take_off_angle_P=directP.takeoff_angle
    
    return(take_off_angle_P)
    
