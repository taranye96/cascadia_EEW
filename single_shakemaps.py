#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 16:02:14 2021

@author: vjs
"""
## Make shakemap for Cascadia downdip ruptures

###### FUNCTIONS ######
import numpy as np
    
###### INTERNAL FUNCTIONS ######
def convert_rupt2geo(rupt_file,vs30_csv,geo_file,slip_percent=0):
    '''
    Convert a MudPy .rupt file to a .geo file 
    Input:
        rupt_file:      Path to input rupture file
        vs30_csv:       Path to Vs30 CSV file
        slip_percent:   Float with percent of max slip to keep (i.e., 15%)
        geo_file:       Path to output .geo file
    Output:
        geo_file:       Writes to this file, does not return anything  
    '''
    from pyproj import Geod        

    #get rupture
    rupt = np.genfromtxt(rupt_file)
    #keep only those with slip
    i = np.where(rupt[:,12]>0)[0]
    rupt = rupt[i,:]
    
    #get coordaintes and vs30 values
    d = np.genfromtxt(vs30_csv)
    lon = d[:,0]
    lat = d[:,1]
    vs30 = d[:,2]
    
    #get Rrupt
    #projection obnject
    p = Geod(ellps='WGS84')
    
    #lon will have as many rows as Vs30 points and as many columns as subfautls in rupture
    Nsubfaults = len(rupt)
    Nvs30 = len(vs30)
    lon_surface = np.tile(lon,(Nsubfaults,1)).T
    lat_surface = np.tile(lat,(Nsubfaults,1)).T
    lon_subfaults = np.tile(rupt[:,1],(Nvs30,1))-360
    lat_subfaults = np.tile(rupt[:,2],(Nvs30,1))
    az,baz,dist = p.inv(lon_surface,lat_surface,lon_subfaults,lat_subfaults)
    dist = dist/1000
    
    #get 3D distance
    z = np.tile(rupt[:,3],(Nvs30,1))
    xyz_dist = (dist**2 + z**2)**0.5
    Rrupt = xyz_dist.min(axis=1)
    
    #get scalar moment
    M0 = sum(rupt[:,10]*rupt[:,11]*rupt[:,13]*((rupt[:,8]**2+rupt[:,9]**2)**0.5))
    Mw = (2./3)*(np.log10(M0)-9.1)
    
    #write to file
    out = np.c_[lon,lat,vs30,Rrupt,np.ones(len(lon))*Mw]
    np.savetxt(geo_file,out,fmt='%.5f,%.5f,%d,%.2f,%.2f',header='lon,lat,vs30(m/s),Rrupt(km),Mw')
    
####################################
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

####################################
def compute_shaking_BCHydro(rupt_df):
    '''
    For a given rupture model file, compute the PGA with BCHydro and 
        its uncertainty, and the MMI with WGRW12. Returns new dataframe with
        these values attached.
    Input:
        rupt_df:     Pandas dataframe with the following columns/header:
                        # lon,lat,vs30(m/s),Rrupt(km),Mw
    Output:
        shaking_df:  Pandas dataframe with two additional columns:
    '''
        ## Open quake stuff:
    from openquake.hazardlib import imt, const
    from openquake.hazardlib.gsim.base import RuptureContext
    from openquake.hazardlib.gsim.base import DistancesContext
    from openquake.hazardlib.gsim.base import SitesContext
    from openquake.hazardlib.gsim.abrahamson_2018 import AbrahamsonEtAl2018SInter
    import pandas as pd
    import numpy as np
             
    ## Define intensity measure and uncertainty
    im_type = imt.PGA()
    uncertaintytype = const.StdDev.TOTAL
    
    ## Set GMPEs:
    bchydro = AbrahamsonEtAl2018SInter()
    
    ## Set the empty arrays:
    median_bchydro = np.array([])
    sd_bchydro = np.array([])
    
    
    ## Make contexts:
    rctx = RuptureContext()
    dctx = DistancesContext()
    sctx = SitesContext()
    
    ## Use  dataframe distance and vs30 in dist and site context as arrays:
    dctx.rrup = rupt_df['Rrupt(km)'].values
    
    ## Magnitude in rupture context is always the same, so grab first one:
    rctx.mag = rupt_df['Mw'].values[0]
    
    ## Site context - seems to now need to be from a "site collection", which seems to be a pandas dataframe.
    ##   Boore 2020 needs vs30, and nothing else - so set them for Nans. These assume the function is used on a single
    ###    value to predict - if you're doing it in array form (many records for a single earthquake), then make these arrays.
    sitecol_dict = {'sids':np.arange(1,len(rupt_df['vs30(m/s)'])+1,1),'vs30':rupt_df['vs30(m/s)'].values,'vs30measured':np.full_like(rupt_df['vs30(m/s)'].values,np.nan),'z1pt0':np.full_like(rupt_df['vs30(m/s)'].values,np.nan),'z2pt5':np.full_like(rupt_df['vs30(m/s)'].values,np.nan)}
    sitecollection = pd.DataFrame(sitecol_dict)
    
    ## Then put into a sites context:
    sctx = SitesContext(sitecol=sitecollection)    

    
    ln_median_bchydro,sd_bchydro = bchydro.get_mean_and_stddevs(sctx, rctx, dctx, im_type, [uncertaintytype])
    
    ## Convert median from ln g to g:
    median_bchydro = np.exp(ln_median_bchydro)
    
    ## Get it in cm/s/s for MMI conversion:
    median_bchydro_cm_s2 = (median_bchydro*9.81)*100
    MMI_bchydro_pga = WGRW12(median_bchydro_cm_s2,0)
    
    ## Add as column in dataframe:
    shaking_df = rupt_df.copy()
    shaking_df['PGA_bchydro2018_g'] = median_bchydro
    shaking_df['PGAsd_bchydro2018_lng'] = sd_bchydro[0]  ## This is in a list
    shaking_df['MMI_bchydro_wgrw12_pga'] = MMI_bchydro_pga
    
    return shaking_df


####################################
def compute_shaking_NGA_Sub(rupt_df):
    '''
    For a given rupture model file, compute the PGA with NGA-Sub and 
        its uncertainty, and the MMI with WGRW12. Returns new dataframe with
        these values attached.
    Input:
        rupt_df:     Pandas dataframe with the following columns/header:
                        # lon,lat,vs30(m/s),Rrupt(km),Mw
    Output:
        shaking_df:  Pandas dataframe with two additional columns:
    '''
    ## Open quake stuff:
    from openquake.hazardlib import imt, const
    from openquake.hazardlib.gsim.base import RuptureContext
    from openquake.hazardlib.gsim.base import DistancesContext
    from openquake.hazardlib.gsim.base import SitesContext
    from openquake.hazardlib.gsim.parker_2020 import ParkerEtAl2020SInter
    import pandas as pd
    import numpy as np
             
    ## Define intensity measure and uncertainty
    IMs = [imt.PGA(), imt.PGV()]
    shaking_df = rupt_df.copy()
    for im_type in IMs:
        
        uncertaintytype = const.StdDev.TOTAL
        
        ## Set GMPEs:
        NGASub = ParkerEtAl2020SInter()
        
        ## Set the empty arrays:
        median_NGASub = np.array([])
        sd_NGASub = np.array([])
        
        
        ## Make contexts:
        rctx = RuptureContext()
        dctx = DistancesContext()
        sctx = SitesContext()
        
        ## Use  dataframe distance and vs30 in dist and site context as arrays:
        dctx.rrup = rupt_df['Rrupt(km)'].values
        
        ## Magnitude in rupture context is always the same, so grab first one:
        rctx.mag = rupt_df['Mw'].values[0]
        
        ## Site context - seems to now need to be from a "site collection", which seems to be a pandas dataframe.
        ##   Boore 2020 needs vs30, and nothing else - so set them for Nans. These assume the function is used on a single
        ###    value to predict - if you're doing it in array form (many records for a single earthquake), then make these arrays.
        sitecol_dict = {'sids':np.arange(1,len(rupt_df['vs30(m/s)'])+1,1),'vs30':rupt_df['vs30(m/s)'].values,'vs30measured':np.full_like(rupt_df['vs30(m/s)'].values,np.nan),'z1pt0':np.full_like(rupt_df['vs30(m/s)'].values,np.nan),'z2pt5':np.full_like(rupt_df['vs30(m/s)'].values,np.nan)}
        sitecollection = pd.DataFrame(sitecol_dict)
        
        ## Then put into a sites context:
        sctx = SitesContext(sitecol=sitecollection)    
    
        
        ln_median_nga_sub,sd_nga_sub = NGASub.get_mean_and_stddevs(sctx, rctx, dctx, im_type, [uncertaintytype])
        
        ## Convert median from ln g to g:
        median_NGASub = np.exp(ln_median_nga_sub)
    
        ## Get it in cm/s/s for MMI conversion:
        if im_type == imt.PGA():
            median_NGASub_cm_s2 = (median_NGASub*9.81)*100
            MMI_NGASub_pga = WGRW12(median_NGASub_cm_s2,0)
        
            ## Add as column in dataframe:
            shaking_df['PGA_NGASub_g'] = median_NGASub
            shaking_df['PGAsd_NGASub_lng'] = sd_nga_sub[0]  ## This is in a list
            shaking_df['MMI_NGASub_wgrw12_pga'] = MMI_NGASub_pga
        else:
            MMI_NGASub_pgv = WGRW12(median_NGASub,1)
            shaking_df['PGV_NGASub_cm/s'] = median_NGASub
            shaking_df['PGVsd_NGASub_lncm/s'] = sd_nga_sub[0]
            shaking_df['MMI_NGASub_wgrw12_pgv'] = MMI_NGASub_pgv
    
    return shaking_df


def make_shake_map(geo_file,output_file):
    '''
    Make a Shakemap file 
    '''

    import pandas as pd
    from shlex import split
    
    #####################################

    
    ## Import the file:
    rupt_df = pd.read_csv(geo_file)
    
    ## Get shaking dataframe:
    bc_df = compute_shaking_BCHydro(rupt_df)
    nga_df = compute_shaking_NGA_Sub(rupt_df)
    
    output_df = pd.concat([bc_df,nga_df.iloc[:,5:]], axis=1)
    
    ## Write output dataframe to specified file:
    output_df.to_csv(output_file,index=False)
    
    
###############################################################################    

def run_shakemaps(batch, rupt_file, vs30_csv, stn_type):

    import os
    
    # batch = rupt_file.split('/')[-3]
    run = rupt_file.split('/')[-1].strip('.rupt')
    
    ## Output directories:
    geo_directory = f'/Users/tnye/ONC/shakemaps/geofiles/{batch}' ### Note no / at end
    shake_directory = f'/Users/tnye/ONC/shakemaps/shakefiles/{batch}' ### Note no / at end
    # geo_directory = f'/Users/tnye/ONC/shakemaps/geofiles_attenuation_test/{batch}' ### Note no / at end
    # shake_directory = f'/Users/tnye/ONC/shakemaps/shakefiles_attenuation_test/{batch}' ### Note no / at end
    
    ## Output file paths:
    geo_file = f'{geo_directory}/{stn_type}/{run}.geo'
    shake_file = f'{shake_directory}/{stn_type}/{run}.shake'
    
    ## Slip percent:
    slippercent = 15
    
    # Make sure output directories exist
    if not os.path.exists(geo_file.strip(f'{run}.geo')):
        os.makedirs(geo_file.strip(f'{run}.geo'))
    if not os.path.exists(shake_file.strip(f'{run}.shake')):
        os.makedirs(shake_file.strip(f'{run}.shake'))
    
    
    #########
    ## Make geo file:
    convert_rupt2geo(rupt_file,vs30_csv,geo_file,slip_percent=slippercent)
        
    #########
    ## Make shake file
    make_shake_map(geo_file,shake_file)

                    