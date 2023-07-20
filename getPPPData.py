#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 10:37:23 2021

script to access all PPP offset data from ONC sites

@author: Angela Schlesinger, Ocean Networks Canada
"""

class PPPdata:
    """ A simple class for retrieving PPP data """
    
    
    
    def __init__(self,Token):
        from onc.onc import ONC
        self.onc = ONC(token=Token)
    
    def createTimeFormat(self,timestamp):
        import pandas as pd
        ts = pd.to_datetime(timestamp)
        self.ts = ts.strftime('%Y-%m-%dT%H:%M:%S'+'.000Z')
        return(self.ts)
        
        

    
    def findPPPdata(self,filters):
        """ Identify locations that have PPP data for a certain time-range"""
        
        return(self.onc.getDeployments(filters))
    
 
            
    def getParsedData(self,filters):
   
        return(self.onc.getDirectByLocation(filters,allPages=True))
    
    

 
if '__main__':
    from getToken import Token
    import pandas as pd
    PPP = PPPdata(Token=Token())

    station = 'TFNO'
    
    dateFrom = '2021-09-05T00:00:00.000Z'
    dateTo='2021-09-08T00:00:00.000Z'

    deviceCategoryCodes=['PPPINT','PPPFLT','PPPORB']
    
    """ (1) get the locations by querying the discovery service """
    d={}
    for DCC in deviceCategoryCodes:
        print(DCC)
        temp=PPP.findPPPdata({'deviceCategoryCode':DCC,
                         'dateFrom':dateFrom,
                         'dateTo':dateTo
                         })
        d[DCC]= temp
        
    """ (2)  extract the locations from currently installed GNSS/PPP sites """
    
    locations = [d['PPPINT'][i]['locationCode'] for i,j in enumerate(d['PPPINT'])]
    
    
    # """ (2b) find the sensors you are interested in:this is a coumbersome step as we
    #     currently have no service to
    #     identify the sensors a user is interested in, hence you need to query 
    #     a few samples of data and return the sensorCategoryCodes that exist"""
         
    # sensors = PPP.getParsedData({'deviceCategoryCode':deviceCategoryCodes[0],
    #                                        'locationCode':locations[0],
    #                                        'dateFrom':PPP.createTimeFormat(dateFrom),
    #                                        'dateTo':PPP.createTimeFormat(dateTo),
    #                                        'rowLimit':2})
    # sensorCategoryCodes = [sensors['sensorData'][i]['sensorCategoryCode'] 
    #                        for i,j in enumerate(sensors['sensorData'])]
    
   
    """  (3) get the Parsed data for the PPP by locaiton and stream
        you are only interested in some data from these devices:
        
        sensorCategoryCodes:['east_offset','north_offset','vertical_offset',
                             'east_offset_stddev','north_offset_stdev',
                             'vertical_offset_stddev']
          
    """
    
    sCodes = ['east_offset',
                           'north_offset',
                           'vertical_offset',
                           'east_offset_stddev',
                           'north_offset_stddev',
                           'vertical_offset_stddev']
    

            
    """ (4) run for one EEW site and one stream (integer)     """
    
    temp = PPP.getParsedData({'locationCode':station,
                        'deviceCategoryCode':'PPPINT',
                        'sensorCategoryCodes':sCodes,
                        'dateFrom':PPP.createTimeFormat(dateFrom),
                        'dateTo':PPP.createTimeFormat(dateTo)#'
                             })
          
    headers = [temp['sensorData'][i]['sensorCategoryCode'] for i,j in enumerate(temp['sensorData'])]
    df = pd.DataFrame(columns=headers)
    
    for header in headers:
        for idx,itemp in enumerate(temp['sensorData']):
            if header == temp['sensorData'][idx]['sensorCategoryCode']:
                df[header]=temp['sensorData'][idx]['data']['values']
     
    
    """ PPlot for quick check of data"""
    df.plot(subplots=True, layout=(len(headers),1))    

    df.to_csv(f'/Users/tnye/ONC/data/station_noise/PPP/{station}.csv')       