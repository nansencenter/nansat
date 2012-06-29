#-------------------------------------------------------------------------------
# Name:        mapper_geostationary.py
# Purpose:     Generic mappier for all geostationary satellites in Eumetcast format
#
# Author:      Knut-Frode
#
# Created:     27.06.2012
# Copyright:   
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import *

satDict = [\
           {'name': 'GOES13', 'wavelengths': [700, 3900, 6600, 10700]}, 
           {'name': 'GOES15', 'wavelengths': [700, 3900, 6600, 10700]},
           {'name': 'MTSAT2', 'wavelengths': [700, 3800, 6800, 10800]}, # last ch missing 
           {'name': 'MET7', 'wavelengths': [795, 6400, 11500]},
           {'name': 'MSG2', 'wavelengths': [600, 800, 1600, 3900, 6200, 7300, 8700, 9700, 10800, 12000, 13400]}
           ];
        

class Mapper(VRT):
    ''' VRT with mapping of WKV for Geostationary satellite data '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30):
        satellite = gdalDataset.GetDescription().split(",")[2]
        
        for sat in satDict:
            if sat['name'] == satellite:
                print 'This is ' + satellite
                wavelengths = sat['wavelengths']

        if wavelengths is None:
            raise AttributeError("No Eumetcast geostationary satellite");
        
        path = gdalDataset.GetDescription().split(",")[0].split("(")[1]
        datestamp = gdalDataset.GetDescription().split(",")[3]
        if satellite == 'MSG2':
            resolution = 'H'
        else:
            resolution = 'L'

        metaDict = []
        for i, wavelength in enumerate(wavelengths):
            print wavelength
            if wavelength > 1000:
                standard_name = 'brightness_temperature'
            else:
                standard_name = 'grayscale'
            metaDict.append(
                {'source': 'MSG('+path+','+resolution+','+satellite+','+datestamp+','+str(i+1)+',Y,F,1,1)',
                 'sourceBand': 1, 'wkv': standard_name, 
                 'parameters': {'wavelength': str(wavelength)}})
                      
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, logLevel=logLevel);

        # Create bands
        self._create_bands(metaDict)

        print "successful"
        return
