# Name:        mapper_obpg_l2
# Purpose:     Mapping for L2 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from vrt import GeolocationArray, VRT, gdal, osr, latlongSRS
from datetime import datetime, timedelta
from math import ceil

class Mapper(VRT):
    ''' Mapper for SeaWIFS/MODIS/MERIS/VIIRS L2 data from OBPG
    
    TODO:
    * Test on SeaWIFS
    * Test on MODIS Terra
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        ''' Create VRT '''
        # number of GCPs along each dimention
        GCP_COUNT = 10
        """
        Title=
        Title=HMODISA Level-2 Data
        Title=MERIS Level-2 Data
        """
        
        titles = ['HMODISA Level-2 Data', 'MERIS Level-2 Data']
        # should raise error in case of not obpg_l2 file
        title = gdalMetadata["Title"]
        assert title in titles, 'obpg_l2 BAD MAPPER'

        # get subdataset and parse to VRT.__init__() for retrieving geo-metadata
        # but NOT from longitude or latitude because it can be smaller!
        subDatasets = gdalDataset.GetSubDatasets()
        for subDataset in subDatasets:
            if 'longitude' not in subDataset[1] and 'latitude' not in subDataset[1]:
                gdalSubDataset = gdal.Open(subDataset[0])
                break
        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalSubDataset)
        
        # parts of dictionary for all Reflectances
        #dictRrs = {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '412'} }
        # dictionary for all possible bands
        allBandsDict = {
        'Rrs':        {'src': {}, 'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air'}},
        'Kd':         {'src': {}, 'dst': {'wkv': 'volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water'}},
        'chlor_a':    {'src': {}, 'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'case': 'I'}},
        'cdom_index': {'src': {}, 'dst': {'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter', 'case': 'II'}},
        'sst':        {'src': {}, 'dst': {'wkv': 'sea_surface_temperature'}},
        'l2_flags':   {'src': {'SourceType': 'SimpleSource', 'DataType': 4},
                       'dst': {'wkv': 'quality_flags', 'dataType': 4}},
        'qual_sst':   {'src': {'SourceType': 'SimpleSource', 'DataType': 4},
                       'dst': {'wkv': 'quality_flags', 'name': 'qual_sst', 'dataType': 4}},
        'latitude':   {'src': {}, 'dst': {'wkv': 'latitude'}},
        'longitude':  {'src': {}, 'dst': {'wkv': 'longitude'}},
        }
        
        # loop through available bands and generate metaDict (non fixed)
        metaDict = []
        bandNo = 0
        for subDataset in subDatasets:
            # get sub dataset name
            subDatasetName = subDataset[1].split(' ')[1]
            self.logger.debug('Subdataset: %s' % subDataset[1])
            self.logger.debug('Subdataset name: "%s"' % subDatasetName)
            # get wavelength if applicable, get dataset name without wavelength
            try:
                wavelength = int(subDatasetName.split('_')[-1])
            except:
                wavelength = None
                subBandName = subDatasetName
            else:
                subBandName = subDatasetName.split('_')[0]

            self.logger.debug('subBandName, wavelength: %s %s' % (subBandName, str(wavelength)))
            
            if subBandName in allBandsDict:
                # get name, slope, intercept
                self.logger.debug('name: %s' % subBandName)
                tmpSubDataset = gdal.Open(subDataset[0])
                tmpSubMetadata = tmpSubDataset.GetMetadata()
                slope = tmpSubMetadata.get('slope', '1')
                intercept = tmpSubMetadata.get('intercept', '0')
                self.logger.debug('slope, intercept: %s %s ' % (slope, intercept))
                # create meta entry
                metaEntry = {'src': {'SourceFilename': subDataset[0],
                                     'sourceBand':  1,
                                     'ScaleRatio': slope,
                                     'ScaleOffset': intercept},
                              'dst': {}}
                # add more to src
                for srcKey in allBandsDict[subBandName]['src']:
                    metaEntry['src'][srcKey] = allBandsDict[subBandName]['src'][srcKey]
                # add dst from allBandsDict
                for dstKey in allBandsDict[subBandName]['dst']:
                    metaEntry['dst'][dstKey] = allBandsDict[subBandName]['dst'][dstKey]
                
                # add wavelength, band name to dst
                if wavelength is not None:
                    metaEntry['dst']['suffix'] = str(wavelength)
                    metaEntry['dst']['wavelength'] = str(wavelength)

                # append band metadata to metaDict
                self.logger.debug('metaEntry: %d => %s' % (bandNo, str(metaEntry)))
                metaDict.append(metaEntry)
                bandNo += 1
                
                if subBandName == 'Rrs':
                    metaEntryRrsw = {
                        'src': [{
                            'SourceFilename': subDataset[0],
                            'SourceBand':  1,
                            'ScaleRatio': slope,
                            'ScaleOffset': intercept,
                            'DataType': 6}],
                        'dst': {
                            'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_water',
                            'suffix': str(wavelength),
                            'wavelength': str(wavelength),
                            'PixelFunctionType': 'NormReflectanceToRemSensReflectance',
                        }}

                    # append band metadata to metaDict
                    self.logger.debug('metaEntry: %d => %s' % (bandNo, str(metaEntryRrsw)))
                    metaDict.append(metaEntryRrsw)
                    bandNo += 1
                    
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        geolocationMetadata = gdalSubDataset.GetMetadata('GEOLOCATION')
        xDatasetSource = geolocationMetadata['X_DATASET']
        xDatasetBand = geolocationMetadata['X_BAND']
        xDataset = gdal.Open(xDatasetSource)
        xVRT = VRT(vrtDataset=xDataset)

        yDatasetSource = geolocationMetadata['Y_DATASET']
        yDatasetBand = geolocationMetadata['Y_BAND']
        yDataset = gdal.Open(yDatasetSource)
        yVRT = VRT(vrtDataset=yDataset)
        
        # estimate pixel/line step
        pixelStep = int(ceil(float(gdalSubDataset.RasterXSize) / float(xDataset.RasterXSize)))
        lineStep = int(ceil(float(gdalSubDataset.RasterYSize) / float(xDataset.RasterYSize)))
        self.logger.debug('pixel/lineStep %f %f' % (pixelStep, lineStep))
        
        # ==== ADD GCPs and Pojection ====        
        # get lat/lon matrices
        longitude = xVRT.dataset.GetRasterBand(1).ReadAsArray()
        latitude = yVRT.dataset.GetRasterBand(1).ReadAsArray()

        # add geolocation array
        if (    longitude.min() >= -180 and longitude.max() <= 180 and
                latitude.min() >= -90 and latitude.max() <= 90):
            self.add_geolocationArray(GeolocationArray(xDatasetSource, yDatasetSource, pixelStep=pixelStep, lineStep=lineStep))
        else:
            self.remove_geolocationArray()

        # estimate step of GCPs
        step0 = max(1, int(float(latitude.shape[0]) / GCP_COUNT))
        step1 = max(1, int(float(latitude.shape[1]) / GCP_COUNT))
        self.logger.debug('gcpCount: %d %d %f %d %d', latitude.shape[0], latitude.shape[1], GCP_COUNT, step0, step1)
        
        # generate list of GCPs
        gcps = []
        k = 0
        for i0 in range(0, latitude.shape[0], step0):
            for i1 in range(0, latitude.shape[1], step1):
                # create GCP with X,Y,pixel,line from lat/lon matrices
                lon = float(longitude[i0, i1])
                lat = float(latitude[i0, i1])
                if (lon >= -180 and lon <= 180 and lat >= -90 and lat <= 90):
                    gcp = gdal.GCP(lon, lat, 0, i1 * pixelStep, i0 * lineStep)
                    self.logger.debug('%d %d %d %f %f', k, gcp.GCPPixel, gcp.GCPLine, gcp.GCPX, gcp.GCPY)
                    gcps.append(gcp)
                    k += 1
        
        # append GCPs and lat/lon projection to the vsiDataset
        self.dataset.SetGCPs(gcps, latlongSRS.ExportToWkt())
        
        # set TIME
        startYear = int(gdalMetadata['Start Year'])
        startDay =  int(gdalMetadata['Start Day'])
        startMillisec =  int(gdalMetadata['Start Millisec'])
        startDate = datetime(startYear, 1, 1) + timedelta(startDay, 0, 0, startMillisec)
        self._set_time(startDate)
