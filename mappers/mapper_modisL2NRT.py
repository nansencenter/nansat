#-------------------------------------------------------------------------------
# Name:        mapper_modisL1
# Purpose:     Mapping for MODIS-L1 data
#
# Author:      antonk
#
# Created:     13.12.2011
# Copyright:   (c) NERSC 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import Geolocation, VRT, gdal, osr
import re
from datetime import datetime, timedelta

class Mapper(VRT):
    ''' VRT with mapping of WKV for MODIS Level 1 (QKM, HKM, 1KM) '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        ''' Create MODIS_L2 VRT '''
        # number of GCPs along each dimention
        GCP_COUNT = 10
       
        # should raise error in case of not MODIS_L2_NRT
        mTitle = gdalMetadata["Title"];
        if mTitle is not 'HMODISA Level-2 Data':
            AttributeError("MODIS_L2_NRT BAD MAPPER");

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
        #dictRrs = {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '412'} }
        # dictionary for all possible bands
        allBandsDict = {
        'Kd_490': {'wkv': 'volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water', 'parameters': {'wavelength': '490'} },
        'chlor_a': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'parameters': {'band_name': 'algal_1', 'case': 'I'} },
        'cdom_index': {'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter', 'parameters': {'band_name': 'yellow_subs', 'case': 'II'} },
        'sst': {'wkv': 'sea_surface_temperature', 'parameters': {'band_name': 'sst'} },        
        'l2_flags': {'wkv': 'quality_flags', 'SourceType':'SimpleSource', 'parameters': {'band_name': 'l2_flags'} },
        'latitude': {'wkv': 'latitude', 'parameters': {'band_name': 'latitude'} },
        'longitude': {'wkv': 'longitude', 'parameters': {'band_name': 'longitude'} },
        }
        
        # loop through available bands and generate metaDict (non fixed)
        metaDict = []
        for subDataset in subDatasets:
            subDatasetName = subDataset[1].split(' ')[1]
            self.logger.debug('Subdataset: %s' % subDataset[1])
            self.logger.debug('Subdataset name: %s' % subDatasetName)
            # try to get Rrs_412 or similar from subdataset name
            # if success - append Reflectance with respective parameters to meta
            rrsBandName = re.findall('Rrs_\d*', subDatasetName)
            metaEntry = None
            if len(rrsBandName) > 0:
                tmpSubDataset = gdal.Open(subDataset[0])
                slope = tmpSubDataset.GetMetadataItem('slope')
                intercept = tmpSubDataset.GetMetadataItem('intercept')
                metaEntry = {'source': subDataset[0],
                                 'sourceBand':  1,
                                 'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                                 'parameters': {'band_name': rrsBandName[0],
                                                'wavelength': rrsBandName[0][-3:],
                                                'scale': slope,
                                                'offset': intercept,
                                                }
                                }
            else:
                # if the band is not Rrs_NNN
                # try to find it (and metadata) in dict of known bands (allBandsDict)
                for bandName in allBandsDict:
                    if bandName == subDatasetName:
                        tmpSubDataset = gdal.Open(subDataset[0])
                        slope = tmpSubDataset.GetMetadataItem('slope')
                        intercept = tmpSubDataset.GetMetadataItem('intercept')
                        metaEntry = {'source': subDataset[0],
                                 'sourceBand':  1,
                                 'wkv': allBandsDict[bandName]['wkv'],
                                 'scale': slope,
                                 'offset': intercept,
                                 'parameters': allBandsDict[bandName]['parameters']
                                }
                        metaEntry['parameters']['scale'] = slope
                        metaEntry['parameters']['offset'] = intercept
                        if 'SourceType' in allBandsDict[bandName]:
                            metaEntry['SourceType'] = allBandsDict[bandName]['SourceType']
                
                self.logger.debug('metaEntry %s' % metaEntry)
                    
            if metaEntry is not None:
                metaDict.append(metaEntry)
        
        self.logger.debug('metaDict: %s' % metaDict)
        
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # make lon/lat VRT objects
        for subDataset in subDatasets:
            if 'longitude' in subDataset[1]:
                xDatasetSource = subDataset[0]
                xVRT = VRT(vrtDataset=gdal.Open(xDatasetSource))
            if 'latitude' in subDataset[1]:
                yDatasetSource = subDataset[0]
                yVRT = VRT(vrtDataset=gdal.Open(yDatasetSource))
        
        # add geolocation
        self.add_geolocation(Geolocation(xDatasetSource, yDatasetSource))

        # ==== ADD GCPs and Pojection ====        
        # get lat/lon matrices
        longitude = xVRT.dataset.GetRasterBand(1).ReadAsArray()
        latitude = yVRT.dataset.GetRasterBand(1).ReadAsArray()

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
                gcp = gdal.GCP(float(longitude[i0, i1]),
                               float(latitude[i0, i1]),
                               0, i1, i0)
                self.logger.debug('%d %d %d %f %f', k, gcp.GCPPixel, gcp.GCPLine, gcp.GCPX, gcp.GCPY)
                gcps.append(gcp)
                k += 1
        
        # append GCPs and lat/lon projection to the vsiDataset
        latlongSRS = osr.SpatialReference()
        latlongSRS.ImportFromProj4("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs")
        latlongSRSWKT = latlongSRS.ExportToWkt()
        self.dataset.SetGCPs(gcps, latlongSRSWKT)
        
        # === ADD GEOLOCATION  ===
        #self.logger.debug('x/y datasets: %s %s', self.lonVRT.fileName, self.latVRT.fileName)
        #self.add_geolocation(self.lonVRT.fileName, 1, self.latVRT.fileName, 1)
        
        # set TIME
        startYear = int(gdalMetadata['Start Year'])
        startDay =  int(gdalMetadata['Start Day'])
        startMillisec =  int(gdalMetadata['Start Millisec'])
        startDate = datetime(startYear, 1, 1) + timedelta(startDay, 0, 0, startMillisec)
        self._set_time(startDate)

