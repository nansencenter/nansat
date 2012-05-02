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
from vrt import VRT, gdal, osr
import re
from matplotlib import pyplot as plt
import numpy as np

class Mapper(VRT):
    ''' VRT with mapping of WKV for MODIS Level 1 (QKM, HKM, 1KM) '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30):
        ''' Create MODIS_L1 VRT '''
        GCP_COUNT = 1000
        # get 1st subdataset and parse to VRT.__init__() for retrieving geo-metadata
        subDatasets = gdalDataset.GetSubDatasets()
       
        # should raise error in case of not MODIS_L2_NRT
        mTitle = gdalMetadata["Title"];
        if mTitle is not 'HMODISA Level-2 Data':
            AttributeError("MODIS_L2_NRT BAD MAPPER");

        # create empty VRT dataset with geolocation only
        gdalSubDataset = gdal.Open(subDatasets[0][0])
        VRT.__init__(self, gdalSubDataset, logLevel=logLevel);
        
        # parts of dictionary for all Reflectances
        #dictRrs = {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'parameters': {'wavelength': '412'} }
        # dictionary for all possible bands
        allBandsDict = {
        'Kd_490': {'wkv': 'volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water', 'parameters': {'wavelength': '490'} },
        'chlor_a': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'parameters': {'band_name': 'algal_1', 'case': 'I'} },
        'cdom_index': {'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter', 'parameters': {'band_name': 'yellow_subs', 'case': 'II'} },
        'l2_flags': {'wkv': 'quality_flags', 'parameters': {'band_name': 'l2_flags'} },
        'latitude': {'wkv': 'latitude', 'parameters': {'band_name': 'latitude'} },
        'longitude': {'wkv': 'longitude', 'parameters': {'band_name': 'longitude'} },        
        }
        
        # loop through available bands and generate metaDict (non fixed)
        metaDict = []
        for subDataset in subDatasets:
            self.logger.debug('Subdataset: %s' % subDataset[1])
            # try to get Rrs_412 or similar from subdataset name
            # if success - append Reflectance with respective parameters to meta
            rrsBandName = re.findall('Rrs_\d*', subDataset[1])
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
                    if bandName in subDataset[1]:
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
            self.logger.debug(metaEntry)
            if metaEntry is not None:
                metaDict.append(metaEntry)
        self.logger.debug(metaDict)
        
        self.logger.debug('metaDict: %s' % metaDict)
        
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # ==== add GCPs and Pojection ====
        for subDataset in subDatasets:
            # read longitude matrix
            if 'longitude' in subDataset[1]:
                ds = gdal.Open(subDataset[0])
                b = ds.GetRasterBand(1)
                longitude = b.ReadAsArray()
            # read latitude matrix
            if 'latitude' in subDataset[1]:
                ds = gdal.Open(subDataset[0])
                b = ds.GetRasterBand(1)
                latitude = b.ReadAsArray()
        self.logger.debug('Lat/Lon grids read')
        
        # estimate step of GCPs
        gcpSize = np.sqrt(GCP_COUNT)
        step0 = max(1, int(float(latitude.shape[0]) / gcpSize))
        step1 = max(1, int(float(latitude.shape[1]) / gcpSize))
        self.logger.debug('gcpCount: %d %d %f %d %d', latitude.shape[0], latitude.shape[1], gcpSize, step0, step1)
        
        # generate list of GCPs
        gcps = []
        k = 0
        for i0 in range(0, latitude.shape[0], step0):
            for i1 in range(0, latitude.shape[1], step1):
                #self.logger.debug('%d %d %f %f', i0, i1, longitude[i0, i1], latitude[i0, i1])
                # create GCP with X,Y,pixel,line from lat/lon matrices
                gcp = gdal.GCP(float(longitude[i0, i1]),
                               float(latitude[i0, i1]),
                               0,
                               i1,
                               i0)
                self.logger.debug('%d %d %d %f %f', k, gcp.GCPPixel, gcp.GCPLine, gcp.GCPX, gcp.GCPY)
                gcps.append(gcp)
                k += 1
        
        # append GCPs and lat/lon projection to the vsiDataset
        latlongSRS = osr.SpatialReference()
        latlongSRS.ImportFromProj4("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs")
        latlongSRSWKT = latlongSRS.ExportToWkt()
        self.dataset.SetGCPs(gcps, latlongSRSWKT)

        # add GEOLOCATION
        for subDataset in subDatasets:
            if 'longitude' in subDataset[1]:
                xSubDatasetName = subDataset[0]
            if 'latitude' in subDataset[1]:
                ySubDatasetName = subDataset[0]
        
        latlongSRS = osr.SpatialReference()
        latlongSRS.ImportFromProj4("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs")
        latlongSRSWKT = latlongSRS.ExportToWkt()
        self.dataset.SetProjection(latlongSRSWKT)
        self.dataset.SetMetadataItem('LINE_OFFSET', '0', 'GEOLOCATION')
        self.dataset.SetMetadataItem('LINE_STEP', '1', 'GEOLOCATION')
        self.dataset.SetMetadataItem('PIXEL_OFFSET', '0', 'GEOLOCATION')
        self.dataset.SetMetadataItem('PIXEL_STEP', '1', 'GEOLOCATION')
        self.dataset.SetMetadataItem('SRS', latlongSRSWKT, 'GEOLOCATION')
        self.dataset.SetMetadataItem('X_BAND', '1', 'GEOLOCATION')
        self.dataset.SetMetadataItem('X_DATASET', xSubDatasetName, 'GEOLOCATION')
        self.dataset.SetMetadataItem('Y_BAND', '1', 'GEOLOCATION')
        self.dataset.SetMetadataItem('Y_DATASET', ySubDatasetName, 'GEOLOCATION')
        
        return
