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
from vrt import Geolocation, VRT, gdal, osr, os
import re
from datetime import datetime, timedelta

import numpy as np
from scipy.io import loadmat


class Mapper(VRT):
    ''' MApper for Matlab files with SMOS data '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        ''' Create SMOS VRT '''
        # check extension
        fName = os.path.split(fileName)[1]
        fExt = os.path.splitext(fileName)[1]
        if fExt == '.MAT' or fExt == '.mat' and 'OSUDP2' in fName:
            # load file
            matFile = loadmat(fileName)
        else:
            AttributeError("SMOS BAD MAPPER");
            
        # get geolocation
        geolocArray = matFile['geolocation'][0]
        srcProj4 = '+proj=stere +lon_0=%f +lat_0=%f +datum=WGS84 +ellps=WGS84 +units=km +no_defs' % (geolocArray[0], geolocArray[1])
        srcProjection = osr.SpatialReference()
        srcProjection.ImportFromProj4(srcProj4)
        srcProjection = srcProjection.ExportToWkt()
        srcGeotransform = (geolocArray[2], geolocArray[4], 0, geolocArray[3], 0, geolocArray[5])
        lon = matFile['longitude']
        #lat = matFile['latitude']
        srcRasterYSize, srcRasterXSize = lon.shape
        # create VRT from lat/lon
        #VRT.__init__(self, lon=lon, lat=lat)
        VRT.__init__(self, srcGeoTransform=srcGeotransform,
                            srcProjection=srcProjection,
                            srcRasterXSize=srcRasterXSize,
                            srcRasterYSize=srcRasterYSize)
        
        # add the following variables
        varNames = ['SSS1', 'SSS2', 'SSS3', 'SST',
                    'Sigma_SSS1', 'Sigma_SSS2', 'Sigma_SSS3', 
                    'Control_Flags_1', 'Control_Flags_2', 'Control_Flags_3', 'Control_Flags_4', 
                    'Science_Flags_1', 'Science_Flags_2', 'Science_Flags_3', 'Science_Flags_4']
        self.varVRTs = {}
        metaDict = []
        for varName in varNames:
            var = matFile[varName]
            self.varVRTs[varName] = VRT(array=var)
            metaDict.append(
            {'src': {'SourceFilename': self.varVRTs[varName].fileName, 'sourceBand':  1},
             'dst': {'BandName': varName}})
        
        # create mask
        cloudBits = [2, 3, 4, 5, 6]
        maxSigma = 3.0
        mask = np.zeros(lon.shape, 'uint16')
        mask[:] = 128
        mask[np.isnan(matFile['SSS1'])] = 0
        mask[matFile['Sigma_SSS1'] > maxSigma] = 1
        mask[matFile['Sigma_SSS2'] > maxSigma] = 1
        mask[matFile['Sigma_SSS3'] > maxSigma] = 1
        for cloudBit in cloudBits:
            for cfi in range(1, 5):
                bitMap = np.bitwise_and(matFile['Control_Flags_%d' % cfi], np.power(2, cloudBit))
                mask[bitMap > 0] = 1

        self.varVRTs['mask'] = VRT(array=mask)
        metaDict.append(
        {'src': {'SourceFilename': self.varVRTs['mask'].fileName, 'sourceBand':  1},
         'dst': {'BandName': 'mask'}})


        self.logger.debug('metaDict: %s' % metaDict)
        
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
        
        # set TIME
            
        """
        
        # parts of dictionary for all Reflectances
        #dictRrs = {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air', 'wavelength': '412'} }
        # dictionary for all possible bands
        allBandsDict = {
        'Kd_490': {'wkv': 'volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water', 'BandName': 'Kd_490', 'wavelength': '490'},
        'chlor_a': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water', 'BandName': 'algal_1', 'case': 'I'},
        'cdom_index': {'wkv': 'volume_absorption_coefficient_of_radiative_flux_in_sea_water_due_to_dissolved_organic_matter', 'BandName': 'yellow_subs', 'case': 'II'},
        'sst': {'wkv': 'sea_surface_temperature', 'BandName': 'sst'},
        'l2_flags': {'wkv': 'quality_flags', 'SourceType': 'SimpleSource', 'BandName': 'l2_flags'},
        'latitude': {'wkv': 'latitude', 'BandName': 'latitude'},
        'longitude': {'wkv': 'longitude', 'BandName': 'longitude'},
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
            metaEntry2 = None
            if len(rrsBandName) > 0:
                tmpSubDataset = gdal.Open(subDataset[0])
                slope = tmpSubDataset.GetMetadataItem('slope')
                intercept = tmpSubDataset.GetMetadataItem('intercept')
                metaEntry = {'src': {'SourceFilename': subDataset[0], 'sourceBand':  1, 'ScaleRatio': slope, 'ScaleOffset': intercept},
                             'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air',
                                     'BandName': rrsBandName[0],
                                     'wavelength': rrsBandName[0][-3:],
                                      }
                             }
                metaEntry2 = {'src': {'SourceFilename': subDataset[0], 'SourceBand':  1},
                              'dst': {'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_water',
                                      'BandName': rrsBandName[0].replace('Rrs', 'Rrsw'),
                                      'wavelength': rrsBandName[0][-3:],
                                      'expression': 'self["%s"] / (0.52 + 1.7 * self["%s"])' % (rrsBandName[0], rrsBandName[0]),
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
                        metaEntry = {'src': {'SourceFilename': subDataset[0], 'SourceBand':  1, 'ScaleRatio': slope, 'ScaleOffset': intercept},
                                     'dst': allBandsDict[bandName]}
                        if 'SourceType' in allBandsDict[bandName]:
                            metaEntry['src']['SourceType'] = allBandsDict[bandName]['SourceType']
                
                self.logger.debug('metaEntry %s' % metaEntry)
                    
            if metaEntry is not None:
                metaDict.append(metaEntry)
            if metaEntry2 is not None:
                metaDict.append(metaEntry2)
        
        self.logger.debug('metaDict: %s' % metaDict)
        
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # make lon/lat VRT objects
        for subDataset in subDatasets:
            if 'longitude' in subDataset[1]:
                xDatasetSource = subDataset[0]
                xDataset = gdal.Open(xDatasetSource)
                xVRT = VRT(vrtDataset=xDataset)
            if 'latitude' in subDataset[1]:
                yDatasetSource = subDataset[0]
                yVRT = VRT(vrtDataset=gdal.Open(yDatasetSource))

        # estimate pixel/line step
        pixelStep = int(float(gdalSubDataset.RasterXSize) / float(xDataset.RasterXSize))
        lineStep = int(float(gdalSubDataset.RasterYSize) / float(xDataset.RasterYSize))
        
        # add geolocation
        self.add_geolocation(Geolocation(xDatasetSource, yDatasetSource, pixelStep=pixelStep, lineStep=lineStep))

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
                               0, i1 * pixelStep, i0 * lineStep)
                self.logger.debug('%d %d %d %f %f', k, gcp.GCPPixel, gcp.GCPLine, gcp.GCPX, gcp.GCPY)
                gcps.append(gcp)
                k += 1
        
        # append GCPs and lat/lon projection to the vsiDataset
        latlongSRS = osr.SpatialReference()
        latlongSRS.ImportFromProj4("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs")
        latlongSRSWKT = latlongSRS.ExportToWkt()
        self.dataset.SetGCPs(gcps, latlongSRSWKT)
        
        # set TIME
        startYear = int(gdalMetadata['Start Year'])
        startDay =  int(gdalMetadata['Start Day'])
        startMillisec =  int(gdalMetadata['Start Millisec'])
        startDate = datetime(startYear, 1, 1) + timedelta(startDay, 0, 0, startMillisec)
        self._set_time(startDate)
        """

