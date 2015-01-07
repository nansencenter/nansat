# Name:        mapper_obpg_l2
# Purpose:     Mapping for L2 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
from datetime import datetime, timedelta
from math import ceil

from nansat.tools import gdal, ogr, WrongMapperError
from nansat.vrt import VRT
from nansat.nsr import NSR

# 2013 x 243 datasets
'''
  SUBDATASET_1_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Area_Mean_Height
  SUBDATASET_1_DESC=[2013x243] //Area_Mean_Height (16-bit integer)
  SUBDATASET_7_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,10.7GHz,H)
  SUBDATASET_7_DESC=[2013x243] //Brightness_Temperature_(res06,10.7GHz,H) (16-bit unsigned integer)
  SUBDATASET_8_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,10.7GHz,V)
  SUBDATASET_8_DESC=[2013x243] //Brightness_Temperature_(res06,10.7GHz,V) (16-bit unsigned integer)
  SUBDATASET_9_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,18.7GHz,H)
  SUBDATASET_9_DESC=[2013x243] //Brightness_Temperature_(res06,18.7GHz,H) (16-bit unsigned integer)
  SUBDATASET_10_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,18.7GHz,V)
  SUBDATASET_10_DESC=[2013x243] //Brightness_Temperature_(res06,18.7GHz,V) (16-bit unsigned integer)
  SUBDATASET_11_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,23.8GHz,H)
  SUBDATASET_11_DESC=[2013x243] //Brightness_Temperature_(res06,23.8GHz,H) (16-bit unsigned integer)
  SUBDATASET_12_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,23.8GHz,V)
  SUBDATASET_12_DESC=[2013x243] //Brightness_Temperature_(res06,23.8GHz,V) (16-bit unsigned integer)
  SUBDATASET_13_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,36.5GHz,H)
  SUBDATASET_13_DESC=[2013x243] //Brightness_Temperature_(res06,36.5GHz,H) (16-bit unsigned integer)
  SUBDATASET_14_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,36.5GHz,V)
  SUBDATASET_14_DESC=[2013x243] //Brightness_Temperature_(res06,36.5GHz,V) (16-bit unsigned integer)
  SUBDATASET_15_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,6.9GHz,H)
  SUBDATASET_15_DESC=[2013x243] //Brightness_Temperature_(res06,6.9GHz,H) (16-bit unsigned integer)
  SUBDATASET_16_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,6.9GHz,V)
  SUBDATASET_16_DESC=[2013x243] //Brightness_Temperature_(res06,6.9GHz,V) (16-bit unsigned integer)
  SUBDATASET_17_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,7.3GHz,H)
  SUBDATASET_17_DESC=[2013x243] //Brightness_Temperature_(res06,7.3GHz,H) (16-bit unsigned integer)
  SUBDATASET_18_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,7.3GHz,V)
  SUBDATASET_18_DESC=[2013x243] //Brightness_Temperature_(res06,7.3GHz,V) (16-bit unsigned integer)
  SUBDATASET_19_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,89.0GHz,H)
  SUBDATASET_19_DESC=[2013x243] //Brightness_Temperature_(res06,89.0GHz,H) (16-bit unsigned integer)
  SUBDATASET_20_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res06,89.0GHz,V)
  SUBDATASET_20_DESC=[2013x243] //Brightness_Temperature_(res06,89.0GHz,V) (16-bit unsigned integer)
  SUBDATASET_21_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,10.7GHz,H)
  SUBDATASET_21_DESC=[2013x243] //Brightness_Temperature_(res10,10.7GHz,H) (16-bit unsigned integer)
  SUBDATASET_22_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,10.7GHz,V)
  SUBDATASET_22_DESC=[2013x243] //Brightness_Temperature_(res10,10.7GHz,V) (16-bit unsigned integer)
  SUBDATASET_23_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,18.7GHz,H)
  SUBDATASET_23_DESC=[2013x243] //Brightness_Temperature_(res10,18.7GHz,H) (16-bit unsigned integer)
  SUBDATASET_24_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,18.7GHz,V)
  SUBDATASET_24_DESC=[2013x243] //Brightness_Temperature_(res10,18.7GHz,V) (16-bit unsigned integer)
  SUBDATASET_25_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,23.8GHz,H)
  SUBDATASET_25_DESC=[2013x243] //Brightness_Temperature_(res10,23.8GHz,H) (16-bit unsigned integer)
  SUBDATASET_26_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,23.8GHz,V)
  SUBDATASET_26_DESC=[2013x243] //Brightness_Temperature_(res10,23.8GHz,V) (16-bit unsigned integer)
  SUBDATASET_27_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,36.5GHz,H)
  SUBDATASET_27_DESC=[2013x243] //Brightness_Temperature_(res10,36.5GHz,H) (16-bit unsigned integer)
  SUBDATASET_28_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,36.5GHz,V)
  SUBDATASET_28_DESC=[2013x243] //Brightness_Temperature_(res10,36.5GHz,V) (16-bit unsigned integer)
  SUBDATASET_29_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,89.0GHz,H)
  SUBDATASET_29_DESC=[2013x243] //Brightness_Temperature_(res10,89.0GHz,H) (16-bit unsigned integer)
  SUBDATASET_30_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res10,89.0GHz,V)
  SUBDATASET_30_DESC=[2013x243] //Brightness_Temperature_(res10,89.0GHz,V) (16-bit unsigned integer)
  SUBDATASET_31_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res23,18.7GHz,H)
  SUBDATASET_31_DESC=[2013x243] //Brightness_Temperature_(res23,18.7GHz,H) (16-bit unsigned integer)
  SUBDATASET_32_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res23,18.7GHz,V)
  SUBDATASET_32_DESC=[2013x243] //Brightness_Temperature_(res23,18.7GHz,V) (16-bit unsigned integer)
  SUBDATASET_33_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res23,23.8GHz,H)
  SUBDATASET_33_DESC=[2013x243] //Brightness_Temperature_(res23,23.8GHz,H) (16-bit unsigned integer)
  SUBDATASET_34_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res23,23.8GHz,V)
  SUBDATASET_34_DESC=[2013x243] //Brightness_Temperature_(res23,23.8GHz,V) (16-bit unsigned integer)
  SUBDATASET_35_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res23,36.5GHz,H)
  SUBDATASET_35_DESC=[2013x243] //Brightness_Temperature_(res23,36.5GHz,H) (16-bit unsigned integer)
  SUBDATASET_36_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res23,36.5GHz,V)
  SUBDATASET_36_DESC=[2013x243] //Brightness_Temperature_(res23,36.5GHz,V) (16-bit unsigned integer)
  SUBDATASET_37_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res23,89.0GHz,H)
  SUBDATASET_37_DESC=[2013x243] //Brightness_Temperature_(res23,89.0GHz,H) (16-bit unsigned integer)
  SUBDATASET_38_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res23,89.0GHz,V)
  SUBDATASET_38_DESC=[2013x243] //Brightness_Temperature_(res23,89.0GHz,V) (16-bit unsigned integer)
  SUBDATASET_39_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res36,36.5GHz,H)
  SUBDATASET_39_DESC=[2013x243] //Brightness_Temperature_(res36,36.5GHz,H) (16-bit unsigned integer)
  SUBDATASET_40_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res36,36.5GHz,V)
  SUBDATASET_40_DESC=[2013x243] //Brightness_Temperature_(res36,36.5GHz,V) (16-bit unsigned integer)
  SUBDATASET_41_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res36,89.0GHz,H)
  SUBDATASET_41_DESC=[2013x243] //Brightness_Temperature_(res36,89.0GHz,H) (16-bit unsigned integer)
  SUBDATASET_42_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Brightness_Temperature_(res36,89.0GHz,V)
  SUBDATASET_42_DESC=[2013x243] //Brightness_Temperature_(res36,89.0GHz,V) (16-bit unsigned integer)
  SUBDATASET_43_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Earth_Azimuth
  SUBDATASET_43_DESC=[2013x243] //Earth_Azimuth (16-bit integer)
  SUBDATASET_44_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Earth_Incidence
  SUBDATASET_44_DESC=[2013x243] //Earth_Incidence (16-bit integer)
  SUBDATASET_45_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Land_Ocean_Flag_6_to_36
  SUBDATASET_45_DESC=[4x2013x243] //Land_Ocean_Flag_6_to_36 (8-bit unsigned character)
  SUBDATASET_55_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Sun_Azimuth
  SUBDATASET_55_DESC=[2013x243] //Sun_Azimuth (16-bit integer)
  SUBDATASET_56_NAME=HDF5:"GW1AM2_201407010010_183D_L1SGRTBR_1110110.h5"://Sun_Elevation
  SUBDATASET_56_DESC=[2013x243] //Sun_Elevation (16-bit integer)
'''

class Mapper(VRT):
    ''' Mapper for SeaWIFS/MODIS/MERIS/VIIRS L2 data from OBPG

    TODO:
    * Test on SeaWIFS
    * Test on MODIS Terra
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 GCP_STEP=20, MAX_LAT=90, MIN_LAT=50, resolution='low',
                 **kwargs):
        ''' Create VRT
        Parameters
        ----------
        GCP_COUNT : int
            number of GCPs along each dimention
        '''
        ifile = os.path.split(fileName)[1]
        if not ifile.startswith('GW1AM2_') or not ifile.endswith('.h5'):
            raise WrongMapperError
        # should raise error in case of not obpg_l2 file
        try:
            ProductName = gdalMetadata['ProductName']
            PlatformShortName = gdalMetadata['PlatformShortName']
            SensorShortName = gdalMetadata['SensorShortName']
        except:
            raise WrongMapperError

        if (not ProductName == 'AMSR2-L1R' or
            not PlatformShortName == 'GCOM-W1' or
            not SensorShortName == 'AMSR2'):
            raise WrongMapperError

        if resolution == 'low':
            subDatasetWidth = 243
        else:
            subDatasetWidth = 486

        # get GCPs from lon/lat grids
        latGrid = gdal.Open('HDF5:"%s"://Latitude_of_Observation_Point_for_89A' % fileName).ReadAsArray()
        lonGrid = gdal.Open('HDF5:"%s"://Longitude_of_Observation_Point_for_89A' % fileName).ReadAsArray()
        if subDatasetWidth == 243:
            latGrid = latGrid[:, ::2]
            lonGrid = lonGrid[:, ::2]

        dx = .5
        dy = .5
        gcps = []
        k = 0
        maxY = 0
        minY = latGrid.shape[0]
        for i0 in range(0, latGrid.shape[0], GCP_STEP):
            for i1 in range(0, latGrid.shape[1], GCP_STEP):
                # create GCP with X,Y,pixel,line from lat/lon matrices
                lon = float(lonGrid[i0, i1])
                lat = float(latGrid[i0, i1])
                if (lon >= -180 and
                    lon <= 180 and
                    lat >= MIN_LAT and
                    lat <= MAX_LAT):
                    gcp = gdal.GCP(lon, lat, 0, i1 + dx, i0 + dy)
                    gcps.append(gcp)
                    k += 1
                    maxY = max(maxY, i0)
                    minY = min(minY, i0)
        yOff = minY
        ySize = maxY - minY

        # remove Y-offset from gcps
        for gcp in gcps:
            gcp.GCPLine -= yOff

        metaDict = []

        subDatasets = gdalDataset.GetSubDatasets()
        metadata = gdalDataset.GetMetadata()
        for subDataset in subDatasets:
            # select subdatasets fro that resolution (width)
            if (subDatasetWidth == int(subDataset[1].split(']')[0].split('x')[-1]) and
                'Latitude' not in subDataset[0] and 'Longitude' not in subDataset[0]):
                name = subDataset[0].split('/')[-1]
                # find scale
                scale = 1
                for meta in metadata:
                    if name + '_SCALE' in meta:
                        scale = float(metadata[meta])
                # create meta entry
                metaEntry = {'src': {'SourceFilename': subDataset[0],
                                     'sourceBand':  1,
                                     'ScaleRatio': scale,
                                     'ScaleOffset': 0,
                                     'yOff': yOff,
                                     'ySize': ySize,},
                             'dst': {'name': name}
                             }
                metaDict.append(metaEntry)

        # create VRT from one of the subdatasets
        gdalSubDataset = gdal.Open(metaEntry['src']['SourceFilename'])
        VRT.__init__(self, srcRasterXSize=subDatasetWidth, srcRasterYSize=ySize)
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)


        # append GCPs and lat/lon projection to the vsiDataset
        self.dataset.SetGCPs(gcps, NSR().wkt)
        self.reproject_GCPs('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
        self.tps = True
