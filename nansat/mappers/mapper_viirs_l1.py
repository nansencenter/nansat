# Name:        mapper_modisL1
# Purpose:     Mapping for MODIS-L1 data
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from __future__ import absolute_import, unicode_literals, division

import os
import glob
from datetime import datetime, timedelta
from math import ceil
try:
    from scipy.ndimage.filters import gaussian_filter
except:
    IMPORT_SCIPY = False
else:
    IMPORT_SCIPY = True

from nansat.nsr import NSR
from nansat.vrt import VRT
from nansat.tools import gdal, ogr
from nansat.exceptions import WrongMapperError, NansatReadError


class Mapper(VRT):
    ''' VRT with mapping of WKV for VIIRS Level 1B '''

    def __init__(self, filename, gdalDataset, gdalMetadata,
                 GCP_COUNT0=5, GCP_COUNT1=20, pixelStep=1,
                 lineStep=1, **kwargs):
        ''' Create VIIRS VRT '''

        if not 'GMTCO_npp_' in filename:
            raise WrongMapperError(filename)
        ifiledir = os.path.split(filename)[0]
        ifiles = glob.glob(ifiledir + 'SVM??_npp_d*_obpg_ops.h5')
        ifiles.sort()

        if not IMPORT_SCIPY:
            raise NansatReadError('VIIRS data cannot be read because scipy is not installed')

        viirsWavelengths = [None, 412, 445, 488, 555, 672, 746, 865, 1240,
                            1378, 1610, 2250, 3700, 4050, 8550, 10736, 12013]

        # create empty VRT dataset with geolocation only
        xDatasetSource = ('HDF5:"%s"://All_Data/VIIRS-MOD-GEO-TC_All/Longitude'
                          % filename)
        xDatasetBand = 1
        xDataset = gdal.Open(xDatasetSource)
        self._init_from_gdal_dataset(xDataset)

        metaDict = []
        for ifile in ifiles:
            ifilename = os.path.split(ifile)[1]
            print(ifilename)
            bNumber = int(ifilename[3:5])
            print(bNumber)
            bWavelength = viirsWavelengths[bNumber]
            print(bWavelength)
            SourceFilename = ('HDF5:"%s"://All_Data/VIIRS-M%d-SDR_All/Radiance'
                              % (ifile, bNumber))
            print(SourceFilename)
            metaEntry = {'src': {'SourceFilename': SourceFilename,
                         'SourceBand': 1},
                         'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                                 'wavelength': str(bWavelength),
                                 'suffix': str(bWavelength)}
                         }
            metaDict.append(metaEntry)

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        xVRTArray = xDataset.ReadAsArray()
        xVRTArray = gaussian_filter(xVRTArray, 5).astype('float32')
        xVRT = VRT.from_array(xVRTArray)

        yDatasetSource = ('HDF5:"%s"://All_Data/VIIRS-MOD-GEO-TC_All/Latitude'
                          % filename)
        yDatasetBand = 1
        yDataset = gdal.Open(yDatasetSource)
        yVRTArray = yDataset.ReadAsArray()
        yVRTArray = gaussian_filter(yVRTArray, 5).astype('float32')
        yVRT = VRT.from_array(yVRTArray)

        # estimate pixel/line step
        self.logger.debug('pixel/lineStep %f %f' % (pixelStep, lineStep))

        # ==== ADD GCPs and Pojection ====
        # get lat/lon matrices
        longitude = xVRT.dataset.GetRasterBand(1).ReadAsArray()
        latitude = yVRT.dataset.GetRasterBand(1).ReadAsArray()

        # estimate step of GCPs
        step0 = max(1, int(float(latitude.shape[0]) / GCP_COUNT0))
        step1 = max(1, int(float(latitude.shape[1]) / GCP_COUNT1))
        self.logger.debug('gcpCount: %d %d %d %d, %d %d',
                          latitude.shape[0], latitude.shape[1],
                          GCP_COUNT0, GCP_COUNT1, step0, step1)

        # generate list of GCPs
        gcps = []
        k = 0
        for i0 in range(0, latitude.shape[0], step0):
            for i1 in range(0, latitude.shape[1], step1):
                # create GCP with X,Y,pixel,line from lat/lon matrices
                lon = float(longitude[i0, i1])
                lat = float(latitude[i0, i1])
                if (lon >= -180 and lon <= 180 and lat >= -90 and lat <= 90):
                    gcp = gdal.GCP(lon, lat, 0,
                                   i1 * pixelStep, i0 * lineStep)
                    self.logger.debug('%d %d %d %f %f',
                                      k, gcp.GCPPixel, gcp.GCPLine,
                                      gcp.GCPX, gcp.GCPY)
                    gcps.append(gcp)
                    k += 1

        # append GCPs and lat/lon projection to the vsiDataset
        self.dataset.SetGCPs(gcps, NSR().wkt)

        # remove geolocation array
        self._remove_geolocation()
        #"""
