# Name:        mapper_globcolour_l3b
# Purpose:     Mapping for GLOBCOLOUR L3B data
# Authors:     Anton Korosov
# Licence:     This file is part of NANSAT. You can redistribute it or modify
#              under the terms of GNU General Public License, v.3
#              http://www.gnu.org/licenses/gpl-3.0.html
from __future__ import print_function
import glob
import os
import datetime
import json
import tempfile
import shutil

import numpy as np

from netCDF4 import Dataset

import pythesint as pti

from nansat.exceptions import WrongMapperError
from nansat.vrt import VRT
from nansat.mappers.globcolour import Globcolour


class Mapper(VRT, Globcolour):
    ''' Create VRT with mapping of WKV for MERIS Level 2 (FR or RR)'''

    def __init__(self, filename, gdalDataset, gdalMetadata, latlonGrid=None,
                 mask='', **kwargs):

        ''' Create MER2 VRT

        Parameters
        -----------
        filename : string
        gdalDataset : gdal dataset
        gdalMetadata : gdal metadata
        latlonGrid : numpy 2 layered 2D array with lat/lons of desired grid
        '''
        # test if input files is GLOBCOLOUR L3B
        iDir, iFile = os.path.split(filename)
        iFileName, iFileExt = os.path.splitext(iFile)
        #print 'idir:', iDir, iFile, iFileName[0:5], iFileExt[0:8]
        if (iFileName[0:4] != 'L3b_'
            or iFileExt != '.nc'
            or not os.path.exists(filename)
            or (gdalDataset is not None
                and (len(gdalDataset.GetSubDatasets()) > 0
                     or gdalDataset.RasterCount > 0))):
            raise WrongMapperError

        # define shape of GLOBCOLOUR grid
        GLOBCOLOR_ROWS = 180 * 24
        GLOBCOLOR_COLS = 360 * 24

        # define lon/lat grids for projected var
        if latlonGrid is None:
            latlonGrid = np.mgrid[90:-90:4320j,
                                  -180:180:8640j].astype('float32')
            #latlonGrid = np.mgrid[80:50:900j, -10:30:1200j].astype('float16')
            #latlonGrid = np.mgrid[47:39:300j, 25:45:500j].astype('float32')

        # create empty VRT dataset with geolocation only
        self._init_from_lonlat(latlonGrid[1], latlonGrid[0])

        # get list of similar (same date) files in the directory
        simFilesMask = os.path.join(iDir, iFileName[0:30] + '*' + mask + '.nc')
        simFiles = glob.glob(simFilesMask)
        simFiles.sort()

        metaDict = []
        self.band_vrts = {'mask': [], 'lonlat': []}
        mask = None
        for simFile in simFiles:
            print('sim: ', simFile)
            # copy simFile to a temporary file
            tmpf = tempfile.mkstemp()
            shutil.copyfile(simFile, tmpf[1])
            f = Dataset(tmpf[1])

            # get iBinned, index for converting from binned into GLOBCOLOR-grid
            colBinned = f.variables['col'][:]
            rowBinned = f.variables['row'][:]
            iBinned = (colBinned.astype('uint32') +
                      (rowBinned.astype('uint32') - 1) * GLOBCOLOR_COLS)
            colBinned = None
            rowBinned = None

            # get iRawPro, index for converting
            # from GLOBCOLOR-grid to latlonGrid
            yRawPro = np.rint(1 + (GLOBCOLOR_ROWS - 1) *
                              (latlonGrid[0] + 90) / 180.)
            lon_step_Mat = 24. * np.cos(np.pi * latlonGrid[0] / 180.)
            xRawPro = np.rint(1 + (latlonGrid[1] + 180) * lon_step_Mat)
            iRawPro = xRawPro.astype('uint32') + (yRawPro.astype('uint32') - 1) * GLOBCOLOR_COLS
            yRawPro = None
            xRawPro = None

            for varName in f.variables:
                # find variable with _mean, eg CHL1_mean
                if '_mean' in varName:
                    var = f.variables[varName]
                    break

            # skip variable if no WKV is give in Globcolour
            if varName not in self.varname2wkv:
                continue

            # get WKV
            varWKV = self.varname2wkv[varName]

            # read binned data
            varBinned = var[:]

            # convert to GLOBCOLOR grid
            varRawPro = np.zeros([GLOBCOLOR_ROWS, GLOBCOLOR_COLS], 'float32')
            varRawPro.flat[iBinned] = varBinned

            # convert to latlonGrid
            varPro = varRawPro.flat[iRawPro.flat[:]].reshape(iRawPro.shape)

            # add mask band
            if mask is None:
                mask = np.zeros(varPro.shape, 'uint8')
                mask[:] = 1
                mask[varPro > 0] = 64

                # add VRT with array with data from projected variable
                self.band_vrts['mask'].append(VRT.from_array(mask))

                # add metadata to the dictionary
                metaDict.append({
                    'src': {'SourceFilename': (self.band_vrts['mask'][-1].
                                               filename),
                            'SourceBand':  1},
                    'dst': {'name': 'mask'}})

            # add VRT with array with data from projected variable
            self.band_vrts['lonlat'].append(VRT.from_array(varPro))

            # add metadata to the dictionary
            metaEntry = {
                'src': {'SourceFilename': self.band_vrts['lonlat'][-1].filename,
                        'SourceBand':  1},
                'dst': {'wkv': varWKV, 'original_name': varName}}

            # add wavelength for nLw
            longName = 'Fully normalised water leaving radiance'
            if longName in f.variables[varName].long_name:
                simWavelength = varName.split('L')[1].split('_mean')[0]
                metaEntry['dst']['suffix'] = simWavelength
                metaEntry['dst']['wavelength'] = simWavelength

            # add all metadata from NC-file
            for attr in var.ncattrs():
                metaEntry['dst'][attr] = var.getncattr(attr)

            metaDict.append(metaEntry)

            # add Rrsw band
            metaEntry2 = self.make_rrsw_meta_entry(metaEntry)
            if metaEntry2 is not None:
                metaDict.append(metaEntry2)

            os.remove(tmpf[1])

        instrument = f.title.strip().split(' ')[-2].split('/')[0]
        mm = pti.get_gcmd_instrument(instrument)
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))

        platform = {
            'MODIS' : 'AQUA',
            'MERIS' : 'ENVISAT',
            'SEAWIFS': 'QUICKBIRD',
            'VIIRS' : 'SUOMI-NPP'}[instrument.upper()]
        pp = pti.get_gcmd_platform(platform)
        self.dataset.SetMetadataItem('platform', json.dumps(pp))

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # add time
        startDate = datetime.datetime(int(iFileName[4:8]),
                                      int(iFileName[8:10]),
                                      int(iFileName[10:12]))

        # Adding valid time to dataset
        self.dataset.SetMetadataItem('time_coverage_start', startDate.isoformat())
        self.dataset.SetMetadataItem('time_coverage_end', startDate.isoformat())
