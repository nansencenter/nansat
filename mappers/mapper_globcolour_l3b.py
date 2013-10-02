# Name:        mapper_globcolour_l3b
# Purpose:     Mapping for GLOBCOLOUR L3B data
# Authors:     Anton Korosov
# Licence:     This file is part of NANSAT. You can redistribute it or modify
#              under the terms of GNU General Public License, v.3
#              http://www.gnu.org/licenses/gpl-3.0.html
import glob
import os.path

from scipy.io.netcdf import netcdf_file
import numpy as np

from vrt import VRT


class Mapper(VRT):
    ''' Create VRT with mapping of WKV for MERIS Level 2 (FR or RR)'''

    def __init__(self, fileName, gdalDataset, gdalMetadata, latlonGrid=None, **kwargs):

        ''' Create MER2 VRT

        Parameters
        -----------
        fileName : string
        gdalDataset : gdal dataset
        gdalMetadata : gdal metadata
        latlonGrid : numpy 2 layered 2D array with lat/lons of desired grid
        '''
        # define lon/lat grids for projected var
        if latlonGrid is None:
            #latlonGrid = np.mgrid[90:-90:4320j, -180:180:8640j].astype('float16')
            #latlonGrid = np.mgrid[80:50:900j, -10:30:1200j].astype('float16')
            latlonGrid = np.mgrid[47:39:300j, 25:45:500j].astype('float32')

        # test if files is GLOBCOLOUR L3B
        iDir, iFile = os.path.split(fileName)
        iFileName, iFileExt = os.path.splitext(iFile)
        print 'idir:', iDir, iFile, iFileName[0:5], iFileExt[0:8]
        assert iFileName[0:4] == 'L3b_' and iFileExt == '.nc'

        # get list of similar (same date) files in the directory
        simFilesMask = os.path.join(iDir, iFileName[0:30] + '*')
        simFiles = glob.glob(simFilesMask)
        print 'simFilesMask, simFiles', simFilesMask, simFiles

        metaDict = []
        for simFile in simFiles:
            print 'simFile', simFile
            f = netcdf_file(simFile)

            # get iBinned, index for converting from binned into raw-gridded
            colBinned = f.variables['col'][:]
            rowBinned = f.variables['row'][:]
            GLOBCOLOR_ROWS = 180 * 24
            GLOBCOLOR_COLS = 360 * 24
            iBinned = colBinned.astype('uint32') + (rowBinned.astype('uint32') - 1) * GLOBCOLOR_COLS
            colBinned = None
            rowBinned = None
            
            # get iRawPro, index for converting from raw-gridded to proj-gridded
            yRawPro = np.rint(1 + (GLOBCOLOR_ROWS - 1) * (latlonGrid[0] + 90) / 180)
            lon_step_Mat = 1 / np.cos(np.pi * latlonGrid[0] / 180.) / 24.
            xRawPro = np.rint(1 + (latlonGrid[1] + 180) / lon_step_Mat)
            iRawPro = xRawPro + (yRawPro - 1) * GLOBCOLOR_COLS
            iRawPro[iRawPro < 0] = 0
            iRawPro = np.rint(iRawPro).astype('uint32')
            yRawPro = None
            xRawPro = None
            
            print iRawPro.shape, iRawPro
