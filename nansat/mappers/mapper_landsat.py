# Name:         mapper_landsat
# Purpose:      Mapping for LANDSAT*.tar.gz
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
import glob
import tarfile
import warnings
import datetime
import json

import pythesint as pti

from nansat.tools import WrongMapperError, parse_time
from nansat.tools import gdal, np
from nansat.vrt import VRT

class Mapper(VRT):
    ''' Mapper for LANDSAT5,6,7,8 .tar.gz or tif files'''

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                       resolution='low', **kwargs):
        ''' Create LANDSAT VRT from multiple tif files or single tar.gz file'''
        mtlFileName = ''
        bandFileNames = []
        bandSizes = []
        bandDatasets = []
        fname = os.path.split(fileName)[1]

        if   (fileName.endswith('.tar') or
              fileName.endswith('.tar.gz') or
              fileName.endswith('.tgz')):
            # try to open .tar or .tar.gz or .tgz file with tar
            try:
                tarFile = tarfile.open(fileName)
            except:
                raise WrongMapperError

            # collect names of bands and corresponding sizes
            # into bandsInfo dict and bandSizes list
            tarNames = sorted(tarFile.getnames())
            for tarName in tarNames:
                # check if TIF files inside TAR qualify
                if   (tarName[0] in ['L', 'M'] and
                      os.path.splitext(tarName)[1] in ['.TIF', '.tif']):
                    # open TIF file from TAR using VSI
                    sourceFilename = '/vsitar/%s/%s' % (fileName, tarName)
                    gdalDatasetTmp = gdal.Open(sourceFilename)
                    # keep name, GDALDataset and size
                    bandFileNames.append(sourceFilename)
                    bandSizes.append(gdalDatasetTmp.RasterXSize)
                    bandDatasets.append(gdalDatasetTmp)
                elif (tarName.endswith('MTL.txt') or
                      tarName.endswith('MTL.TXT')):
                    # get mtl file
                    mtlFileName = tarName

        elif ((fname.startswith('L') or fname.startswith('M')) and
              (fname.endswith('.tif') or
               fname.endswith('.TIF') or
               fname.endswith('._MTL.txt'))):

            # try to find TIF/tif files with the same name as input file
            path, coreName = os.path.split(fileName)
            coreName = os.path.splitext(coreName)[0].split('_')[0]
            coreNameMask = coreName+'*[tT][iI][fF]'
            tifNames = sorted(glob.glob(os.path.join(path, coreNameMask)))
            for tifName in tifNames:
                sourceFilename = tifName
                gdalDatasetTmp = gdal.Open(sourceFilename)
                # keep name, GDALDataset and size
                bandFileNames.append(sourceFilename)
                bandSizes.append(gdalDatasetTmp.RasterXSize)
                bandDatasets.append(gdalDatasetTmp)

            # get mtl file
            mtlFiles = glob.glob(coreName+'*[mM][tT][lL].[tT][xX][tT]')
            if len(mtlFiles) > 0:
                mtlFileName = mtlFiles[0]
        else:
            raise WrongMapperError

        # if not TIF files found - not appropriate mapper
        if not bandFileNames:
            raise WrongMapperError

        # get appropriate band size based on number of unique size and
        # required resoltuion
        if resolution == 'low':
            bandXSise = min(bandSizes)
        elif resolution in ['high', 'hi']:
            bandXSise = max(bandSizes)
        else:
            raise OptionError('Wrong resolution %s for file %s' % (resolution, fileName))

        # find bands with appropriate size and put to metaDict
        metaDict = []
        for bandFileName, bandSize, bandDataset in zip(bandFileNames,
                                                       bandSizes,
                                                       bandDatasets):
            if bandSize == bandXSise:
                # let last part of file name be suffix
                bandSuffix = os.path.splitext(bandFileName)[0].split('_')[-1]

                metaDict.append({
                    'src': {'SourceFilename': bandFileName,
                            'SourceBand':  1,
                            'ScaleRatio': 0.1},
                    'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                            'suffix': bandSuffix}})
                gdalDataset4Use = bandDataset

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset4Use)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        if len(mtlFileName) > 0:
            mtlFileName = os.path.join(os.path.split(bandFileNames[0])[0],
                                        mtlFileName)
            mtlFileLines = [line.strip() for line in
                            self.read_xml(mtlFileName).split('\n')]
            dateString = [line.split('=')[1].strip()
                          for line in mtlFileLines
                            if ('DATE_ACQUIRED' in line or
                              'ACQUISITION_DATE' in line)][0]
            timeStr = [line.split('=')[1].strip()
                        for line in mtlFileLines
                            if ('SCENE_CENTER_TIME' in line or
                                'SCENE_CENTER_SCAN_TIME' in line)][0]
            time_start = parse_time(dateString + 'T' + timeStr).isoformat()
            time_end = (parse_time(dateString + 'T' + timeStr) +
                        datetime.timedelta(microseconds=60000000)).isoformat()

        self.dataset.SetMetadataItem('time_coverage_start', time_start)
        self.dataset.SetMetadataItem('time_coverage_end', time_end)

        # set platform
        platform = 'LANDSAT'
        if fname[2].isdigit():
            platform += '-'+fname[2]
        ee = pti.get_gcmd_platform(platform)
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

        # set instrument
        instrument = {
        'LANDSAT' : 'MSS',
        'LANDSAT-1' : 'MSS',
        'LANDSAT-2' : 'MSS',
        'LANDSAT-3' : 'MSS',
        'LANDSAT-4' : 'TM',
        'LANDSAT-5' : 'TM',
        'LANDSAT-7' : 'ETM+',
        'LANDSAT-8' : 'OLI'}[platform]
        ee = pti.get_gcmd_instrument(instrument)
        self.dataset.SetMetadataItem('instrument', json.dumps(ee))

