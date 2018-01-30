# Name:        mapper_aster_l1a
# Purpose:     Mapping for L2 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from datetime import datetime, timedelta
from math import ceil
from dateutil.parser import parse

from nansat.tools import gdal, ogr
from nansat.exceptions import WrongMapperError
from nansat.vrt import VRT
from nansat.nsr import NSR


class Mapper(VRT):
    ''' Mapper for ASTER L1A VNIR data'''

    def __init__(self, filename, gdalDataset, gdalMetadata,
                 GCP_COUNT=10,
                 bandNames=['VNIR_Band1', 'VNIR_Band2', 'VNIR_Band3N'],
                 bandWaves=[560, 660, 820], **kwargs):
        ''' Create VRT

        Parameters
        -----------
        GCP_COUNT : int
            number of GCPs along each dimention
        bandNames : list of string (band name)
        bandWaves : list of integer (waves corresponding to band name)

        Band name and waves
        --------------------
        'VNIR_Band3B' : 820, 'SWIR_Band4' : 1650, 'SWIR_Band5' : 2165,
        'SWIR_Band6' : 2205, 'SWIR_Band7' : 2260, 'SWIR_Band8' : 2330,
        'SWIR_Band9' : 2395, 'TIR_Band10' : 8300, 'TIR_Band11' : 8650,
        'TIR_Band12' : 9100, 'TIR_Band13' : 10600, 'TIR_Band14' : 11300

        '''
        # check if it is ASTER L1A
        try:
            assert 'AST_L1A_' in filename
            shortName = gdalMetadata['INSTRUMENTSHORTNAME']
            assert shortName == 'ASTER'
        except:
            raise WrongMapperError

        subDatasets = gdalDataset.GetSubDatasets()

        # find datasets for each band and generate metaDict
        metaDict = []
        bandDatasetMask = 'HDF4_EOS:EOS_SWATH:"%s":%s:ImageData'
        for bandName, bandWave in zip(bandNames, bandWaves):
            metaEntry = {'src': {'SourceFilename': (bandDatasetMask
                                                    % (filename, bandName)),
                                 'SourceBand': 1,
                                 'DataType': 6,
                                 },
                         'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                                 'wavelength': str(bandWave),
                                 'suffix': str(bandWave),
                                 }}
            metaDict.append(metaEntry)

        # create empty VRT dataset with geolocation only
        gdalSubDataset = gdal.Open(metaDict[0]['src']['SourceFilename'])
        self._init_from_gdal_dataset(gdalSubDataset, metadata=gdalSubDataset.GetMetadata())

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # find largest lon/lat subdatasets
        latShape0 = 0
        for subDataset in subDatasets:
            if 'Latitude' in subDataset[1]:
                ls = int(subDataset[1].strip().split('[')[1].split('x')[0])
                if ls >= latShape0:
                    latShape0 = ls
                    latSubDS = subDataset[0]
            if 'Longitude' in subDataset[1]:
                ls = int(subDataset[1].strip().split('[')[1].split('x')[0])
                if ls >= latShape0:
                    latShape0 = ls
                    lonSubDS = subDataset[0]
        self.logger.debug(latSubDS)
        self.logger.debug(lonSubDS)

        # get lat/lon matrices
        xDataset = gdal.Open(lonSubDS)
        yDataset = gdal.Open(latSubDS)

        longitude = xDataset.ReadAsArray()
        latitude = yDataset.ReadAsArray()

        step0 = longitude.shape[0] / GCP_COUNT
        step1 = longitude.shape[1] / GCP_COUNT

        # estimate pixel/line step
        pixelStep = int(ceil(float(gdalSubDataset.RasterXSize) /
                             float(xDataset.RasterXSize)))
        lineStep = int(ceil(float(gdalSubDataset.RasterYSize) /
                            float(xDataset.RasterYSize)))
        self.logger.debug('steps: %d %d %d %d' % (step0, step1,
                                                  pixelStep, lineStep))

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
                    self.logger.debug('%d %d %d %f %f' % (k, gcp.GCPPixel,
                                                          gcp.GCPLine,
                                                          gcp.GCPX, gcp.GCPY))
                    gcps.append(gcp)
                    k += 1
        # append GCPs and lat/lon projection to the vsiDataset
        self.dataset.SetGCPs(gcps, NSR().wkt)

        # Adding valid time to dataset
        self.dataset.SetMetadataItem('time_coverage_start',
                                     parse(gdalMetadata['FIRSTPACKETTIME']).isoformat())
        self.dataset.SetMetadataItem('time_coverage_end',
                                     parse(gdalMetadata['LASTPACKETTIME']).isoformat())

        mm = pti.get_gcmd_instrument('ASTER')
        ee = pti.get_gcmd_platform('TERRA')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
