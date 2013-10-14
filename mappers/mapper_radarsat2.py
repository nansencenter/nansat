# Name:         mapper_radarsat2
# Purpose:      Nansat mapping for Radarsat2 data
# Authors:      Morten W. Hansen, Knut-Frode Dagestad, Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from numpy import mod
from math import asin
import tarfile
import zipfile

from vrt import *
from domain import Domain
from nansat_tools import Node

try:
    from osgeo import gdal
except ImportError:
    import gdal

import pdb


class Mapper(VRT):
    ''' Create VRT with mapping of WKV for Radarsat2 '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create Radarsat2 VRT '''
        fPathName, fExt = os.path.splitext(fileName)

        if zipfile.is_zipfile(fileName):
            # Open zip file using VSI
            fPath, fName = os.path.split(fPathName)
            fileName = '/vsizip/%s/%s' % (fileName, fName)
            gdalDataset = gdal.Open(fileName)
            gdalMetadata = gdalDataset.GetMetadata()

        #if it is not RADARSAT-2, return
        product = gdalMetadata.get("SATELLITE_IDENTIFIER", "Not_RADARSAT-2")
        if product != 'RADARSAT-2':
            raise AttributeError("RADARSAT-2 BAD MAPPER")
        
        # read product.xml
        productXmlName = os.path.join(fileName, 'product.xml')
        productXml = self.read_xml(productXmlName)

        # Get additional metadata from product.xml
        rs2_0 = Node.create(productXml)
        rs2_1 = rs2_0.node('sourceAttributes')
        rs2_2 = rs2_1.node('radarParameters')
        antennaPointing = 90 if rs2_2['antennaPointing'].lower() =='right' else -90
        rs2_3 = rs2_1.node('orbitAndAttitude').node('orbitInformation')     
        passDirection = rs2_3['passDirection']

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        #define dictionary of metadata and band specific parameters
        pol = []
        metaDict = []

        # Get the subdataset with calibrated sigma0 only
        for dataset in gdalDataset.GetSubDatasets():
            if dataset[1] == 'Sigma Nought calibrated':
                s0dataset = gdal.Open(dataset[0])
                s0datasetName = dataset[0][:]
                s0datasetPol = s0dataset.GetRasterBand(1).\
                    GetMetadata()['POLARIMETRIC_INTERP']
                for i in range(1, s0dataset.RasterCount+1):
                    polString = s0dataset.GetRasterBand(i).\
                        GetMetadata()['POLARIMETRIC_INTERP']
                    suffix = polString
                    # The nansat data will be complex if the SAR data is of type 10
                    dtype = s0dataset.GetRasterBand(i).DataType
                    if dtype==10:
                        # add intensity band
                        metaDict.append({
                            'src': {
                                'SourceFilename': 'RADARSAT_2_CALIB:SIGMA0:' + fileName + '/product.xml',
                                'SourceBand': i,
                                'DataType': dtype,
                            },
                            'dst': {
                                'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                                'PixelFunctionType': 'intensity',
                                'SourceTransferType': \
                                    gdal.GetDataTypeName(dtype),
                                'suffix': suffix,
                                'polarization': polString,
                                'dataType': 6,
                            }
                        })
                        # modify suffix for adding the compled band below
                        suffix = polString+'_complex'
                    pol.append(polString)
                    metaDict.append({
                        'src': {
                            'SourceFilename': ('RADARSAT_2_CALIB:SIGMA0:'
                                                    + fileName +
                                                    '/product.xml'),
                            'SourceBand': i,
                            'DataType': dtype,
                        },
                        'dst': {
                            'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                            'suffix': suffix,
                            'polarization': polString
                        }
                    })
            if dataset[1] == 'Beta Nought calibrated':
                b0dataset = gdal.Open(dataset[0])
                b0datasetName = dataset[0][:]
                for j in range(1, b0dataset.RasterCount+1):
                    polString = b0dataset.GetRasterBand(j).\
                        GetMetadata()['POLARIMETRIC_INTERP']
                    if polString == s0datasetPol:
                        b0datasetBand = j

        # Add Sigma0 bands with metadata
        self._create_bands(metaDict)

        # Add derived band (incidence angle) calculated using pixel function
        # "BetaSigmaToIncidence":
        src = [
            {
                'SourceFilename': b0datasetName,
                'SourceBand':  b0datasetBand,
                'DataType': dtype
            },
            {
                'SourceFilename': s0datasetName,
                'SourceBand': 1,
                'DataType': dtype
            }
        ]
        dst = {
            'wkv': 'angle_of_incidence',
            'PixelFunctionType': 'BetaSigmaToIncidence',
            'SourceTransferType': gdal.GetDataTypeName(dtype),
            '_FillValue': -10000,   # NB: this is also hard-coded in
                                    #     pixelfunctions.c
            'dataType': 6,
            'name': 'incidence_angle'
        }

        self._create_band(src, dst)
        self.dataset.FlushCache()

        ###################################################################
        # Add sigma0_VV - pixel function of sigma0_HH and beta0_HH
        # incidence angle is calculated within pixel function
        # It is assummed that HH is the first band in sigma0 and
        # beta0 sub datasets
        ###################################################################
        if 'VV' not in pol and 'HH' in pol:
            s0datasetNameHH = pol.index('HH')+1
            src = [{'SourceFilename': s0datasetName,
                    'SourceBand': s0datasetNameHH,
                    'DataType': 6},
                   {'SourceFilename': b0datasetName,
                    'SourceBand': b0datasetBand,
                    'DataType': 6}]
            dst = {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                   'PixelFunctionType': 'Sigma0HHBetaToSigma0VV',
                   'polarization': 'VV',
                   'suffix': 'VV'}
            self._create_band(src, dst)
            self.dataset.FlushCache()

        ############################################
        # Add SAR look direction to metadata domain
        ############################################
        self.dataset.SetMetadataItem('SAR_center_look_direction', str(mod(
            Domain(ds=gdalDataset).upwards_azimuth_direction( orbit_direction =
            str(passDirection) ) + antennaPointing, 360)))

        # Set time
        validTime = gdalDataset.GetMetadata()['ACQUISITION_START_TIME']
        self.logger.info('Valid time: %s', str(validTime))
        self._set_time(parse(validTime))
