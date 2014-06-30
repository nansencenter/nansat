# Name:         mapper_radarsat2
# Purpose:      Nansat mapping for Radarsat2 data
# Authors:      Morten W. Hansen, Knut-Frode Dagestad, Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
import tarfile
import zipfile
from dateutil.parser import parse

import numpy as np
import scipy.ndimage
from math import asin

from nansat.vrt import VRT
from nansat.domain import Domain
from nansat.node import Node
from nansat.tools import Error, initial_bearing

try:
    from osgeo import gdal
except ImportError:
    import gdal


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
        if rs2_2['antennaPointing'].lower() == 'right':
            antennaPointing = 90
        else:
            antennaPointing = -90
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
                band = s0dataset.GetRasterBand(1)
                s0datasetPol = band.GetMetadata()['POLARIMETRIC_INTERP']
                for i in range(1, s0dataset.RasterCount+1):
                    iBand = s0dataset.GetRasterBand(i)
                    polString = iBand.GetMetadata()['POLARIMETRIC_INTERP']
                    suffix = polString
                    # The nansat data will be complex if the SAR data is of type 10
                    dtype = iBand.DataType
                    if dtype == 10:
                        # add intensity band
                        metaDict.append(
                            {'src': {'SourceFilename': ('RADARSAT_2_CALIB:SIGMA0:'
                                                        + fileName
                                                        + '/product.xml'),
                                     'SourceBand': i,
                                     'DataType': dtype},
                             'dst': {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                                     'PixelFunctionType': 'intensity',
                                     'SourceTransferType': gdal.GetDataTypeName(dtype),
                                     'suffix': suffix,
                                     'polarization': polString,
                                     'dataType': 6}})
                        # modify suffix for adding the compled band below
                        suffix = polString+'_complex'
                    pol.append(polString)
                    metaDict.append(
                        {'src': {'SourceFilename': ('RADARSAT_2_CALIB:SIGMA0:'
                                                    + fileName
                                                    + '/product.xml'),
                                 'SourceBand': i,
                                 'DataType': dtype},
                         'dst': {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                                 'suffix': suffix,
                                 'polarization': polString}})

            if dataset[1] == 'Beta Nought calibrated':
                b0dataset = gdal.Open(dataset[0])
                b0datasetName = dataset[0][:]
                for j in range(1, b0dataset.RasterCount+1):
                    jBand = b0dataset.GetRasterBand(j)
                    polString = jBand.GetMetadata()['POLARIMETRIC_INTERP']
                    if polString == s0datasetPol:
                        b0datasetBand = j

        ###############################
        # Add SAR look direction
        ###############################
        d = Domain(ds=gdalDataset)
        lon, lat = d.get_geolocation_grids(100)

        # Calculate SAR look direction (assuming right-looking)
        SAR_look_direction = initial_bearing(lon[:, :-1], lat[:, :-1],
                  lon[:, 1:], lat[:, 1:]) + antennaPointing + 90.0
        # Interpolate to regain lost row
        SAR_look_direction = np.mod(SAR_look_direction, 360)
        SAR_look_direction = scipy.ndimage.interpolation.zoom(
                                SAR_look_direction, (1, 11./10.))
        # Decompose, to avoid interpolation errors around 0 <-> 360
        SAR_look_direction_u = np.sin(np.deg2rad(SAR_look_direction))
        SAR_look_direction_v = np.cos(np.deg2rad(SAR_look_direction))
        look_u_VRT = VRT(array=SAR_look_direction_u, lat=lat, lon=lon)
        look_v_VRT = VRT(array=SAR_look_direction_v, lat=lat, lon=lon)

        # Note: If incidence angle and look direction are stored in
        #       same VRT, access time is about twice as large
        lookVRT = VRT(lat=lat, lon=lon)
        lookVRT._create_band(
                [{'SourceFilename': look_u_VRT.fileName, 'SourceBand': 1},
                 {'SourceFilename': look_v_VRT.fileName, 'SourceBand': 1}],
                 {'PixelFunctionType': 'UVToDirectionTo'})

        # Blow up to full size
        lookVRT = lookVRT.get_resized_vrt(
                    gdalDataset.RasterXSize, gdalDataset.RasterYSize)
        # Store VRTs so that they are accessible later
        self.subVRTs = {'look_u_VRT': look_u_VRT,
                        'look_v_VRT': look_v_VRT,
                        'lookVRT': lookVRT}

        # Add band to full sized VRT
        lookFileName = self.subVRTs['lookVRT'].fileName
        metaDict.append({'src':
                            {'SourceFilename': lookFileName,
                             'SourceBand': 1},
                         'dst':
                            {'wkv': 'sensor_azimuth_angle',
                             'name': 'SAR_look_direction'}})

        ###############################
        # Create bands
        ###############################
        self._create_bands(metaDict)

        ###################################################
        # Add derived band (incidence angle) calculated
        # using pixel function "BetaSigmaToIncidence":
        ###################################################
        src = [{'SourceFilename': b0datasetName,
                'SourceBand':  b0datasetBand,
                'DataType': dtype},
               {'SourceFilename': s0datasetName,
                'SourceBand': 1,
                'DataType': dtype}]
        dst = {'wkv': 'angle_of_incidence',
               'PixelFunctionType': 'BetaSigmaToIncidence',
               'SourceTransferType': gdal.GetDataTypeName(dtype),
               '_FillValue': -10000,   # NB: this is also hard-coded in
                                       #     pixelfunctions.c
               'dataType': 6,
               'name': 'incidence_angle'}

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
        # Add SAR metadata
        ############################################
        if antennaPointing == 90:
            self.dataset.SetMetadataItem('ANTENNA_POINTING', 'RIGHT')
        if antennaPointing == -90:
            self.dataset.SetMetadataItem('ANTENNA_POINTING', 'LEFT')
        self.dataset.SetMetadataItem('ORBIT_DIRECTION', str(passDirection).upper())

        # Set time
        validTime = gdalDataset.GetMetadata()['ACQUISITION_START_TIME']
        self.logger.info('Valid time: %s', str(validTime))
        self._set_time(parse(validTime))
