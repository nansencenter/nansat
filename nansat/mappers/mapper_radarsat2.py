# Name:         mapper_radarsat2
# Purpose:      Nansat mapping for Radarsat2 data
# Authors:      Morten W. Hansen, Knut-Frode Dagestad, Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from __future__ import unicode_literals, division, absolute_import
import os
import tarfile
import zipfile
from dateutil.parser import parse
from math import asin

import json
import numpy as np
try:
    import scipy.ndimage
except:
    IMPORT_SCIPY = False
else:
    IMPORT_SCIPY = True

import pythesint as pti

from nansat.nsr import NSR
from nansat.vrt import VRT
from nansat.domain import Domain
from nansat.node import Node
from nansat.utils import initial_bearing, gdal, ogr
from nansat.exceptions import WrongMapperError, NansatReadError


class Mapper(VRT):
    ''' Create VRT with mapping of WKV for Radarsat2 '''

    def __init__(self, inputFileName, gdalDataset, gdalMetadata,
                 xmlonly=False,  **kwargs):
        ''' Create Radarsat2 VRT '''
        fPathName, fExt = os.path.splitext(inputFileName)

        if zipfile.is_zipfile(inputFileName):
            # Open zip file using VSI
            fPath, fName = os.path.split(fPathName)
            filename = '/vsizip/%s/%s' % (inputFileName, fName)
            if not 'RS' in fName[0:2]:
                raise WrongMapperError('%s: Provided data is not Radarsat-2'
                        %fName)
            gdalDataset = gdal.Open(filename)
            gdalMetadata = gdalDataset.GetMetadata()
        else:
            filename = inputFileName

        #if it is not RADARSAT-2, return
        if (not gdalMetadata or not 'SATELLITE_IDENTIFIER' in list(gdalMetadata.keys())):
            raise WrongMapperError(filename)
        elif gdalMetadata['SATELLITE_IDENTIFIER'] != 'RADARSAT-2':
            raise WrongMapperError(filename)

        if zipfile.is_zipfile(inputFileName):
            # Open product.xml to get additional metadata
            zz = zipfile.ZipFile(inputFileName)
            productXmlName = os.path.join(os.path.basename(inputFileName).split('.')[0],'product.xml')
            productXml = zz.open(productXmlName).read()
        else:
            # product.xml to get additionali metadata
            productXmlName = os.path.join(filename,'product.xml')
            if not os.path.isfile(productXmlName):
                raise WrongMapperError(filename)
            productXml = open(productXmlName).read()

        if not IMPORT_SCIPY:
            raise NansatReadError('Radarsat-2 data cannot be read because scipy is not installed')

        # parse product.XML
        rs2_0 = Node.create(productXml)

        if xmlonly:
            self.init_from_xml(rs2_0)
            return

        # Get additional metadata from product.xml
        rs2_1 = rs2_0.node('sourceAttributes')
        rs2_2 = rs2_1.node('radarParameters')
        if rs2_2['antennaPointing'].lower() == 'right':
            antennaPointing = 90
        else:
            antennaPointing = -90
        rs2_3 = rs2_1.node('orbitAndAttitude').node('orbitInformation')
        passDirection = rs2_3['passDirection']

        # create empty VRT dataset with geolocation only
        self._init_from_gdal_dataset(gdalDataset)

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
                    # The nansat data will be complex
                    # if the SAR data is of type 10
                    dtype = iBand.DataType
                    if dtype == 10:
                        # add intensity band
                        metaDict.append(
                            {'src': {'SourceFilename':
                                     ('RADARSAT_2_CALIB:SIGMA0:'
                                      + filename + '/product.xml'),
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
                                                    + filename
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

        '''
        (GDAL?) Radarsat-2 data is stored with maximum latitude at first
        element of each column and minimum longitude at first element of each
        row (e.g. np.shape(lat)=(59,55) -> latitude maxima are at lat[0,:],
        and longitude minima are at lon[:,0])

        In addition, there is an interpolation error for direct estimate along
        azimuth. We therefore estimate the heading along range and add 90
        degrees to get the "satellite" heading.

        '''
        if str(passDirection).upper() == 'DESCENDING':
            sat_heading = initial_bearing(lon[:, :-1], lat[:, :-1],
                                          lon[:, 1:], lat[:, 1:]) + 90
        elif str(passDirection).upper() == 'ASCENDING':
            sat_heading = initial_bearing(lon[:, 1:], lat[:, 1:],
                                          lon[:, :-1], lat[:, :-1]) + 90
        else:
            print('Can not decode pass direction: ' + str(passDirection))

        # Calculate SAR look direction
        look_direction = sat_heading + antennaPointing
        # Interpolate to regain lost row
        look_direction = np.mod(look_direction, 360)
        look_direction = scipy.ndimage.interpolation.zoom(
            look_direction, (1, 11./10.))
        # Decompose, to avoid interpolation errors around 0 <-> 360
        look_direction_u = np.sin(np.deg2rad(look_direction))
        look_direction_v = np.cos(np.deg2rad(look_direction))
        look_u_VRT = VRT.from_array(look_direction_u)
        look_v_VRT = VRT.from_array(look_direction_v)

        # Note: If incidence angle and look direction are stored in
        #       same VRT, access time is about twice as large
        lookVRT = VRT.from_lonlat(lon, lat)
        lookVRT.create_band(
            [{'SourceFilename': look_u_VRT.filename, 'SourceBand': 1},
             {'SourceFilename': look_v_VRT.filename, 'SourceBand': 1}],
            {'PixelFunctionType': 'UVToDirectionTo'})

        # Blow up to full size
        lookVRT = lookVRT.get_resized_vrt(gdalDataset.RasterXSize, gdalDataset.RasterYSize)
        # Store VRTs so that they are accessible later
        self.band_vrts['look_u_VRT'] = look_u_VRT
        self.band_vrts['look_v_VRT'] = look_v_VRT
        self.band_vrts['lookVRT'] = lookVRT

        # Add band to full sized VRT
        lookFileName = self.band_vrts['lookVRT'].filename
        metaDict.append({'src': {'SourceFilename': lookFileName,
                                 'SourceBand': 1},
                         'dst': {'wkv': 'sensor_azimuth_angle',
                                 'name': 'look_direction'}})

        ###############################
        # Create bands
        ###############################
        self.create_bands(metaDict)

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

        self.create_band(src, dst)
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
            self.create_band(src, dst)
            self.dataset.FlushCache()

        ############################################
        # Add SAR metadata
        ############################################
        if antennaPointing == 90:
            self.dataset.SetMetadataItem('ANTENNA_POINTING', 'RIGHT')
        if antennaPointing == -90:
            self.dataset.SetMetadataItem('ANTENNA_POINTING', 'LEFT')
        self.dataset.SetMetadataItem('ORBIT_DIRECTION',
                                     str(passDirection).upper())

        # set valid time
        self.dataset.SetMetadataItem('time_coverage_start',
                                     (parse(gdalMetadata['FIRST_LINE_TIME']).
                                      isoformat()))
        self.dataset.SetMetadataItem('time_coverage_end',
                                     (parse(gdalMetadata['LAST_LINE_TIME']).
                                      isoformat()))

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('sar')
        ee = pti.get_gcmd_platform('radarsat-2')

        # TODO: Validate that the found instrument and platform are indeed what we
        # want....

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

    def init_from_xml(self, productXml):
        ''' Fast init from metada in XML only '''
        numberOfLines = int(productXml
                        .node('imageAttributes')
                        .node('rasterAttributes')
                        .node('numberOfLines')
                        .value)
        numberOfSamples = int(productXml
                        .node('imageAttributes')
                        .node('rasterAttributes')
                        .node('numberOfSamplesPerLine')
                        .value)

        VRT.__init__(self, srcRasterXSize=numberOfSamples, srcRasterYSize=numberOfLines)

        gcps = []
        geogrid = productXml.node('imageAttributes').node('geographicInformation').node('geolocationGrid')
        for child in geogrid.children:
            pix = float(child.node('imageCoordinate').node('pixel').value)
            lin = float(child.node('imageCoordinate').node('line').value)
            lon = float(child.node('geodeticCoordinate').node('longitude').value)
            lat = float(child.node('geodeticCoordinate').node('latitude').value)
            gcps.append(gdal.GCP(lon, lat, 0, pix, lin))

        self.dataset.SetGCPs(gcps, NSR().wkt)

        dates = list(map(parse, [child.node('timeStamp').value for child in
                             (productXml.node('sourceAttributes')
                                        .node('orbitAndAttitude')
                                        .node('orbitInformation')
                                        .nodeList('stateVector'))]))

        self.dataset.SetMetadataItem('time_coverage_start', min(dates).isoformat())
        self.dataset.SetMetadataItem('time_coverage_end', max(dates).isoformat())
        self.dataset.SetMetadataItem('platform', json.dumps(pti.get_gcmd_platform('radarsat-2')))
        self.dataset.SetMetadataItem('instrument', json.dumps(pti.get_gcmd_instrument('SAR')))
        self.dataset.SetMetadataItem('Entry Title', 'Radarsat-2 SAR')
        self.dataset.SetMetadataItem('Data Center', 'CSA')
        self.dataset.SetMetadataItem('ISO Topic Category', 'Oceans')
        self.dataset.SetMetadataItem('Summary', 'Radarsat-2 SAR data')
