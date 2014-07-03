# Name:         mapper_asar.py
# Purpose:      Mapper for Envisat/ASAR data
# Authors:      Asuka Yamakava, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

import numpy as np
import scipy.ndimage

from nansat.vrt import VRT
from envisat import Envisat
from nansat.domain import Domain
from nansat.tools import initial_bearing


class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for ASAR Level 1

        See Also
        --------
            http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):

        '''
        Parameters
        -----------
        fileName : string

        gdalDataset : gdal dataset

        gdalMetadata : gdal metadata

        '''

        product = gdalMetadata.get("MPH_PRODUCT")
        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER")

        Envisat.__init__(self, fileName, product[0:4])

        # get channel string (remove '/', since NetCDF
        # does not support that in metadata)
        polarization = [{'channel': gdalMetadata['SPH_MDS1_TX_RX_POLAR']
                        .replace("/", ""), 'bandNum': 1}]
        # if there is the 2nd band, get channel string
        if 'SPH_MDS2_TX_RX_POLAR' in gdalMetadata.keys():
            channel = gdalMetadata['SPH_MDS2_TX_RX_POLAR'].replace("/", "")
            if not(channel.isspace()):
                polarization.append({'channel': channel,
                                     'bandNum': 2})

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # get calibration constant
        gotCalibration = True
        try:
            for iPolarization in polarization:
                metaKey = ('MAIN_PROCESSING_PARAMS_ADS_CALIBRATION_FACTORS.%d.EXT_CAL_FACT'
                           % (iPolarization['bandNum']))
                iPolarization['calibrationConst'] = float(
                    gdalDataset.GetMetadataItem(metaKey, 'records'))
        except:
            try:
                for iPolarization in polarization:
                    # Apparently some ASAR files have calibration
                    # constant stored in another place
                    metaKey = ('MAIN_PROCESSING_PARAMS_ADS_0_CALIBRATION_FACTORS.%d.EXT_CAL_FACT'
                               % (iPolarization['bandNum']))
                    iPolarization['calibrationConst'] = float(
                        gdalDataset.GetMetadataItem(metaKey, 'records'))
            except:
                self.logger.warning('Cannot get calibrationConst')
                gotCalibration = False

        # add dictionary for raw counts
        metaDict = []
        for iPolarization in polarization:
            metaDict.append({'src': {'SourceFilename': fileName,
                                     'SourceBand': iPolarization['bandNum']},
                             'dst': {'name': 'raw_counts_%s'
                                     % iPolarization['channel']}})

        #####################################################################
        # Add incidence angle and look direction through small VRT objects
        #####################################################################
        lon = self.get_array_from_ADS('first_line_longs')
        lat = self.get_array_from_ADS('first_line_lats')
        inc = self.get_array_from_ADS('first_line_incidence_angle')

        # Calculate SAR look direction (ASAR is always right-looking)
        SAR_look_direction = initial_bearing(lon[:, :-1], lat[:, :-1],
                                             lon[:, 1:], lat[:, 1:])
        # Interpolate to regain lost row
        SAR_look_direction = scipy.ndimage.interpolation.zoom(
                                SAR_look_direction, (1, 11./10.))
        # Decompose, to avoid interpolation errors around 0 <-> 360
        SAR_look_direction_u = np.sin(np.deg2rad(SAR_look_direction))
        SAR_look_direction_v = np.cos(np.deg2rad(SAR_look_direction))
        look_u_VRT = VRT(array=SAR_look_direction_u, lat=lat, lon=lon)
        look_v_VRT = VRT(array=SAR_look_direction_v, lat=lat, lon=lon)

        # Note: If incidence angle and look direction are stored in
        #       same VRT, access time is about twice as large
        incVRT = VRT(array=inc, lat=lat, lon=lon)
        lookVRT = VRT(lat=lat, lon=lon)
        lookVRT._create_band(
                [{'SourceFilename': look_u_VRT.fileName, 'SourceBand': 1},
                 {'SourceFilename': look_v_VRT.fileName, 'SourceBand': 1}],
                 {'PixelFunctionType': 'UVToDirectionTo'})

        # Blow up bands to full size
        incVRT = incVRT.get_resized_vrt(
                    gdalDataset.RasterXSize, gdalDataset.RasterYSize)
        lookVRT = lookVRT.get_resized_vrt(
                    gdalDataset.RasterXSize, gdalDataset.RasterYSize)
        # Store VRTs so that they are accessible later
        self.subVRTs = {'incVRT': incVRT,
                        'look_u_VRT': look_u_VRT,
                        'look_v_VRT': look_v_VRT,
                        'lookVRT': lookVRT}

        # Add band to full sized VRT
        incFileName = self.subVRTs['incVRT'].fileName
        lookFileName = self.subVRTs['lookVRT'].fileName
        metaDict.append({'src':
                            {'SourceFilename': incFileName,
                             'SourceBand': 1},
                         'dst':
                            {'wkv': 'angle_of_incidence',
                             'name': 'incidence_angle'}})
        metaDict.append({'src':
                            {'SourceFilename': lookFileName,
                             'SourceBand': 1},
                         'dst':
                            {'wkv': 'sensor_azimuth_angle',
                             'name': 'SAR_look_direction'}})

        ####################
        # Add Sigma0-bands
        ####################
        if gotCalibration:
            for iPolarization in polarization:
                # add dictionary for sigma0, ice and water
                short_names = ['sigma0', 'sigma0_normalized_ice',
                               'sigma0_normalized_water']
                wkt = [
                    'surface_backwards_scattering_coefficient_of_radar_wave',
                    'surface_backwards_scattering_coefficient_of_radar_wave_normalized_over_ice',
                    'surface_backwards_scattering_coefficient_of_radar_wave_normalized_over_water']
                sphPass = [gdalMetadata['SPH_PASS'], '', '']
                sourceFileNames = [fileName, incFileName]

                pixelFunctionTypes = ['RawcountsIncidenceToSigma0',
                                      'Sigma0NormalizedIce']
                if iPolarization['channel'] == 'HH':
                    pixelFunctionTypes.append('Sigma0HHNormalizedWater')
                elif iPolarization['channel'] == 'VV':
                    pixelFunctionTypes.append('Sigma0VVNormalizedWater')

                # add pixelfunction bands to metaDict
                for iPixFunc in range(len(pixelFunctionTypes)):
                    srcFiles = []
                    for j, jFileName in enumerate(sourceFileNames):
                        sourceFile = {'SourceFilename': jFileName}
                        if j == 0:
                            sourceFile['SourceBand'] = iPolarization['bandNum']
                            # if ASA_full_incAng, set 'ScaleRatio' into source file dict
                            sourceFile['ScaleRatio'] = np.sqrt(1.0 / iPolarization['calibrationConst'])
                        else:
                            sourceFile['SourceBand'] = 1
                        srcFiles.append(sourceFile)

                    metaDict.append({'src': srcFiles,
                                     'dst': {'short_name': short_names[iPixFunc],
                                             'wkv': wkt[iPixFunc],
                                             'PixelFunctionType': pixelFunctionTypes[iPixFunc],
                                             'polarization': iPolarization['channel'],
                                             'suffix': iPolarization['channel'],
                                             'pass': sphPass[iPixFunc],
                                             'dataType': 6}})

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # Add oribit and look information to metadata domain
        self.dataset.SetMetadataItem('ANTENNA_POINTING', 'RIGHT')  # ASAR is always right-looking
        self.dataset.SetMetadataItem('ORBIT_DIRECTION',
                                     gdalMetadata['SPH_PASS'].upper())

        ###################################################################
        # Add sigma0_VV
        ###################################################################
        polarizations = []
        for pp in polarization:
            polarizations.append(pp['channel'])
        if 'VV' not in polarizations and 'HH' in polarizations:
            srcFiles = []
            for j, jFileName in enumerate(sourceFileNames):
                sourceFile = {'SourceFilename': jFileName}
                if j == 0:
                    sourceFile['SourceBand'] = iPolarization['bandNum']
                    # if ASA_full_incAng, set 'ScaleRatio' into source file dict
                    sourceFile['ScaleRatio'] = np.sqrt(
                        1.0 / iPolarization['calibrationConst'])
                else:
                    sourceFile['SourceBand'] = 1
                srcFiles.append(sourceFile)
            dst = {'wkv':
                      'surface_backwards_scattering_coefficient_of_radar_wave',
                   'PixelFunctionType': 'Sigma0HHToSigma0VV',
                   'polarization': 'VV',
                   'suffix': 'VV'}
            self._create_band(srcFiles, dst)
            self.dataset.FlushCache()

        # set time
        self._set_envisat_time(gdalMetadata)
