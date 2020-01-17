# Name:         mapper_asar.py
# Purpose:      Mapper for Envisat/ASAR data
# Authors:      Asuka Yamakava, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import numpy as np
from osgeo import gdal
from dateutil.parser import parse
import json
try:
    import scipy.ndimage
except:
    IMPORT_SCIPY = False
else:
    IMPORT_SCIPY = True

import pythesint as pti

from nansat.vrt import VRT
from nansat.mappers.envisat import Envisat
from nansat.domain import Domain
from nansat.utils import initial_bearing
from nansat.exceptions import WrongMapperError, NansatReadError


class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for ASAR Level 1

        See Also
        --------
            http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
    '''

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):

        '''
        Parameters
        -----------
        filename : string

        gdalDataset : gdal dataset

        gdalMetadata : gdal metadata

        '''

        self.setup_ads_parameters(filename, gdalMetadata)

        if self.product[0:4] != "ASA_":
            raise WrongMapperError

        if not IMPORT_SCIPY:
            raise NansatReadError('ASAR data cannot be read because scipy is not installed')

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
        self._init_from_gdal_dataset(gdalDataset, metadata=gdalMetadata)

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
            iBand = gdalDataset.GetRasterBand(iPolarization['bandNum'])
            dtype = iBand.DataType
            shortName = 'RawCounts_%s' %iPolarization['channel']
            bandName = shortName
            dstName = 'raw_counts_%s' % iPolarization['channel']
            if (8 <= dtype and dtype < 12):
                bandName = shortName+'_complex'
                dstName = dstName + '_complex'

            metaDict.append({'src': {'SourceFilename': filename,
                                     'SourceBand': iPolarization['bandNum']},
                             'dst': {'name': dstName}})


            '''
            metaDict.append({'src': {'SourceFilename': filename,
                                     'SourceBand': iPolarization['bandNum']},
                             'dst': {'name': 'raw_counts_%s'
                                     % iPolarization['channel']}})
            '''
            # if raw data is complex, add the intensity band
            if (8 <= dtype and dtype < 12):
                # choose pixelfunction type
                if (dtype == 8 or dtype == 9):
                    pixelFunctionType = 'IntensityInt'
                else:
                    pixelFunctionType = 'intensity'
                # get data type of the intensity band
                intensityDataType = {'8': 3, '9': 4,
                                     '10': 5, '11': 6}.get(str(dtype), 4)
                # add intensity band
                metaDict.append(
                    {'src': {'SourceFilename': filename,
                             'SourceBand': iPolarization['bandNum'],
                             'DataType': dtype},
                     'dst': {'name': 'raw_counts_%s'
                                     % iPolarization['channel'],
                             'PixelFunctionType': pixelFunctionType,
                             'SourceTransferType': gdal.GetDataTypeName(dtype),
                             'dataType': intensityDataType}})

        #####################################################################
        # Add incidence angle and look direction through small VRT objects
        #####################################################################
        lon = self.get_array_from_ADS('first_line_longs')
        lat = self.get_array_from_ADS('first_line_lats')
        inc = self.get_array_from_ADS('first_line_incidence_angle')

        # Calculate SAR look direction (ASAR is always right-looking)
        look_direction = initial_bearing(lon[:, :-1], lat[:, :-1],
                                             lon[:, 1:], lat[:, 1:])
        # Interpolate to regain lost row
        look_direction = scipy.ndimage.interpolation.zoom(
            look_direction, (1, 11./10.))
        # Decompose, to avoid interpolation errors around 0 <-> 360
        look_direction_u = np.sin(np.deg2rad(look_direction))
        look_direction_v = np.cos(np.deg2rad(look_direction))
        look_u_VRT = VRT.from_array(look_direction_u)
        look_v_VRT = VRT.from_array(look_direction_v)

        # Note: If incidence angle and look direction are stored in
        #       same VRT, access time is about twice as large
        incVRT = VRT.from_array(inc)
        lookVRT = VRT.from_lonlat(lon, lat)
        lookVRT.create_band([{'SourceFilename': look_u_VRT.filename,
                               'SourceBand': 1},
                              {'SourceFilename': look_v_VRT.filename,
                               'SourceBand': 1}],
                             {'PixelFunctionType': 'UVToDirectionTo'})

        # Blow up bands to full size
        incVRT = incVRT.get_resized_vrt(gdalDataset.RasterXSize, gdalDataset.RasterYSize)
        lookVRT = lookVRT.get_resized_vrt(gdalDataset.RasterXSize, gdalDataset.RasterYSize)
        # Store VRTs so that they are accessible later
        self.band_vrts = {'incVRT': incVRT,
                        'look_u_VRT': look_u_VRT,
                        'look_v_VRT': look_v_VRT,
                        'lookVRT': lookVRT}

        # Add band to full sized VRT
        incFileName = self.band_vrts['incVRT'].filename
        lookFileName = self.band_vrts['lookVRT'].filename
        metaDict.append({'src': {'SourceFilename': incFileName,
                                 'SourceBand': 1},
                         'dst': {'wkv': 'angle_of_incidence',
                                 'name': 'incidence_angle'}})
        metaDict.append({'src': {'SourceFilename': lookFileName,
                                 'SourceBand': 1},
                         'dst': {'wkv': 'sensor_azimuth_angle',
                                 'name': 'look_direction'}})

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
                sourceFileNames = [filename, incFileName]

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
                            # if ASA_full_incAng,
                            # set 'ScaleRatio' into source file dict
                            sourceFile['ScaleRatio'] = np.sqrt(
                                1.0 / iPolarization['calibrationConst'])
                        else:
                            sourceFile['SourceBand'] = 1
                        srcFiles.append(sourceFile)

                    metaDict.append({
                        'src': srcFiles,
                        'dst': {'short_name': short_names[iPixFunc],
                                'wkv': wkt[iPixFunc],
                                'PixelFunctionType': (
                                pixelFunctionTypes[iPixFunc]),
                                'polarization': iPolarization['channel'],
                                'suffix': iPolarization['channel'],
                                'pass': sphPass[iPixFunc],
                                'dataType': 6}})

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # Add oribit and look information to metadata domain
        # ASAR is always right-looking
        self.dataset.SetMetadataItem('ANTENNA_POINTING', 'RIGHT')
        self.dataset.SetMetadataItem('ORBIT_DIRECTION',
                        gdalMetadata['SPH_PASS'].upper().strip())

        ###################################################################
        # Estimate sigma0_VV from sigma0_HH
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
                    # if ASA_full_incAng,
                    # set 'ScaleRatio' into source file dict
                    sourceFile['ScaleRatio'] = np.sqrt(
                        1.0 / iPolarization['calibrationConst'])
                else:
                    sourceFile['SourceBand'] = 1
                srcFiles.append(sourceFile)
            dst = {'wkv': (
                   'surface_backwards_scattering_coefficient_of_radar_wave'),
                   'PixelFunctionType': 'Sigma0HHToSigma0VV',
                   'polarization': 'VV',
                   'suffix': 'VV'}
            self.create_band(srcFiles, dst)
            self.dataset.FlushCache()

        # set time
        self._set_envisat_time(gdalMetadata)

        # When using TPS for reprojection, use only every 3rd GCP
        # to improve performance (tradeoff vs accuracy)
        self.dataset.SetMetadataItem('skip_gcps', '3')

        self.dataset.SetMetadataItem('time_coverage_start',
                                     (parse(gdalMetadata['MPH_SENSING_START']).
                                      isoformat()))
        self.dataset.SetMetadataItem('time_coverage_end',
                                     (parse(gdalMetadata['MPH_SENSING_STOP']).
                                      isoformat()))

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('asar')
        ee = pti.get_gcmd_platform('envisat')

        # TODO: Validate that the found instrument and platform are indeed what
        # we want....

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
