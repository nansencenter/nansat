# Name:         mapper_asar.py
# Purpose:      Mapper for Envisat/ASAR data
# Authors:      Asuka Yamakava, Anton Korosov, Knut-Frode Dagestad
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from nansat import Nansat
from nansat.vrt import VRT
from envisat import Envisat
from nansat.domain import Domain
from nansat_tools import initial_bearing
import numpy as np
import scipy.ndimage


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

        # Add incidence angle through small Nansat object
        lon = self.get_array_from_ADS('first_line_longs')
        lat = self.get_array_from_ADS('first_line_lats')
        inc = self.get_array_from_ADS('first_line_incidence_angle')
        ADS = Nansat(domain=Domain(lon=lon, lat=lat), array=inc)

        # Calculate SAR look direction (ASAR is always right-looking)
        SAR_look_direction = initial_bearing(lon[:, :-1], lat[:, :-1],
                                             lon[:, 1:], lat[:, 1:])
        # Interpolate to regain lost row
        SAR_look_direction = scipy.ndimage.interpolation.zoom(
                                SAR_look_direction, (1, 11./10.))
        ADS.add_band(array=SAR_look_direction)

        # Reproject small array onto full SAR domain
        ADS.reproject(Domain(ds=gdalDataset), eResampleAlg=1)
        self.subVRTs = {}
        self.subVRTs['ADS_arrays'] = ADS.vrt
        # Add band to VRT
        ADSFileName = self.subVRTs['ADS_arrays'].fileName
        metaDict.append({'src':
                            {'SourceFilename': ADSFileName,
                             'SourceBand': 1},
                         'dst':
                            {'wkv': 'angle_of_incidence'}})
        metaDict.append({'src':
                            {'SourceFilename': ADSFileName,
                             'SourceBand': 2},
                         'dst':
                            {'wkv': 'sensor_azimuth_angle',
                             'name': 'SAR_look_direction'}})

        if gotCalibration:
            for iPolarization in polarization:
                # add dicrtionary for sigma0, ice and water
                short_names = ['sigma0', 'sigma0_normalized_ice',
                               'sigma0_normalized_water']
                wkt = [
                    'surface_backwards_scattering_coefficient_of_radar_wave',
                    'surface_backwards_scattering_coefficient_of_radar_wave_normalized_over_ice',
                    'surface_backwards_scattering_coefficient_of_radar_wave_normalized_over_water']
                sphPass = [gdalMetadata['SPH_PASS'], '', '']
                sourceFileNames = [fileName, ADSFileName]

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

        # Add SAR look direction to metadata domain
        self.dataset.SetMetadataItem('ANTENNA_POINTING', 'RIGHT')  # ASAR is always right-looking
        self.dataset.SetMetadataItem('ORBIT_DIRECTION',
                                     gdalMetadata['SPH_PASS'].upper())

        # "SAR_center_look_direction" below is obsolete, and may soon be deleted
        #
        # Note that this is the look direction in the center of the domain. For
        # longer domains, especially at high latitudes, the azimuth direction
        # may vary a lot over the domain, and using the center angle will be a
        # coarse approximation.
        self.dataset.SetMetadataItem('SAR_center_look_direction',
                                     str(np.mod(Domain(ds=gdalDataset).
                                         upwards_azimuth_direction() + 90,
                                                360)))

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
