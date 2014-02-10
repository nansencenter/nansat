# Name:         mapper_asar.py
# Purpose:      Mapper for Envisat/ASAR data
# Authors:      Asuka Yamakava, Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from vrt import VRT
from envisat import Envisat
from domain import Domain
import numpy as np


class Mapper(VRT, Envisat):
    ''' VRT with mapping of WKV for ASAR Level 1

        See Also
        --------
            http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata,
                 full_incAng=True, geolocation=False, zoomSize=500,
                 step=1, **kwargs):
        '''
        Parameters
        -----------
        fileName : string

        gdalDataset : gdal dataset

        gdalMetadata : gdal metadata

        full_incAng : bool (default is True)
            if True, add full size incedence angle

        geolocation : bool (default is False)
            if True, add gdal geolocation

        zoomSize: int (used in envisat.py)
            size, to which the ADS array will be zoomed using scipy
            array of this size will be stored in memory

        step: int (used in envisat.py)
            step of pixel and line in GeolocationArrays. lat/lon grids are
            generated at that step
        '''

        product = gdalMetadata.get("MPH_PRODUCT")
        if product[0:4] != "ASA_":
            raise AttributeError("ASAR_L1 BAD MAPPER")

        Envisat.__init__(self, fileName, product[0:4])

        # get channel string (remove '/', since NetCDF doesnt support that in metadata)
        polarization = [{'channel' : gdalMetadata['SPH_MDS1_TX_RX_POLAR'].replace("/", ""),
                        'bandNum' : 1}]
        # if there is the 2nd band, get channel string
        if 'SPH_MDS2_TX_RX_POLAR' in gdalMetadata.keys():
            channel = gdalMetadata['SPH_MDS2_TX_RX_POLAR'].replace("/", "")
            if not(channel.isspace()):
                polarization.append({'channel' : channel,
                                     'bandNum' : 2})

        # Create VRTdataset with small VRTRawRasterbands
        self.subVRTs = {'adsVRTs' : self.get_ads_vrts(gdalDataset,
                                         ["first_line_incidence_angle"],
                                         zoomSize=zoomSize, step=step, **kwargs)}

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # get calibration constant
        gotCalibration = True
        try:
            for iPolarization in polarization:
                metaKey = ('MAIN_PROCESSING_PARAMS_ADS_CALIBRATION_FACTORS.%d.EXT_CAL_FACT'
                           %(iPolarization['bandNum']))
                iPolarization['calibrationConst'] = float(
                    gdalDataset.GetMetadataItem(metaKey, 'records'))
        except:
            try:
                for iPolarization in polarization:
                    # Apparently some ASAR files have calibration constant stored in another place
                    metaKey = ('MAIN_PROCESSING_PARAMS_ADS_0_CALIBRATION_FACTORS.%d.EXT_CAL_FACT'
                               %(iPolarization['bandNum']))
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
                             'dst': {'short_name': 'RawCounts_%s'
                                     %iPolarization['channel']}})

        if full_incAng:
            for adsVRT in self.subVRTs['adsVRTs']:
                metaDict.append({
                    'src': {
                        'SourceFilename': adsVRT.fileName,
                        'SourceBand': 1},
                    'dst': {
                        'name': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('name').replace('last_line_', ''),
                        'short_name': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('name').replace('last_line_', ''),
                        'units': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('units')}
                    })
        if gotCalibration:
            for iPolarization in polarization:
                # add dicrtionary for sigma0, ice and water
                short_names = ['sigma0', 'sigma0_normalized_ice',
                               'sigma0_normalized_water']
                wkt = ['surface_backwards_scattering_coefficient_of_radar_wave',
                       'surface_backwards_scattering_coefficient_of_radar_wave_normalized_over_ice',
                       'surface_backwards_scattering_coefficient_of_radar_wave_normalized_over_water']
                sphPass = [gdalMetadata['SPH_PASS'], '', '']

                sourceFileNames = [fileName,
                                   self.subVRTs['adsVRTs'][0].fileName]

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

        # add geolocation arrays

        if geolocation:
            self.add_geolocation_from_ads(gdalDataset,
                                          zoomSize=zoomSize, step=step)

        # Add SAR look direction to metadata domain
        self.dataset.SetMetadataItem('ANTENNA_POINTING', 'RIGHT') # ASAR is always right-looking
        self.dataset.SetMetadataItem('ORBIT_DIRECTION', gdalMetadata['SPH_PASS'].upper())

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
                    sourceFile['ScaleRatio'] = np.sqrt(1.0 / iPolarization['calibrationConst'])
                else:
                    sourceFile['SourceBand'] = 1
                srcFiles.append(sourceFile)
            dst = {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                   'PixelFunctionType': 'Sigma0HHToSigma0VV',
                   'polarization': 'VV',
                   'suffix': 'VV'}
            self._create_band(srcFiles, dst)
            self.dataset.FlushCache()

        # set time
        self._set_envisat_time(gdalMetadata)
