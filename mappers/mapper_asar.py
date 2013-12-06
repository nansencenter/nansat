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

        # get polarization string (remove '/', since NetCDF doesnt support that in metadata)
        polarization = gdalMetadata['SPH_MDS1_TX_RX_POLAR'].replace("/", "")

        # Create VRTdataset with small VRTRawRasterbands
        self.subVRTs = {'adsVRTs' : self.get_ads_vrts(gdalDataset,
                                         ["first_line_incidence_angle"],
                                         zoomSize=zoomSize, step=step, **kwargs)}

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset)

        # get calibration constant
        gotCalibration = True
        try:
            calibrationConst = float(gdalDataset.GetMetadataItem(
                "MAIN_PROCESSING_PARAMS_ADS_CALIBRATION_FACTORS.1.EXT_CAL_FACT", "records"))
        except:
            try:
                # Apparently some ASAR files have calibration constant stored in another place
                calibrationConst = float(gdalDataset.GetMetadataItem(
                    "MAIN_PROCESSING_PARAMS_ADS_1_CALIBRATION_FACTORS.1.EXT_CAL_FACT", "records"))
            except:
                self.logger.warning('Cannot get calibrationConst')
                gotCalibration = False

        # add dictionary for raw counts
        metaDict = [{'src': {'SourceFilename': fileName, 'SourceBand': 1},
                     'dst': {'short_name': 'RawCounts'}}]

        if full_incAng:
            for adsVRT in self.subVRTs['adsVRTs']:
                metaDict.append({'src': {'SourceFilename': adsVRT.fileName,
                                         'SourceBand': 1},
                                 'dst': {'name': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('name').replace('last_line_', ''),
                                         'units': adsVRT.dataset.GetRasterBand(1).GetMetadataItem('units')}})
        if gotCalibration:
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
            if polarization == 'HH':
                pixelFunctionTypes.append('Sigma0HHNormalizedWater')
            elif polarization == 'VV':
                pixelFunctionTypes.append('Sigma0VVNormalizedWater')

            # add pixelfunction bands to metaDict
            for iPixFunc in range(len(pixelFunctionTypes)):
                srcFiles = []
                for iFileName in sourceFileNames:
                    sourceFile = {'SourceFilename': iFileName,
                                  'SourceBand': 1}
                    # if ASA_full_incAng, set 'ScaleRatio' into source file dict
                    if iFileName == fileName:
                        sourceFile['ScaleRatio'] = np.sqrt(1.0/calibrationConst)
                    srcFiles.append(sourceFile)

                metaDict.append({'src': srcFiles,
                                 'dst': {'short_name': short_names[iPixFunc],
                                         'wkv': wkt[iPixFunc],
                                         'PixelFunctionType': pixelFunctionTypes[iPixFunc],
                                         'polarization': polarization,
                                         'suffix': polarization,
                                         'pass': sphPass[iPixFunc],
                                         'dataType': 6}})

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set time
        self._set_envisat_time(gdalMetadata)

        # add geolocation arrays

        if geolocation:
            self.add_geolocation_from_ads(gdalDataset,
                                          zoomSize=zoomSize, step=step)

        # Add SAR look direction to metadata domain
        # Note that this is the look direction in the center of the domain. For
        # longer domains, especially at high latitudes, the azimuth direction
        # may vary a lot over the domain, and using the center angle will be a
        # coarse approximation.
        self.dataset.SetMetadataItem('SAR_center_look_direction',
                                     str(np.mod(Domain(ds=gdalDataset).
                                         upwards_azimuth_direction() + 90,
                                                360)))
