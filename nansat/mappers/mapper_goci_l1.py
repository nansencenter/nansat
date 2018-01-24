# Name:        mapper_goci_l1
# Purpose:     Mapping for GOCI-L1B data
# Authors:     Anton Korosov
# Licence:     This file is part of NANSAT. You can redistribute it or modify
#              under the terms of GNU General Public License, v.3
#              http://www.gnu.org/licenses/gpl-3.0.html

from nansat.tools import gdal, ogr
from nansat.exceptions import WrongMapperError
from nansat.vrt import VRT


class Mapper(VRT):
    ''' VRT with mapping of WKV for MODIS Level 1 (QKM, HKM, 1KM) '''

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):
        ''' Create MODIS_L1 VRT '''

        # raise error in case of not GOCI L1B
        try:
            title = gdalMetadata['HDFEOS_POINTS_Scene_Header_Scene_Title']
        except (TypeError, KeyError):
            raise WrongMapperError
        if not title == 'GOCI Level-1B Data':
            raise WrongMapperError

        # set GOCI projection parameters
        lat_0 = gdalMetadata['HDFEOS_POINTS_Map_Projection_Central_Latitude_(parallel)']
        lon_0 = gdalMetadata['HDFEOS_POINTS_Map_Projection_Central_Longitude_(meridian)']
        rasterXSize = int(gdalMetadata['HDFEOS_POINTS_Scene_Header_number_of_columns'].split(' ')[0])
        rasterYSize = int(gdalMetadata['HDFEOS_POINTS_Scene_Header_number_of_rows'].split(' ')[0])
        proj4 = '+proj=ortho +lat_0=%s +lon_0=%s units=m +ellps=WGS84 +datum=WGS84 +no_defs' % (lat_0, lon_0)
        srs = osr.SpatialReference()
        srs.ImportFromProj4(proj4)
        projection = srs.ExportToWkt()
        geoTransform = (-1391500.0, 500.0, 0.0, 1349500.0, 0.0, -500.0)

        # create empty VRT dataset with georeference only
        self._init_from_dataset_params(rasterXSize, rasterYSize, geoTransform, projection)

        # add bands from subdatasets
        subDatasets = gdalDataset.GetSubDatasets()
        metaDict = []
        wavelengths = ['412', '443', '490', '555', '660', '680', '745', '865']
        for subDSI in range(8):
            metaEntry = {'src': {'SourceFilename': subDatasets[subDSI][0],
                                 'SourceBand': 1,
                                 'ScaleRatio': 1e-6},
                         'dst': {'wkv': 'toa_outgoing_spectral_radiance',
                                 'wavelength': wavelengths[subDSI],
                                 'suffix': wavelengths[subDSI]}}
            metaDict.append(metaEntry)

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)
