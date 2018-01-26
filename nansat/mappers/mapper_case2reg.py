# Name:         mapper_case2reg.py
# Purpose:      Mapper for the BEAM/Visat output of Case2Regional algorithm
# Authors:      Asuka Yamakava, Anton Korosov, Morten Wergeland Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os.path
import numpy as np

from nansat.vrt import VRT
from nansat.node import Node
from nansat.nsr import NSR
from nansat.exceptions import WrongMapperError

import pythesint as pti

from nansat.mappers import mapper_generic as mg

class Mapper(mg.Mapper):
    '''Mapping for the BEAM/Visat output of Case2Regional algorithm'''
    def __init__(self, filename, gdalDataset, gdalMetadata,
                 wavelengths=[None, 413, 443, 490, 510, 560, 620, 665,
                              681, 709, 753, None, 778, 864], **kwargs):

        fPathName, fExt = os.path.splitext(filename)
        fPath, fName = os.path.split(fPathName)
        if fExt != '.nc' or 'MER_' not in fName or 'N1_C2IOP' not in fName:
            raise WrongMapperError

        # get all metadata using the GENERIC Mapper
        mg.Mapper.__init__(self, filename, gdalDataset, gdalMetadata)

        #add metadata for Rrs bands
        rrsDict = pti.get_wkv_variable('surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air')

        for bi in range(1, 1+self.dataset.RasterCount):
            b = self.dataset.GetRasterBand(bi)
            bMetadata = b.GetMetadata()
            rawName = bMetadata.get('name', '')
            if 'reflec_' in rawName:
                refNumber = int(rawName.split('_')[1])
                wavelength = wavelengths[refNumber]
                b.SetMetadataItem('name', 'Rrs_' + str(wavelength))
                b.SetMetadataItem('wavelength', str(wavelength))
                for rrsKey in rrsDict:
                    b.SetMetadata(rrsKey, rrsDict[rrsKey])

                src = [{
                    'SourceFilename': b.GetMetadataItem('SourceFilename'),
                    'SourceBand':  b.GetMetadataItem('SourceBand'),
                    'DataType': 6}]
                dst = {
                    'wkv': 'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_water',
                    'suffix': str(wavelength),
                    'wavelength': str(wavelength),
                    'PixelFunctionType': 'NormReflectanceToRemSensReflectance'}
                self.create_band(src, dst)

        self.dataset.FlushCache()
