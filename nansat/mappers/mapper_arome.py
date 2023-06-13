from dateutil.parser import parse

import json
import pythesint as pti

from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF
from nansat.exceptions import WrongMapperError

class Mapper(NetcdfCF):

    def __init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs):

        if not gdal_metadata:
            raise WrongMapperError
        if 'source' not in list(gdal_metadata.keys()):
            raise WrongMapperError
        if not 'arome' in gdal_metadata['title'].lower() and \
                not 'meps' in gdal_metadata['title'].lower():
            raise WrongMapperError
    
        super(Mapper, self).__init__(filename, gdal_dataset, gdal_metadata, *args, **kwargs)

        #self.dataset.SetMetadataItem('time_coverage_start', parse(
        #    gdal_metadata['NC_GLOBAL#min_time'], ignoretz=True, fuzzy=True).isoformat())
        #self.dataset.SetMetadataItem('time_coverage_end', parse(
        #    gdal_metadata['NC_GLOBAL#max_time'], ignoretz=True, fuzzy=True).isoformat()))

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('computer')
        ee = pti.get_gcmd_platform('models')

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))
