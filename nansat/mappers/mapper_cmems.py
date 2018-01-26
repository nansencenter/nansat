import datetime
from dateutil.parser import parse

import json
import pythesint as pti

from nansat.vrt import VRT
from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF
from nansat.exceptions import WrongMapperError


# Dictionary with mappings of GCMD metadata keywords:
def get_gcmd_keywords_mapping():
    gcmd_keywords_mapping = {
        'MERCATOR PSY4QV3R1': {
            'instrument': json.dumps(pti.get_gcmd_instrument('computer')),
            'platform': json.dumps(pti.get_gcmd_platform('models')),
        },
        'NERSC-HYCOM model fields': {
            'instrument': json.dumps(pti.get_gcmd_instrument('computer')),
            'platform': json.dumps(pti.get_gcmd_platform('models')),
        },
        'MERCATOR BIOMER4V1R2': {
            'instrument': json.dumps(pti.get_gcmd_instrument('computer')),
            'platform': json.dumps(pti.get_gcmd_platform('models')),
        },
    }
    return gcmd_keywords_mapping

class Mapper(NetcdfCF):

    def __init__(self, *args, **kwargs):

        filename = args[0]
        gdal_metadata = VRT._remove_strings_in_metadata_keys(args[2],
                                                            ['NC_GLOBAL#', 'NANSAT_', 'GDAL_'])

        gcmd_keywords_mapping = get_gcmd_keywords_mapping()
        for key, val in list(gcmd_keywords_mapping.items()):
            if 'source' in list(gdal_metadata.keys()) and key in gdal_metadata['source']:
                instrument = gcmd_keywords_mapping[key]['instrument']
                platform = gcmd_keywords_mapping[key]['platform']

        if not 'instrument' in locals():
            raise WrongMapperError

        super(Mapper, self).__init__(*args, **kwargs)

        time_coverage_start, time_coverage_end = self.time_coverage()

        self.dataset.SetMetadataItem('time_coverage_start',
                (time_coverage_start.isoformat()))
        self.dataset.SetMetadataItem('time_coverage_end',
                (time_coverage_end.isoformat()))
        self.dataset.SetMetadataItem('instrument', instrument)
        self.dataset.SetMetadataItem('platform', platform)

    def time_coverage(self):
        times = self.times()
        return times[0].astype(datetime.datetime), \
                times[-1].astype(datetime.datetime)
