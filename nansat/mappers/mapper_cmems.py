import datetime
from dateutil.parser import parse

import json
import pythesint as pti

from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF
from nansat.tools import WrongMapperError


# Dictionary with mappings of GCMD metadata keywords:
gcmd_keywords_mapping = {
    'MERCATOR PSY4QV3R1': {
        'instrument': json.dumps(pti.get_gcmd_instrument('computer')),
        'platform': json.dumps(pti.get_gcmd_platform('models')),
    },
    'NERSC-HYCOM model fields': {
        'instrument': json.dumps(pti.get_gcmd_instrument('computer')),
        'platform': json.dumps(pti.get_gcmd_platform('models')),
    },
}

class Mapper(NetcdfCF):

    def __init__(self, *args, **kwargs):

        filename = args[0]
        gdal_metadata = self._remove_strings_in_metadata_keys(args[2])

        for key, val in gcmd_keywords_mapping.iteritems():
            if (gdal_metadata.has_key('source') and key in
                    gdal_metadata['source']):
                instrument = gcmd_keywords_mapping[key]['instrument']
                platform = gcmd_keywords_mapping[key]['platform']

        if not instrument:
            raise WrongMapperError

        super(Mapper, self).__init__(*args, **kwargs)

        time_coverage_start, time_coverage_end = self.time_coverage(filename)

        self.dataset.SetMetadataItem('time_coverage_start',
                (time_coverage_start.isoformat()))
        self.dataset.SetMetadataItem('time_coverage_end',
                (time_coverage_end.isoformat()))
        self.dataset.SetMetadataItem('instrument', instrument)
        self.dataset.SetMetadataItem('platform', platform)

    def time_coverage(self, filename):
        times = self.times(filename)
        return times[0].astype(datetime.datetime), \
                times[-1].astype(datetime.datetime)
