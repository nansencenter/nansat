import datetime
from dateutil.parser import parse

import json
import pythesint as pti

from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF
from nansat.tools import WrongMapperError


class Mapper(NetcdfCF):

    # Dictionary with mappings of GCMD metadata keywords:
    gcmd_keywords_mapping = {
        'MERCATOR PSY4QV3R1': {
            'time_coverage_start': time_coverage_start_mercator,
            'time_coverage_end': time_coverage_end_mercator,
        },
        'NERSC-HYCOM model fields': {
            'time_coverage_start': time_coverage_start_nersc,
            'time_coverage_end': time_coverage_end_nersc,
        },
    }

    def __init__(self, *args, **kwargs):

        filename = args[1]
        gdal_metadata = self._remove_strings_in_metadata_keys(args[2])
        gdal_metadata['filename'] = filename

        #import ipdb
        #ipdb.set_trace()
        #raise WrongMapperError
        
        for key, val in self.gcmd_keywords_mapping.iteritems():
            if not (gdal_metadata.has_key('source') and key in
                    gdal_metadata['source']):
                raise WrongMapperError
            else:
                time_coverage_start = \
                        self.gcmd_keywords_mapping[key]['time_coverage_start'](gdal_metadata)
                time_coverage_end = \
                        self.gcmd_keywords_mapping[key]['time_coverage_end'](gdal_metadata)

        super(Mapper, self).__init__(*args, **kwargs)

        self.dataset.SetMetadataItem('time_coverage_start', (time_coverage_start.isoformat()))
        self.dataset.SetMetadataItem('time_coverage_end', (time_coverage_end.isoformat()))

        mm = pti.get_gcmd_instrument('computer')
        self.dataset.SetMetadataItem('instrument', json.dumps(mm))

        ee = pti.get_gcmd_platform('models')
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

    def time_coverage_start_mercator(metadata):
        return parse(metadata['julian_day_unit'], fuzzy=True) + \
                datetime.timedelta(hours=int(metadata['time_min']))
    
    def time_coverage_end_mercator(metadata):
        return parse(metadata['julian_day_unit'], fuzzy=True) + \
                datetime.timedelta(hours=int(metadata['time_max']))
    
    def time_coverage_start_nersc(metadata):
        filename = ''
        times = self.times(metadata['filename'])
    
        import ipdb
        ipdb.set_trace()
        print 'hei'
    
    def time_coverage_end_nersc(metadata):
        filename = ''
        times = self.times(metadata['filename'])
    
