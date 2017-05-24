from dateutil.parser import parse

import json
import pythesint as pti

from nansat.mappers.mapper_netcdf_cf import Mapper as NetcdfCF
from nansat.tools import WrongMapperError

# Dictionary with mappings of GCMD metadata keywords:
gmcd_keywords = {
        'global-analysis-forecast-phy': ''
    }

class Mapper(NetcdfCF):

    def __init__(self, *args, **kwargs):
        raise WrongMapperError
