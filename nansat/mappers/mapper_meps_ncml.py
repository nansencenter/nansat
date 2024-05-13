import os
import netCDF4
import datetime

import numpy as np

from nansat.exceptions import WrongMapperError
from nansat.mappers.mapper_meps import Mapper as NCMapper

class Mapper(NCMapper):


    def __init__(self, ncml_url, gdal_dataset, gdal_metadata, netcdf_dim=None, *args, **kwargs):

        if not ncml_url.endswith(".ncml"):
            raise WrongMapperError

        dt = 0
        if netcdf_dim is not None and "time" in netcdf_dim.keys():
            ds = netCDF4.Dataset(ncml_url)
            time = netcdf_dim["time"]
            dt = time - np.datetime64(
                datetime.datetime.fromisoformat(ds.time_coverage_start.replace("Z", "+00:00")))
            import ipdb
            ipdb.set_trace()
        url = self._get_odap_url(ncml_url, np.round(dt))

        super(Mapper, self).__init__(url, gdal_dataset, gdal_metadata, *args, **kwargs)


    def _get_odap_url(self, fn, file_num=0):
        """ Get the opendap url to file number 'file_num'. The
        default file number is 0, and yields the forecast time.
        """
        url = (
                "" + os.path.split(fn)[0] + "/member_%02d"
                "/meps_" + os.path.basename(fn).split("_")[2] +
                "_%02d_" + os.path.basename(fn).split("_")[3][:-2]
            ) % (int(os.path.basename(fn)[8:11]), file_num)
        return url
