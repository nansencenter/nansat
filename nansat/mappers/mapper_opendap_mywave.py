from nansat.mappers.opendap import Opendap
from nansat.exceptions import WrongMapperError
from nansat.nsr import NSR
from netCDF4 import Dataset
import numpy as np


class Mapper(Opendap):

    baseURLs = ['http://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive']
    timeVarName = 'time'
    xName = 'rlon'
    yName = 'rlat'
    timeCalendarStart = '1970-01-01'

    def __init__(self, filename, gdal_dataset, gdal_metadata, date=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        timestamp = date if date else self.get_date(filename)
        ds = Dataset(filename)
        proj4_str = Mapper.assemble_proj4_str(ds.variables['projection_3'])
        try:
            self.srcDSProjection = NSR(proj4_str)
        except KeyError:
            raise WrongMapperError

        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)

    @staticmethod
    def assemble_proj4_str(ds_proj_var):
        """ Generate a GDAL accepted proj4 string """
        proj4_pattern = "+proj=ob_tran +o_proj=longlat +lon_0=%s +o_lat_p=%s +a=6367470 " \
                        "+b=6367470 +to_meter=0.0174532925199 +wktext"
        ds_proj = proj4_pattern % (ds_proj_var.grid_north_pole_longitude - 180,
                                   ds_proj_var.grid_north_pole_latitude)
        return ds_proj

    def convert_dstime_datetimes(self, ds_time):
        """Convert time variable to np.datetime64"""
        ds_datetimes = np.array(
            [(np.datetime64(self.timeCalendarStart).astype('M8[s]')
              + np.timedelta64(int(sec), 's').astype('m8[s]')) for sec in ds_time]).astype('M8[s]')
        return ds_datetimes
