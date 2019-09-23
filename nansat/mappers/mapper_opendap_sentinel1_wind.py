import os
import numpy as np
from dateutil.parser import parse
from datetime import datetime
from netCDF4 import Dataset

from nansat.nsr import NSR
from nansat.mappers.opendap import Opendap

class Mapper(Opendap):

    baseURLs = ['http://thredds.nersc.no/thredds/dodsC/sarvind/SarVind',]
    timeVarName = 'time'
    xName = 'x'
    yName = 'y'

    def __init__(self, filename, gdal_dataset, gdal_metadata, 
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        ds = Dataset(filename)
        self.srcDSProjection = ds.variables['stereographic'].spatial_ref
        self.timeCalendarStart = parse(ds.variables[self.timeVarName].units, fuzzy=True)
        self.create_vrt(filename, gdal_dataset, gdal_metadata, self.get_date(filename), ds, bands,
                cachedir)

    @staticmethod
    def get_date(filename):
        """Extract date and time parameters from filename and return
        it as a formatted (isoformat) string

        Parameters
        ----------

        filename: str
            nn

        Returns
        -------
            str, YYYY-mm-ddThh:MMZ

        """
        _, filename = os.path.split(filename)
        t = datetime.strptime(filename.split('_')[5], '%Y%m%dT%H%M%S')
        return datetime.strftime(t, '%Y-%m-%dT%H:%M:%SZ')

    def convert_dstime_datetimes(self, ds_time):
        """ Convert time variable to np.datetime64

        The time unit is days.
        """
        hours = (ds_time[0] - np.floor(ds_time[0]))*24
        minutes = (hours - np.floor(hours))*60
        secs = (minutes - np.floor(minutes))*60
        millisecs = int(np.round((secs - np.floor(secs))*10**3))
        ds_datetimes = np.array(
            [np.datetime64(self.timeCalendarStart) + 
                np.timedelta64(int(np.floor(ds_time[0])), 'D') +
                np.timedelta64(int(np.floor(hours)), 'h') +
                np.timedelta64(int(np.floor(minutes)), 'm') +
                np.timedelta64(int(np.floor(secs)), 's') +
                np.timedelta64(millisecs, 'ms')])
        return ds_datetimes
