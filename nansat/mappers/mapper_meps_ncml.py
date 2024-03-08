import os

from nansat.mappers.mapper_meps import Mapper as NCMapper
from nansat.exceptions import WrongMapperError

class Mapper(NCMapper):


    def __init__(self, ncml_url, gdal_dataset, gdal_metadata, file_num=0, *args, **kwargs):

        if not ncml_url.endswith(".ncml"):
            raise WrongMapperError

        url = self._get_odap_url(ncml_url, file_num)

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
