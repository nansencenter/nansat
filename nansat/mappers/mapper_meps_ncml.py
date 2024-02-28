import os

from nansat.mappers.mapper_meps import Mapper as NCMapper
from nansat.exceptions import WrongMapperError

class Mapper(NCMapper):


    def __init__(self, ncml_url, gdal_dataset, gdal_metadata, file_num=0, *args, **kwargs):

        if not ncml_url.endswith(".ncml"):
            raise WrongMapperError

        url = self._get_odap_url(ncml_url, file_num)
        # Dette fungerer ikke fordi gdal ikke finner gds.GetSubDatasets() = []
        # Men ncml-filen inneholder lustre-stier, så vi kan evt bruke de direkte,
        # og lese som netcdf på ppi - men dette blir en workaround for MET. Opendap burde fungere...

        super(Mapper, self).__init__(url, gdal_dataset, gdal_metadata, *args, **kwargs)


    def _get_odap_url(self, fn, file_num):
        url = os.path.split(fn)[0] + "/member_%02d" % int(os.path.basename(fn)[8:11]) + "/meps_" + os.path.basename(fn).split("_")[2] + "_%02d_" % file_num + os.path.basename(fn).split("_")[3][:-2]
        return url
