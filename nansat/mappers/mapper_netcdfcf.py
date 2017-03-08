from nansat.vrt import VRT
from nansat.tools import WrongMapperError

class Mapper(VRT):

    def __init__(filename, gdalDataset, gdalMetadata, *args, **kwargs):

        if not gdalMetadata:
            raise WrongMapperError

    def add_band_at_time(self, band, time):
        pass
