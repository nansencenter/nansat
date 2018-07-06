# ------------------------------------------------------------------------------
# Name:     mapper_sentinel1_iw.py
# Purpose:
#
# Author:       Artem Moiseev
#
# Created:  05.07.2018
# Copyright:    (c) NERSC
# License: GPL V3
# ------------------------------------------------------------------------------

import os
import re

from nansat import Domain
from nansat.exceptions import WrongMapperError
from nansat.mappers.mapper_generic import Mapper as GenericMapper


class Mapper(GenericMapper):
    """Create VRT with mapping os Sentinel-1 (A and B) IW mode (S1A_IW)

    Parameters
    ----------
    filename : str
        name of input Sentinel-1 IW file
    gdal_dataset : None
    gdal_metadata : None

    Note
    ----
    Reads data with <generic> mapper and then adds GCPs (Acquired from nansat.Domain
    object generated form lon and lat arrays)
    """

    def __init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs):
        # Check if an input file can be read with this mapper, else raise error
        Mapper.check_input(filename)
        # Read the file with <generic> mapper
        GenericMapper.__init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs)
        dom = self.generate_domain()
        # Get GCPs from the created domain and set them for the input file
        self.dataset.SetGCPs(dom.vrt.dataset.GetGCPs(), dom.vrt.dataset.GetProjection())

    @staticmethod
    def parse_filename(file_src):
        """Parse input file src to extract file name and format

        Parameters
        ----------
        file_src : str
            path of input Sentinel-1 IW file

        Returns
        -------
        filename_base: str
            filename without format (i.e. without '.nc')
        file_format: str
            format of the file (i.e '.nc')
        """
        _, filename = os.path.split(file_src)
        filename_base, file_format = os.path.splitext(filename)
        return filename_base, file_format

    @staticmethod
    def check_input(file_src):
        """Check that input file (according to file_src) can be read by the mapper

        Parameters
        ----------
        file_src : str
            path of input Sentinel-1 IW file
        """
        # Extract filename and file format for the src
        filename_base, file_format = Mapper.parse_filename(file_src.strip())

        if file_format != '.nc':
            raise WrongMapperError('Only netCDF format is supported')
        elif not re.fullmatch(r'S1._IW_RVL.*', filename_base):
            raise WrongMapperError()

    def generate_domain(self):
        """Generate nansat.Domain object based on longitudes and latitudes form
        the input file

        Returns
        -------
        dom: nansat.Domain
            Domain generated from longitudes and latitudes
        """
        lon = self.dataset.GetRasterBand(1).ReadAsArray()
        lat = self.dataset.GetRasterBand(2).ReadAsArray()
        dom = Domain(lon=lon, lat=lat)
        return dom
