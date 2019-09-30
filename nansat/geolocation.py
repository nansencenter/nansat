# Name:    geolocation.py
# Purpose: Container of Geolocation class
# Authors:      Anton Korosov
# Created:      14.01.2018
# Copyright:    (c) NERSC 2011 - 2018
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
from __future__ import absolute_import

import gdal, osr

from nansat.nsr import NSR


class Geolocation(object):
    """Container for GEOLOCATION data

    Keeps references to bands with X and Y coordinates, offset and step
    of pixel and line. All information is stored in dictionary self.data

    Instance of Geolocation is used in VRT and ususaly created in
    a Mapper.
    """
    # instance attributes
    data = None
    x_vrt = None
    y_vrt = None

    def __init__(self, x_vrt, y_vrt, **kwargs):
        """Create Geolocation object from input VRT objects

        Parameters
        -----------
        x_vrt_name : VRT
            dataset with X-coordinates
        y_vrt_name : VRT
            dataset with Y-coordinates
        **kwargs : dict
            parameters for self._init_data()

        Note
        ----
        Modifies, self.x_vrt, self.y_vrt, self.data

        """
        # dictionary with all metadata
        self.data = dict()
        # VRT objects
        self.x_vrt = x_vrt
        self.y_vrt = y_vrt

        self._init_data(x_vrt.filename, y_vrt.filename, **kwargs)

    def _init_data(self, x_filename, y_filename, x_band=1, y_band=1, srs=None,
                                                line_offset=0, line_step=1,
                                                pixel_offset=0, pixel_step=1):
        """Init data of Geolocation object from input parameters
        Parameters
        -----------
        x_filename : str
            name of file for X-dataset
        y_filename : str
            name of file for Y-dataset
        x_band : int
            number of the band in the X-dataset
        y_band : int
            number of the band in the Y-dataset
        srs : str
            WKT
        line_offset : int
            offset of first line
        line_step : int
            step of lines
        pixel_offset : int
            offset of first pixel
        pixel_step : int
            step of pixels

        Notes
        -----
        Saves everything in self.data dict

        """
        if srs is None:
            srs = NSR().wkt
        self.data['SRS'] = srs
        self.data['X_DATASET'] = x_filename
        self.data['Y_DATASET'] = y_filename
        self.data['X_BAND'] = str(x_band)
        self.data['Y_BAND'] = str(y_band)
        self.data['LINE_OFFSET'] = str(line_offset)
        self.data['LINE_STEP'] = str(line_step)
        self.data['PIXEL_OFFSET'] = str(pixel_offset)
        self.data['PIXEL_STEP'] = str(pixel_step)

    @classmethod
    def from_dataset(cls, dataset):
        """Create geolocation from GDAL dataset
        Parameters
        ----------
        dataset : gdal.Dataset
            input dataset to copy Geolocation metadata from
        """
        self = cls.__new__(cls) # empty object
        self.x_vrt = None
        self.y_vrt = None
        self.data = dataset.GetMetadata('GEOLOCATION')
        return self

    @classmethod
    def from_filenames(cls, x_filename, y_filename, **kwargs):
        """Create geolocation from names of files with geolocation
        Parameters
        ----------
        x_filename : str
            name of file for X-dataset
        y_filename : str
            name of file for Y-dataset
        **kwargs : dict
            parameters for self._init_data()
        """
        self = cls.__new__(cls) # empty object
        self.x_vrt = None
        self.y_vrt = None
        self.data = {}
        self._init_data(x_filename, y_filename, **kwargs)
        return self

    def get_geolocation_grids(self):
        """Read values of geolocation grids"""
        lon_dataset = gdal.Open(self.data['X_DATASET'])
        lon_grid = lon_dataset.GetRasterBand(int(self.data['X_BAND'])).ReadAsArray()
        lat_dataset = gdal.Open(self.data['Y_DATASET'])
        lat_grid = lat_dataset.GetRasterBand(int(self.data['Y_BAND'])).ReadAsArray()
        return lon_grid, lat_grid
