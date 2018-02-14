# Name:         exceptions.py
# Purpose:      Definitions of nansat exceptions
# Authors:      Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from __future__ import absolute_import

class NansatProjectionError(Exception):
    """ Cannot get the projection """
    pass


class NansatGDALError(Exception):
    """ Error from GDAL  """
    pass


class NansatReadError(Exception):
    """ Exception if a file cannot be read with Nansat """
    pass

class NansatGeolocationError(Exception):
    """ Exception if geolocation is wrong (e.g., all lat/lon values are 0) """
    pass

class NansatMissingProjectionError(Exception):
    """ Exception raised if no (sub-) dataset has projection """

class WrongMapperError(Exception):
    """ Error for handling data that does not fit a given mapper """
    pass

