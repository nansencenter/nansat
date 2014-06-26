# Name:    nansat_tools.py
# Purpose: collection of libraries used in NANSAT modules
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2013
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
'''
# import standard modules

## used in domain and nansat_tools
import os.path
import re

## used in nansat, vrt, figure
import os

## used in vrt and nansat_tools
import logging

## used in domain
from math import atan2, sin, pi, cos, acos, radians, degrees, copysign
import string
from xml.etree.ElementTree import ElementTree

## used in nansat
import dateutil.parser
import glob
import inspect
import pdb
import sys

## used in vrt
import datetime
from dateutil.parser import parse
from random import choice
from string import Template, ascii_uppercase, digits

## used in figure, nansat_map
from math import floor, log10, pow

## used in nansat_tools
import copy
import warnings
import xml.dom.minidom as xdm

# try to import additional modules
## used in domain, nansat, vrt and nansat_tools
try:
    from osgeo import gdal, osr, ogr
except ImportError:
    try:
        import gdal
        import osr
        import ogr
    except ImportError:
        warnings.warn('Cannot import GDAL!'
                      'Nansat, Vrt, Domain and nansat_tools will not work'
                      'Try installing GDAL.')

## used in nansat, domain, vrt and figure
try:
    import numpy as np
except ImportError:
    warnings.warn('Cannot import numpy!'
                  'Domain, Figure and Vrt will not work.'
                  'Try installing numpy.')

## used in domain and figure
try:
    import matplotlib.pyplot as plt
except:
    warnings.warn('Cannot import matplotlib.pyplot!'
                  'Domain.write_map() and Figure will not work'
                  'Try installing matplotlib.')

## used in nansat and figure
try:
    from matplotlib import cm
except:
    warnings.warn('Cannot import matplotlib.cm!'
                  'Nansat.write_geotiffimage and Figure will not work.'
                  'Try installing matplotlib.')

## used in domain
try:
    from matplotlib.patches import Polygon
except:
    warnings.warn('Cannot import matplotlib.patches.Polygon!'
                  'Domain.write_map() will not work'
                  'Try installing matplotlib.')

try:
    from mpl_toolkits.basemap import Basemap
except:
    warnings.warn('Cannot import mpl_toolkits.basemap.Basemap!'
                  'Domain.write_map() will not work'
                  'Try installing Basemap.')

## used in nansat
try:
    from numpy import arange
except ImportError:
    warnings.warn('Cannot import numpy.arange!'
                  'Nansat.write_geotiffimage will not work.'
                  'Try installing numpy.')

## used in figure
try:
    from numpy import outer
except ImportError:
    warnings.warn('Cannot import numpy.outer!'
                  'Figure.create_legend will not work.'
                  'Try installing numpy.')

try:
    import Image
    import ImageDraw
    import ImageFont
except ImportError:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError:
        warnings.warn('Cannot import PIL!'
                      'Figure will not work'
                      'Try installing PIL.')

## used in nansat_tools and nansatmap
try:
    from scipy import mod, ndimage
except ImportError:
    warnings.warn('Cannot import scipy.mod!'
                  'nansat_toolds will not work'
                  'Try installing scipy.')
'''
obpg = {
'red': [  (0.00, 0.56, 0.56),
          (0.19, 0.00, 0.00),
          (0.38, 0.00, 0.00),
          (0.50, 0.00, 0.00),
          (0.63, 1.00, 1.00),
          (0.88, 1.00, 1.00),
          (1.00, 0.40, 0.40)],
'green': [(0.00, 0.00, 0.00),
          (0.19, 0.00, 0.00),
          (0.38, 1.00, 1.00),
          (0.50, 1.00, 1.00),
          (0.63, 1.00, 1.00),
          (0.88, 0.00, 0.00),
          (1.00, 0.00, 0.00)],
'blue': [ (0.00, 0.43, 0.43),
          (0.19, 1.00, 1.00),
          (0.38, 1.00, 1.00),
          (0.50, 0.00, 0.00),
          (0.63, 0.00, 0.00),
          (0.88, 0.00, 0.00),
          (1.00, 0.00, 0.00)],
}

ak01 = {

'red': [  (0,0.1,0.1,),
(0.1,0.56,0.56,),
(0.22,0,0,),
(0.27,0,0,),
(0.37,0.3,0.3,),
(0.47,0,0,),
(0.52,0,0,),
(0.64,1,1,),
(0.76,1,1,),
(0.88,0.4,0.4,),
(1,1,1,)],


'green': [(0,0,0,),
(0.1,0,0,),
(0.22,0,0,),
(0.27,0,0,),
(0.37,0.6,0.6,),
(0.47,0.6,0.6,),
(0.52,1,1,),
(0.64,1,1,),
(0.76,0,0,),
(0.88,0,0,),
(1,0.5,0.5,)],


'blue': [ (0,0.1,0.1,),
(0.1,0.5,0.5,),
(0.22,0.5,0.5,),
(0.27,1,1,),
(0.37,1,1,),
(0.47,0,0,),
(0.52,0,0,),
(0.64,0,0,),
(0.76,0,0,),
(0.88,0,0,),
(1,0.5,0.5,)],



}
try:
    cm.register_cmap(name='obpg', data=obpg, lut=256)
    cm.register_cmap(name='ak01', data=ak01, lut=256)
except:
    warnings.warn('Cannot generate and register the OBPG colormap!')

class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass


class OptionError(Error):
    '''Error for improper options (arguments) '''
    pass


class ProjectionError(Error):
    '''Cannot get the projection'''
    pass


class GDALError(Error):
    '''Error from GDAL '''
    pass


def initial_bearing(lon1, lat1, lon2, lat2):
        '''Initial bearing when traversing from point1 (lon1, lat1)
        to point2 (lon2, lat2)

        See http://www.movable-type.co.uk/scripts/latlong.html

        Parameters
        ----------
        lon1, lat1 : float
            longitude and latitude of start point
        lon2, lat2 : float
            longitude and latitude of end point

        Returns
        -------
        initial_bearing : float
            The initial bearing (azimuth direction) when heading out
            from the start point towards the end point along a great circle.

        '''
        rlon1 = np.radians(lon1)
        rlat1 = np.radians(lat1)
        rlon2 = np.radians(lon2)
        rlat2 = np.radians(lat2)
        deltalon = rlon2 - rlon1
        bearing = np.arctan2(np.sin(rlon2 - rlon1) * np.cos(rlat2),
                             np.cos(rlat1) * np.sin(rlat2) -
                             np.sin(rlat1) * np.cos(rlat2) *
                             np.cos(rlon2 - rlon1))
        return mod(np.degrees(bearing) + 360, 360)

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the spherical earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    distance_meters = 6367000 * c
    return distance_meters

def add_logger(logName='', logLevel=None):
    ''' Creates and returns logger with default formatting for Nansat

    Parameters
    -----------
    logName : string, optional
        Name of the logger

    Returns
    --------
    logging.logger

    See also
    --------
    http://docs.python.org/howto/logging.html

    '''
    if logLevel is not None:
        os.environ['LOG_LEVEL'] = str(logLevel)
    # create (or take already existing) logger
    # with default logging level WARNING
    logger = logging.getLogger(logName)
    logger.setLevel(int(os.environ['LOG_LEVEL']))

    # if logger already exits, default stream handler has been already added
    # otherwise create and add a new handler
    if len(logger.handlers) == 0:
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        # create formatter
        formatter = logging.Formatter('%(asctime)s|%(levelno)s|%(module)s|'
                                      '%(funcName)s|%(message)s',
                                      datefmt='%I:%M:%S')
        # add formatter to ch
        ch.setFormatter(formatter)
        # add ch to logger
        logger.addHandler(ch)

    logger.handlers[0].setLevel(int(os.environ['LOG_LEVEL']))

    return logger


def set_defaults(dictionary, newParm):
        '''Check input params and set defaut values

        Look throught default parameters (self.d) and given parameters (dict)
        and paste value from input if the key matches

        Parameters
        ----------
        dict : dictionary
            parameter names and values

        Modifies
        ---------
        self.d

        '''
        for key in newParm:
            if key in dictionary:
                dictionary[key] = newParm[key]

        return dictionary
