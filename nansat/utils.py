# Name:    utils.py
# Purpose: collection of data and funcs used in NANSAT modules
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
from __future__ import absolute_import

import os
import warnings
import logging
from dateutil.parser import parse

try:
    import matplotlib
except ImportError:
    MATPLOTLIB_IS_INSTALLED = False
else:
    MATPLOTLIB_IS_INSTALLED = True
    if 'DISPLAY' not in os.environ:
        matplotlib.use('Agg')
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import hex2color
    from matplotlib import cm

import numpy as np

try:
    import gdal, ogr, osr
except:
    from osgeo import gdal, ogr, osr
gdal.UseExceptions()

NUMPY_TO_GDAL_TYPE_MAP = {
        'uint8': 1,
        'int8': 1,
        'uint16': 2,
        'int16': 3,
        'uint32': 4,
        'int32': 5,
        'float32': 6,
        'float64': 7,
        'complex64': 10,
        'complex128': 11
    }

def remove_keys(dict, keys):
    if keys is None:
        keys = []
    for key in keys:
        dict.pop(key, None)
    return dict

def register_colormaps():
    ''' Create custom colormaps and register them '''
    obpg = {'red': [(0.00, 0.56, 0.56),
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

            'blue': [(0.00, 0.43, 0.43),
                     (0.19, 1.00, 1.00),
                     (0.38, 1.00, 1.00),
                     (0.50, 0.00, 0.00),
                     (0.63, 0.00, 0.00),
                     (0.88, 0.00, 0.00),
                     (1.00, 0.00, 0.00)],
            }


    ak01 = {'red': [(0, 0.1, 0.1,),
                    (0.1, 0.56, 0.56,),
                    (0.22, 0, 0,),
                    (0.27, 0, 0,),
                    (0.37, 0.3, 0.3,),
                    (0.47, 0, 0,),
                    (0.52, 0, 0,),
                    (0.64, 1, 1,),
                    (0.76, 1, 1,),
                    (0.88, 0.4, 0.4,),
                    (1, 1, 1,)],

            'green': [(0, 0, 0,),
                      (0.1, 0, 0,),
                      (0.22, 0, 0,),
                      (0.27, 0, 0,),
                      (0.37, 0.6, 0.6,),
                      (0.47, 0.6, 0.6,),
                      (0.52, 1, 1,),
                      (0.64, 1, 1,),
                      (0.76, 0, 0,),
                      (0.88, 0, 0,),
                      (1, 0.5, 0.5,)],

            'blue': [(0, 0.1, 0.1,),
                     (0.1, 0.5, 0.5,),
                     (0.22, 0.5, 0.5,),
                     (0.27, 1, 1,),
                     (0.37, 1, 1,),
                     (0.47, 0, 0,),
                     (0.52, 0, 0,),
                     (0.64, 0, 0,),
                     (0.76, 0, 0,),
                     (0.88, 0, 0,),
                     (1, 0.5, 0.5,)],
            }

    if MATPLOTLIB_IS_INSTALLED:
        cm.register_cmap(cmap=LinearSegmentedColormap('obpg', obpg, 256))
        cm.register_cmap(cmap=LinearSegmentedColormap('ak01', ak01, 256))

def initial_bearing(lon1, lat1, lon2, lat2):
        """Initial bearing when traversing from point1 (lon1, lat1)
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

        """
        rlon1 = np.radians(lon1)
        rlat1 = np.radians(lat1)
        rlon2 = np.radians(lon2)
        rlat2 = np.radians(lat2)
        bearing = np.arctan2(np.sin(rlon2 - rlon1) * np.cos(rlat2),
                             np.cos(rlat1) * np.sin(rlat2) -
                             np.sin(rlat1) * np.cos(rlat2) *
                             np.cos(rlon2 - rlon1))
        return np.mod(np.degrees(bearing) + 360, 360)


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
    """ Creates and returns logger with default formatting for Nansat

    Parameters
    -----------
    logName : string, optional
        Name of the logger

    Returns
    --------
    logging.logger

    See also
    --------
    `<http://docs.python.org/howto/logging.html>`_

    """
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


def get_random_color(c0=None, minDist=100, low=0, high=255):
    """Create random color which is far enough from the input color

    Parameters
    ----------
        c0 : str
            hexademical representation of the color (e.g. '#ff0000' for red)
        minDist : int
            minimal distance to input color

    Returns
    -------
        c0 : str
            hexademical representation of the new random color
    """
    if not MATPLOTLIB_IS_INSTALLED:
        raise ImportError('Matplotlib is not installed')

    # check inputs
    if c0 is None:
        c0 = '#000000'
    # convert input color to tuple of R,G,B
    c0rgb = np.array(hex2color(c0))

    # create new random color
    c1rgb = np.array([np.random.randint(low, high),
                      np.random.randint(low, high),
                      np.random.randint(low, high)])

    # calculate distance
    d = np.sum((c0rgb - c1rgb)**2)**0.5

    # if distance is small, create new random color
    if d < minDist:
        c1 = get_random_color(c0, minDist)
    else:
        # convert to HEX code
        c1 = '#%02x%02x%02x' % tuple(c1rgb)

    return c1


def parse_time(time_string):
    ''' Parse time string accounting for possible wrong formatting
    Parameters
    ----------
    time_string : str
        string with date and time
    Returns
    -------
        time_value : datetime object

    '''
    time_string = time_string.strip()
    # To account for datasets on the format YYYY-MM-DDZ which is
    # invalid since it has no time, but a timezone
    try:
        time_value = parse(time_string)
    except ValueError:
        if (len(time_string) == 11 and
                time_string.endswith('Z')):
            time_value = parse(time_string[:10])

    return time_value

register_colormaps()

numpy_to_gdal_type = {
    'uint8': 'Byte',
    'int8': 'Byte',
    'uint16': 'UInt16',
    'int16': 'Int16',
    'uint32': 'UInt32',
    'int32': 'Int32',
    'float32': 'Float32',
    'float64': 'Float64',
    'complex64': 'CFloat32',
    'complex128': 'CFloat64'}

gdal_type_to_offset = {
    'Byte': '1',
    'UInt16': '2',
    'Int16': '2',
    'UInt32': '4',
    'Int32': '4',
    'Float32': '4',
    'Float64': '8',
    'CFloat32': '8',
    'CFloat64': '16'}
