# Name:    tools.py
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
from nansat.warnings import NansatFutureWarning
import logging
from dateutil.parser import parse

try:
    if 'DISPLAY' not in os.environ:
        import matplotlib; matplotlib.use('Agg')
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from matplotlib.colors import hex2color
    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Polygon
except ImportError:
    BASEMAP_LIB_IS_INSTALLED = False
    MATPLOTLIB_IS_INSTALLED = False

else:
    BASEMAP_LIB_IS_INSTALLED = True
    MATPLOTLIB_IS_INSTALLED = True

import numpy as np

try:
    import gdal, ogr, osr
except:
    from osgeo import gdal, ogr, osr
gdal.UseExceptions()

EXCEPTION_WARNING = (
                             'Use Nansat(filename).')

class OptionError(Exception):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            'nansat.tools.OptionError will be disabled in Nansat 1.1. Use ValueError instead.',
            NansatFutureWarning)


class ProjectionError(Exception):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            'nansat.tools.ProjectionError will be disabled in Nansat 1.1. Use ' \
                'nansat.exceptions.NansatProjectionError instead.',
            NansatFutureWarning)


class GDALError(Exception):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            'nansat.tools.GDALError will be disabled in Nansat 1.1. Use ' \
                'nansat.exceptions.NansatGDALError instead.',
            NansatFutureWarning)


class NansatReadError(Exception):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            'nansat.tools.NansatReadError will be disabled in Nansat 1.1. Use ' \
                'nansat.exceptions.NansatReadError instead.',
        NansatFutureWarning)


class GeolocationError(Exception):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            'nansat.tools.GeolocationError will be disabled in Nansat 1.1. Use ' \
                'nansat.exceptions.NansatGeolocationError instead.',
            NansatFutureWarning)


class WrongMapperError(Exception):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            'nansat.tools.WrongMapperError will be disabled in Nansat 1.1. Use ' \
                'nansat.exceptions.WrongMapperError instead.',
        NansatFutureWarning)



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
        cm.register_cmap(name='obpg', data=obpg, lut=256)
        cm.register_cmap(name='ak01', data=ak01, lut=256)

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


def get_random_color(c0=None, minDist=100, low=0, high=255):
    ''' Create random color which is far enough from the input color

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
    '''
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


def test_openable(fname):
    try:
        f = open(fname, 'r')
    except IOError:
        raise
    f.close()


def write_domain_map(border, out_filename, lon_vec=None, lat_vec=None, lon_border=10.,
                     lat_border=10.,figure_size=(6, 6), dpi=50, projection='cyl', resolution='c',
                     continets_color='coral',meridians=10, parallels=10, p_color='r', p_line='k',
                     p_alpha=0.5, padding=0., mer_labels=[False, False, False, False],
                     par_labels=[False, False, False, False], pltshow=False, labels=None):
        """Create an image with a map of the domain

        Uses Basemap to create a World Map
        Adds a semitransparent patch with outline of the Domain
        Writes to an image file

        Parameters
        -----------
        out_filename : string
            name of the output file name
        lon_vec : [floats] or [[floats]]
            longitudes of patches to display
        lat_vec : [floats] or [[floats]]
            latitudes of patches to display
        lon_border : float
            10, horisontal border around patch (degrees of longitude)
        lat_border : float
            10, vertical border around patch (degrees of latitude)
        figure_size : tuple of two integers
            (6, 6), size of the generated figure in inches
        dpi: int
            50, resolution of the output figure (size 6,6 and dpi 50
            produces 300 x 300 figure)
        projection : string, one of Basemap projections
            'cyl', projection of the map
        resolution : string, resolution of the map
            'c', crude
            'l', low
            'i', intermediate
            'h', high
            'f', full
        continets_color : string or any matplotlib color representation
            'coral', color of continets
        meridians : int
            10, number of meridians to draw
        parallels : int
            10, number of parallels to draw
        p_color : string or any matplotlib color representation
            'r', color of the Domain patch
        p_line : string or any matplotlib color representation
            'k', color of the Domain outline
        p_alpha : float 0 - 1
            0.5, transparency of Domain patch
        padding : float
            0., width of white padding around the map
        mer_labels : list of 4 booleans
            where to put meridian labels, see also Basemap.drawmeridians()
        par_lables : list of 4 booleans
            where to put parallel labels, see also Basemap.drawparallels()
        labels : list of str
            labels to print on top of patches
        """
        if not BASEMAP_LIB_IS_INSTALLED:
            raise ImportError(' Basemap is not installed. Cannot use Domain.write_map. '
                              ' Enable by: conda install -c conda forge basemap ')

        # if lat/lon vectors are not given as input
        if lon_vec is None or lat_vec is None or len(lon_vec) != len(lat_vec):
            lon_vec, lat_vec = border

        # convert vectors to numpy arrays
        lon_vec = np.array(lon_vec)
        lat_vec = np.array(lat_vec)

        # estimate mean/min/max values of lat/lon of the shown area
        # (real lat min max +/- latBorder) and (real lon min max +/- lonBorder)
        min_lon = max(-180, lon_vec.min() - lon_border)
        max_lon = min(180, lon_vec.max() + lon_border)
        min_lat = max(-90, lat_vec.min() - lat_border)
        max_lat = min(90, lat_vec.max() + lat_border)
        mean_lon = lon_vec.mean()
        mean_lat = lat_vec.mean()

        # generate template map (can be also tmerc)
        plt.figure(num=1, figsize=figure_size, dpi=dpi)
        bmap = Basemap(projection=projection,
                       lat_0=mean_lat, lon_0=mean_lon,
                       llcrnrlon=min_lon, llcrnrlat=min_lat,
                       urcrnrlon=max_lon, urcrnrlat=max_lat,
                       resolution=resolution)

        # add content: coastline, continents, meridians, parallels
        bmap.drawcoastlines()
        bmap.fillcontinents(color=continets_color)
        bmap.drawmeridians(np.linspace(min_lon, max_lon, meridians),
                           labels=mer_labels, fmt='%2.1f')
        bmap.drawparallels(np.linspace(min_lat, max_lat, parallels),
                           labels=par_labels, fmt='%2.1f')

        # convert input lat/lon vectors to arrays of vectors with one row
        # if only one vector was given
        if len(lon_vec.shape) == 1:
            lon_vec = [lon_vec]
            lat_vec = [lat_vec]

        for i in range(len(lon_vec)):
            # convert lat/lons to map units
            map_x, map_y = bmap(list(lon_vec[i].flat), list(lat_vec[i].flat))

            # from x/y vectors create a Patch to be added to map
            boundary = Polygon(list(zip(map_x, map_y)),
                               alpha=p_alpha, ec=p_line, fc=p_color)

            # add patch to the map
            plt.gca().add_patch(boundary)
            plt.gca().set_aspect('auto')

            if labels is not None and labels[i] is not None:
                plt.text(np.mean(map_x), np.mean(map_y), labels[i],
                         va='center', ha='right', alpha=0.5, fontsize=10)

        # save figure and close
        plt.savefig(out_filename, bbox_inches='tight',
                    dpi=dpi, pad_inches=padding)
        if pltshow:
            plt.show()
        else:
            plt.close('all')

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
