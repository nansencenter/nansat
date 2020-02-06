# Name:    tools.py
# Purpose: Collection of methods which use core Nansat classes for
#          scientific data analysis and visualization
# Authors:      Artem Moiseev
# Created:      17.01.2020
# Copyright:    (c) NERSC 2020
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

import os
import warnings
import functools

from nansat.nansat import Nansat
from nansat import utils

def distance2coast(dst_domain, distance_src=None):
    """ Estimate distance to the nearest coast (in km) for each pixcel in the
    domain of interest. The method utilizes NASA's OBPG group Distance to the Nearest Coast
    product: https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/. The product is stored in GeoTiff
    format with pixcelsize of 0.01x0.01 degree.

    Parameters
    -----------
    dst_domain : Domain
        destination domain
    distance_src : str
        path to the NASA Distance to the Nearest coast GeoTIFF product

    Returns
    --------
    distance : Nansat object with distance to the coast mask in current projection

    See Also
    ---------
    `<https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/>`_
    `<http://nansat.readthedocs.io/en/latest/source/features.html#differentiating-between-land-and-water>`

    """
    # Get path to the auxilary dataset predefined in enviromental variable
    if distance_src is None:
        distance_src = os.getenv('DIST2COAST')
    # If path to the distance data source was not specified or directly provided raise an error
    if distance_src is None or not os.path.exists(distance_src):
        raise IOError('Distance to the nearest coast product does not exist - see Nansat '
                      'documentation to get it (the path is % s)' % distance_src)
    distance = Nansat(distance_src)
    # Reproject the source file on the domain of interest
    distance.reproject(dst_domain, addmask=False)
    return distance


def deprecated(func):
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn('The function <{}> was moved to the nansat.utils module and soon will not '
                      'be available from the nansat.tools'.format(func.__name__),
                      stacklevel=2)
        return func(*args, **kwargs)
    return new_func


@deprecated
def remove_keys(dict, keys):
    return utils.remove_keys(dict, keys)


@deprecated
def register_colormaps():
    utils.register_colormaps()


@deprecated
def initial_bearing(lon1, lat1, lon2, lat2):
    return utils.initial_bearing(lon1, lat1, lon2, lat2)


@deprecated
def haversine(lon1, lat1, lon2, lat2):
    return utils.haversine(lon1, lat1, lon2, lat2)


@deprecated
def add_logger(logName='', logLevel=None):
    return utils.add_logger(logName=logName, logLevel=logLevel)


@deprecated
def get_random_color(c0=None, minDist=100, low=0, high=255):
    return utils.get_random_color(c0=c0, minDist=minDist, low=low, high=high)


@deprecated
def parse_time(time_string):
    return utils.parse_time(time_string)


@deprecated
def write_domain_map(border, out_filename, lon_vec=None, lat_vec=None, lon_border=10.,
                     lat_border=10., figure_size=(6, 6), dpi=50, projection='cyl', resolution='c',
                     continets_color='coral',meridians=10, parallels=10, p_color='r', p_line='k',
                     p_alpha=0.5, padding=0., mer_labels=[False, False, False, False],
                     par_labels=[False, False, False, False], pltshow=False, labels=None):
    utils.write_domain_map(
        border, out_filename, lon_vec=lon_vec, lat_vec=lat_vec, lon_border=lon_border,
        lat_border=lat_border, figure_size=figure_size, dpi=dpi, projection=projection,
        resolution=resolution, continets_color=continets_color, meridians=meridians,
        parallels=parallels, p_color=p_color, p_line=p_line, p_alpha=p_alpha, padding=padding,
        mer_labels=mer_labels, par_labels=par_labels, pltshow=pltshow, labels=labels)
