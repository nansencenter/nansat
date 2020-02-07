# Name:    tools.py
# Purpose: Collection of functions which use core Nansat classes for
#          scientific data analysis and visualization
# Authors:      Artem Moiseev, Anton Korosov
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

try:
    import matplotlib.pyplot as plt
except ImportError:
    MATPLOTLIB_IS_INSTALLED = False
else:
    MATPLOTLIB_IS_INSTALLED = True

try:
    import cartopy.crs as ccrs
except ImportError:
    CARTOPY_IS_INSTALLED = False
else:
    CARTOPY_IS_INSTALLED = True

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

def get_domain_map(domain,
                   crs=None,
                   lon_margin=10.,
                   lat_margin=10.,
                   lw=1,
                   linestyle='b.-',
                   fill_color='coral',
                   fill_alpha=0.5,
                   draw_gridlines=True,
                   draw_labels=True,
                   grid_lw=2,
                   grid_color='gray',
                   grid_alpha=0.5,
                   grid_linestyle='--',
                   **kwargs):
    """
    Create a pyplot figure axis with Domain map, Cartopy projection, coastlines

    Parameters
    -----------
    domain : Domain
        the desired Domain to plot
    crs : cartopy CRS or None
        projection of the map, cartopy.crs.PlateCaree by default
    lon_margin : float
        10, horisontal border around patch (degrees of longitude)
    lat_margin : float
        10, vertical border around patch (degrees of latitude)
    linestyle : str
        domain line style
    lw : float
        domain line width
    fill_color : str
        domain fill color
    fill_alpha : float
        domain fill transparency
    draw_gridlines : bool
        Add gridlines to the plot?
    draw_labels : bool
        Add labels to the plot?
    grid_lw : float
        gridlines line width
    grid_color : str,
        gridlines color
    grid_alpha : float
        gridlines transparency
    grid_linestyle : str
        gridlines style

    Returns
    -------
    ax : pyplot/cartopy axes
        Axes with the map and domain patch

    """
    if not CARTOPY_IS_INSTALLED:
        raise ImportError(' Cartopy is not installed. Cannot use get_domain_map '
                          ' Enable by: conda install -c conda-forge cartopy ')

    if crs is None:
        crs = ccrs.PlateCarree()

    blon, blat = domain.get_border()

    ax = plt.subplot(111, projection=crs)
    ax.stock_img()
    ax.coastlines()
    ax.plot(blon, blat, linestyle,
                        lw=lw,
                        transform=ccrs.PlateCarree())
    ax.fill(blon, blat, color=fill_color,
                        alpha=fill_alpha,
                        transform=crs)
    plt.xlim([blon.min()-lon_margin, blon.max()+lon_margin])
    plt.ylim([blat.min()-lat_margin, blat.max()+lat_margin])
    if draw_gridlines:
        gl = ax.gridlines(draw_labels=draw_labels,
                          linewidth=grid_lw,
                          color=grid_color,
                          alpha=grid_alpha,
                          linestyle=grid_linestyle,
                          crs=crs)
    return ax


def show_domain_map(domain, **kwargs):
    """ Show Domain map interactively

    Parameters
    ----------
    domain : Domain
        the Domain to show
    **kwargs : dict
        parameters for nansat.tools.get_domain_map

    """
    ax = get_domain_map(domain, **kwargs)
    plt.show()

def save_domain_map(domain,
                    filename,
                    figsize=(5,5),
                    dpi=150,
                    bbox_inches='tight',
                    pad_inches=0,
                    **kwargs):
    """ Save Domain to PNG file

    Parameters
    ----------
    domain : Domain
        the Domain to show
    filename : str
        destination filename
    figsize : (int, int)
        size of the figure in inches
    dpi : PNG resolution
    bbox_inches : str or int
        see pyplot.savefig
    pad_inches : int
        see pyplot.savefig
    **kwargs : dict
        parameters for nansat.tools.get_domain_map

    """
    plt.figure(figsize=figsize)
    ax = get_domain_map(domain, **kwargs)
    plt.savefig(filename, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close()

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
