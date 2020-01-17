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


def distance2coast(dst_domain, distance_src=None):
    """Method is estmating distance to the nearest coast (in km) for each pixcel in the 
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
    `<http://nansat.readthedocs.io/en/latest/source/features.html#differentiating-between-land-and-water>`_

    """
    # Get path to theauxilary dataset predefined in enviromental variable 
    if not distance_src:
        distance_src = os.getenv('DIST2COAST')
    try:
        distance = Nansat(distance_src)
    except TypeError:
        # Raise an error if the enviromental variable does not exist
        raise IOError('Distance to the nearest coast product is not specified. Add the path to the'
                      'directory with this data to an environment variable named DIST2COAST or explisetly'
                      'provide to the function as <distance_src>')
    except OSError:
        # Raise an error if the enviromental variable was incorectly defined
        raise IOError('Distance to the nearest coast product does not exist - see Nansat'
                      'documentation to get it (the path is % s)' % distance_src)
    # Reproject the source file on the domain of interest
    distance.reproject(dst_domain, addmask=False)
    return distance
