# Name:    nansat_tools.py
# Purpose: General functions/tools used within the Nansat toolbox
#
# Authors:      Asuka Yamakava, Anton Korosov, Knut-Frode Dagestad
#
# Created:     18.02.2012
# Copyright:   (c) NERSC 2012
# Licence:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details:
# http://www.gnu.org/licenses/

from math import atan2, sin, cos, radians, degrees
from scipy import mod
import logging

def initial_bearing(lon1, lat1, lon2, lat2):
        '''Initial bearing when traversing from point1 (lon1, lat1)
        to point2 (lon2, lat2)
        
        See http://www.movable-type.co.uk/scripts/latlong.html
        
        Parameters
        ----------
        lon1, lat1: float
            longitude and latitude of start point
        lon2, lat2: float
            longitude and latitude of end point

        Returns
        -------
        initial_bearing: float
            The initial bearing (azimuth direction) when heading out
            from the start point towards the end point along a 
            great circle.'''
        rlon1 = radians(lon1)
        rlat1 = radians(lat1)
        rlon2 = radians(lon2)
        rlat2 = radians(lat2)
        deltalon = rlon2-rlon1
        bearing = atan2(sin(rlon2-rlon1)*cos(rlat2),
                        cos(rlat1)*sin(rlat2) -
                        sin(rlat1)*cos(rlat2)*cos(rlon2-rlon1))
        return mod(degrees(bearing)+360, 360)

def add_logger(logName='', logLevel=30):
    ''' Creates and returns logger with default formatting for Nansat
    
    Parameters:
    -----------
        logName: string, optional
            Name of the logger
        logLevel: int, option, default=30
            level at which logging is performed
    
    Returns:
    --------
        logging.logger
        See also: http://docs.python.org/howto/logging.html
    '''

    # create (or take already existing) logger
    # with default logging level WARNING
    logger = logging.getLogger(logName)
    logger.setLevel(logLevel)
    
    # if logger already exits, default stream handler was already added
    # otherwise create and add a new handler
    if len(logger.handlers) == 0:
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        # create formatter
        formatter = logging.Formatter('%(asctime)s|%(name)s|%(funcName)s|%(levelname)s|%(message)s', datefmt='%I:%M:%S')
        # add formatter to ch
        ch.setFormatter(formatter)
        # add ch to logger
        logger.addHandler(ch)

    logger.handlers[0].setLevel(logLevel)    
    
    return logger

