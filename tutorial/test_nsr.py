#!/usr/bin/env python
# Name:    test_nsr.py
# Purpose: Tutorial for NSR class
# Authors:      Anton Korosov, Asuka Yamakawa, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      13.02.2014
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

from nansat import NSR

print 'NSR(None)', NSR(None)
print 'NSR()', NSR()
print 'NSR(+proj=longlat)', NSR('+proj=longlat')
print 'NSR(crap)', NSR('crap')
print 'NSR(osr.SRS_WKT_WGS84)', NSR(osr.SRS_WKT_WGS84)
print 'NSR(4326)', NSR(4326)
print 'NSR(NSR(4326))', NSR(NSR(4326))
print 'NSR(+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=75 +lon_0=10 +no_defs)', NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=75 +lon_0=10 +no_defs')

print '\n*** nsr_test completed successfully. Output files are found here:' + oFileName
