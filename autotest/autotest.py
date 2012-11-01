#!/usr/bin/env python

# Name:     autotest.py
# Purpose:  Test basic functionality of Nansat, useful e.g. when 
#               upgrading Nansat or local system 
#
# Authors:   Knut-Frode Dagestad
#
# Created:     1.11.2012
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

from nansat import *

# Test generation of Domain from strings
srsString = "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs"
extentString = "-lle 1 40 65 66 -ts 1000 1000"
d = Domain(srsString, extentString)

