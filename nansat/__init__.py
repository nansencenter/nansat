# Name: __init__.py
# Purpose: Use the current folder as a package
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov, Aleksander Vines
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2015
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
import sys
import warnings
import importlib
pixfun_module_name = 'nansat._pixfun_py{0}'.format(sys.version_info[0])

# check if pixel functions were compiled using setup_tools
try:
    pixfun = importlib.import_module(pixfun_module_name)
    pixfun.registerPixelFunctions()
except ImportError as e:
    print(e)
    warnings.warn('''Cannot register C pixel functions!
                     Either nansat was not installed using setup.py or
                     pixel functions were not compiled automatically.
                     For development, use "python setup.py build_ext --inplace"
                     to compile pixel functions manually into the source tree.
                     ''')
from nansat.nsr import NSR
from nansat.domain import Domain
from nansat.nansat import Nansat
from nansat.figure import Figure

__all__ = ['NSR', 'Domain', 'Nansat', 'Figure']

os.environ['LOG_LEVEL'] = '30'
