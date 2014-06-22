# Name: __init__.py
# Purpose: Use the current folder as a package
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2014
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

try:
    from ._pixfun import registerPixelFunctions
    registerPixelFunctions()
except Exception as e:
    print repr(e)
    warnings.warn('''Cannot register C pixel functions! Nansat will have reduced functionality''')

try:
    from .nsr import NSR
except ImportError:
    warnings.warn('''Cannot import NSR! Nansat will not work''')

try:
    from .domain import Domain
except ImportError:
    warnings.warn('''Cannot import Domain! Nansat will not work''')

try:
    from .nansat import Nansat
except ImportError:
    warnings.warn('''Cannot import VRT! Nansat will not work''')

try:
    from .figure import Figure
except ImportError:
    warnings.warn('''Cannot import Figure! Nansat will not work''')

try:
    from .nansatmap import Nansatmap
except ImportError:
    warnings.warn('''Cannot import Nansatmap! Nansat will not work''')

try:
    from .nansatshape import Nansatshape
except ImportError:
    warnings.warn('''Cannot import NansatOGR! Nansat will not work''')

try:
    from .nansat_tools import np, plt, Basemap, osr, ogr, gdal
except ImportError:
    warnings.warn('''Cannot import Numpy, Matplotlib! Nansat will not work''')

try:
    from .mosaic import Mosaic
except ImportError:
    warnings.warn('''Cannot import Mosaic! Mosaic will not work''')

os.environ['LOG_LEVEL'] = '30'

__all__ = ['NSR', 'Nansat',  'Nansatshape', 'Domain', 'Figure', 'Nansatmap',
           'Mosaic', 'np', 'plt', 'Basemap', 'gdal', 'ogr', 'osr']
