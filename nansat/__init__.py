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
import os, sys
import warnings

from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)

# check if pixel functions were compiled using setup_tools
try:
    from ._pixfun import registerPixelFunctions
    registerPixelFunctions()
except Exception as e:
    print repr(e)
    warnings.warn('''Cannot register C pixel functions!
                     Either nansat was not installed using setup.py or
                     pixel functions were not compiled automatically.
                     For development, use "python setup.py build_ext --inplace"
                     to compile pixel functions manually into the source tree.
                     ''')

from .nsr import NSR
from .domain import Domain
from .nansat import Nansat

__all__ = ['NSR', 'Domain', 'Nansat']

try:
    from .figure import Figure
except ImportError:
    warnings.warn('''Cannot import Figure! Nansat will not make figures!''')
else:
    __all__.append('Figure')

try:
    from .nansatmap import Nansatmap
except ImportError:
    warnings.warn('''Cannot import Nansatmap! Nansat will not make maps!''')
else:
    __all__.append('Nansatmap')

try:
    from .mosaic import Mosaic
except ImportError:
    warnings.warn('''Cannot import Mosaic! Nansat will not mosaic files!''')
else:
    __all__.append('Mosaic')

os.environ['LOG_LEVEL'] = '30'

# import some libraries for convenience
from .tools import gdal, ogr
import numpy as np
import matplotlib.pyplot as plt
__all__ += ['gdal', 'ogr', 'np', 'plt']
