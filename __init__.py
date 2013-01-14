# Name: __init__.py
# Purpose: Use the current folder as a package
# Authors:      Asuka Yamakava, Anton Korosov, Knut-Frode Dagestad,
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

import warnings

try:
    from domain import Domain
except ImportError:
    warnings.warn('''Cannot import Domain! Nansat will not work''')

try:
    from nansat import Nansat
except ImportError:
    warnings.warn('''Cannot import VRT! Nansat will not work''')

try:
    from figure import Figure
except ImportError:
    warnings.warn('''Cannot import Figure! Nansat will not work''')

__all__ = ["Nansat", "Domain", "Figure"]

