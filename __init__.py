# Name: __init__.py
# Purpose: Use this folder as a package
#
# Author:      Asuka Yamakawa
#
# Created:     08.01.2013
# Copyright:   (c) NERSC 2013
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

