#-------------------------------------------------------------------------------
# Name: __init__
# Purpose:
#
# Author:      asumak
#
# Created:     08.01.2013
# Copyright:   (c) asumak 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import warnings

try:
    from domain import *
except ImportError:
    warnings.warn('''Cannot import Domain! Nansat will not work''')

try:
    from nansat import *
except ImportError:
    warnings.warn('''Cannot import VRT! Nansat will not work''')

try:
    from vrt import *
except ImportError:
    warnings.warn('''Cannot import VRT! Nansat will not work''')

try:
    from figure import *
except ImportError:
    warnings.warn('''Cannot import Figure! Nansat will not work''')

__all__ = ['nansat', 'Nansat', 'Figure', 'Domain']

