#------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:16.04.2015 10:48
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------
from __future__ import unicode_literals, absolute_import
import os
import unittest
import datetime
import warnings

from mock import patch

try:
    import matplotlib
except ImportError:
    MATPLOTLIB_IS_INSTALLED = False
else:
    MATPLOTLIB_IS_INSTALLED = True
    import matplotlib.pyplot as plt
    from matplotlib.colors import hex2color

from nansat.utils import get_random_color, parse_time, register_colormaps
from nansat.tests import nansat_test_data as ntd


class UtilsTest(unittest.TestCase):
    @unittest.skipUnless(MATPLOTLIB_IS_INSTALLED, 'Matplotlib is required')
    def test_get_random_color(self):
        ''' Should return HEX code of random color '''
        c0 = get_random_color()
        c1 = get_random_color(c0)
        c2 = get_random_color(c1, 300)

        self.assertEqual(type(hex2color(c0)), tuple)
        self.assertEqual(type(hex2color(c1)), tuple)
        self.assertEqual(type(hex2color(c2)), tuple)

    @patch('nansat.utils.MATPLOTLIB_IS_INSTALLED', False)
    def test_get_random_color__matplotlib_missing(self):
        with self.assertRaises(ImportError):
            c0 = get_random_color()

    def test_parse_time(self):
        dt = parse_time('2016-01-19')

        self.assertEqual(type(dt), datetime.datetime)

    def test_parse_time_incorrect(self):
        dt = parse_time('2016-01-19Z')

        self.assertEqual(type(dt), datetime.datetime)

    @unittest.skipUnless(MATPLOTLIB_IS_INSTALLED, 'Matplotlib is required')
    def test_register_colormaps(self):
        register_colormaps()
        self.assertIn('obpg', plt.colormaps())
        self.assertIn('ak01', plt.colormaps())
