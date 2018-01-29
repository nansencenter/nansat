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

try:
    if 'DISPLAY' not in os.environ:
        import matplotlib; matplotlib.use('Agg')
    from matplotlib.colors import hex2color
except ImportError:
    MATPLOTLIB_IS_INSTALLED = False
else:
    MATPLOTLIB_IS_INSTALLED = True

from nansat.warnings import NansatFutureWarning
from nansat.tools import get_random_color, parse_time
from nansat.tools import (OptionError,
                            ProjectionError,
                            GDALError,
                            NansatReadError,
                            GeolocationError,
                            WrongMapperError)

class ToolsTest(unittest.TestCase):
    @unittest.skipUnless(MATPLOTLIB_IS_INSTALLED, 'Matplotlib is required')
    def test_get_random_color(self):
        ''' Should return HEX code of random color '''
        c0 = get_random_color()
        c1 = get_random_color(c0)
        c2 = get_random_color(c1, 300)

        self.assertEqual(type(hex2color(c0)), tuple)
        self.assertEqual(type(hex2color(c1)), tuple)
        self.assertEqual(type(hex2color(c2)), tuple)

    def test_parse_time(self):
        dt = parse_time('2016-01-19')

        self.assertEqual(type(dt), datetime.datetime)

    def test_parse_time_incorrect(self):
        dt = parse_time('2016-01-19Z')

        self.assertEqual(type(dt), datetime.datetime)

    def test_OptionError_warning(self):
        with warnings.catch_warnings(record=True) as recorded_warnings:
            with self.assertRaises(OptionError):
                raise OptionError
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)

    def test_ProjectionError_warning(self):
        with warnings.catch_warnings(record=True) as recorded_warnings:
            with self.assertRaises(ProjectionError):
                raise ProjectionError
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)

    def test_GDALError_warning(self):
        with warnings.catch_warnings(record=True) as recorded_warnings:
            with self.assertRaises(GDALError):
                raise GDALError
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)

    def test_NansatReadError_warning(self):
        with warnings.catch_warnings(record=True) as recorded_warnings:
            with self.assertRaises(NansatReadError):
                raise NansatReadError
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)

    def test_GeolocationError_warning(self):
        with warnings.catch_warnings(record=True) as recorded_warnings:
            with self.assertRaises(GeolocationError):
                raise GeolocationError
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)

    def test_WrongMapperError_warning(self):
        with warnings.catch_warnings(record=True) as recorded_warnings:
            with self.assertRaises(WrongMapperError):
                raise WrongMapperError
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)

