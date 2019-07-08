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
    import matplotlib.pyplot as plt
    from matplotlib.colors import hex2color
except ImportError:
    MATPLOTLIB_IS_INSTALLED = False
else:
    MATPLOTLIB_IS_INSTALLED = True

try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    BASEMAP_LIB_IS_INSTALLED = False
else:
    BASEMAP_LIB_IS_INSTALLED = True

from nansat.figure import Image
from nansat.domain import Domain
from nansat.tools import get_random_color, parse_time, write_domain_map
from nansat.tests import nansat_test_data as ntd

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

    @patch('nansat.tools.MATPLOTLIB_IS_INSTALLED', False)
    def test_get_random_color__matplotlib_missing(self): 
        with self.assertRaises(ImportError):
            c0 = get_random_color()

    def test_parse_time(self):
        dt = parse_time('2016-01-19')

        self.assertEqual(type(dt), datetime.datetime)

    def test_parse_time_incorrect(self):
        dt = parse_time('2016-01-19Z')

        self.assertEqual(type(dt), datetime.datetime)

    @unittest.skipUnless(BASEMAP_LIB_IS_INSTALLED, 'Basemap is required')
    def test_write_domain_map(self):
        plt.switch_backend('Agg')
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        border = d.get_border()
        tmpfilename = os.path.join(ntd.tmp_data_path, 'domain_write_map.png')
        write_domain_map(border, tmpfilename, labels=['Patch1'])
        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (50, 50))

    @patch('nansat.tools.BASEMAP_LIB_IS_INSTALLED', False)
    def test_write_domain_map__basemap_missing(self):
        plt.switch_backend('Agg')
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        border = d.get_border()
        tmpfilename = os.path.join(ntd.tmp_data_path, 'domain_write_map.png')
        with self.assertRaises(ImportError):
            write_domain_map(border, tmpfilename, labels=['Patch1'])

    @unittest.skipUnless(BASEMAP_LIB_IS_INSTALLED, 'Basemap is required')
    def test_write_domain_map_dpi100(self):
        plt.switch_backend('Agg')
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        border = d.get_border()
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'domain_write_map_dpi100.png')
        write_domain_map(border, tmpfilename, dpi=100)
        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (100, 100))

    @unittest.skipUnless(BASEMAP_LIB_IS_INSTALLED, 'Basemap is required')
    def test_write_domain_map_labels(self):
        plt.switch_backend('Agg')
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        border = d.get_border()
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'domain_write_map_labels.png')
        write_domain_map(border, tmpfilename,
                    mer_labels=[False, False, False, True],
                    par_labels=[True, False, False, False])
        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
