#------------------------------------------------------------------------------
# Name:         test_tools.py
# Purpose:      Test tools from nansat.tools
# Author:       Artem Moiseev
# Created:      17.01.2020
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------

import os
import unittest
from mock import patch, DEFAULT

import numpy as np
try:
    import matplotlib.pyplot as plt
    plt.switch_backend('Agg')
except ImportError:
    MATPLOTLIB_IS_INSTALLED = False
else:
    MATPLOTLIB_IS_INSTALLED = True
try:
    import cartopy
except ImportError:
    CARTOPY_IS_INSTALLED = False
else:
    CARTOPY_IS_INSTALLED = True

import nansat
from nansat.domain import Domain
from nansat.figure import Image
from nansat.nansat import Nansat as NANSAT
from nansat.tests import nansat_test_data as ntd
from nansat.tools import (distance2coast,
                            register_colormaps,
                            get_domain_map,
                            show_domain_map,
                            save_domain_map,
                            get_domain_map)


class ToolsTest(unittest.TestCase):
    def setUp(self):
        self.d = Domain(4326, "-te 25 70 35 72 -ts 100 100")
        # define a test Nansat object
        test_domain = Domain(4326, "-lle -180 -90 180 90 -ts 500 500")
        self.n = NANSAT.from_domain(test_domain, array=np.ones([500,500]))

    @patch('nansat.tools.os.getenv')
    def test_distance2coast_source_not_exists_envvar(self, mock_getenv):
        mock_getenv.return_value='/path/dos/not/exist'
        with self.assertRaises(IOError) as err:
            distance2coast(self.d)
        self.assertEqual('Distance to the nearest coast product does not exist - '
                         'see Nansat documentation to get it (the path is '
                         '/path/dos/not/exist)', str(err.exception))

    def test_distance2coast_source_not_exists_attribute(self):
        with self.assertRaises(IOError) as err:
            distance2coast(self.d, distance_src='/path/dos/not/exist')
        self.assertEqual('Distance to the nearest coast product does not exist - '
                         'see Nansat documentation to get it (the path is '
                         '/path/dos/not/exist)', str(err.exception))

    @patch.multiple('nansat.tools', Nansat=DEFAULT, os=DEFAULT)
    def test_distance2coast_integration(self, Nansat, os):
        Nansat.return_value = self.n
        os.path.exists.return_value=True
        result = distance2coast(self.d)
        self.assertEqual(type(result), NANSAT)

    def test_warning(self):
        register_colormaps()
        with self.assertWarns(UserWarning) as w:
            register_colormaps()

    @unittest.skipUnless(CARTOPY_IS_INSTALLED, 'Cartopy is required')
    def test_get_domain_map(self):
        ax = get_domain_map(self.d)
        self.assertIsInstance(ax, plt.Axes)

    def test_get_domain_map_no_cartopy(self):
        nansat.tools.CARTOPY_IS_INSTALLED = False
        with self.assertRaises(ImportError) as err:
            ax = get_domain_map(self.d)
        self.assertIn('Cartopy is not installed', str(err.exception))
        nansat.tools.CARTOPY_IS_INSTALLED = CARTOPY_IS_INSTALLED

    @unittest.skipUnless(CARTOPY_IS_INSTALLED, 'Cartopy is required')
    def test_save_domain_map(self):
        tmpfilename = os.path.join(ntd.tmp_data_path, 'domain_save_map.png')
        save_domain_map(self.d, tmpfilename)
        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (100, 100))
