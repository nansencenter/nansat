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
from nansat.domain import Domain
from nansat.nansat import Nansat
from nansat.tools import distance2coast, register_colormaps
from mock import patch
import numpy as np


class ToolsTest(unittest.TestCase):
    def setUp(self):
        self.d = Domain(4326, "-te 25 70 35 72 -ts 100 100")
        # define a test Nansat object
        test_domain = Domain(4326, "-lle -180 -90 180 90 -ts 500 500")
        self.n = Nansat.from_domain(test_domain, array=np.ones([500,500]))

    @patch.dict(os.environ,{'DIST2COAST':'/path/dos/not/exist'})
    def test_distance2coast_source_not_exists_envvar(self):
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

    @patch('nansat.tools.Nansat')
    def test_distance2coast_integration(self, MockNansat):
        MockNansat.return_value = self.n
        result = distance2coast(self.d)
        self.assertEqual(type(result), Nansat)

    def test_warning(self):
        register_colormaps()
        with self.assertWarns(UserWarning) as w:
            register_colormaps()
