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
from nansat.tools import distance2coast
from mock import patch


class ToolsTest(unittest.TestCase):
    def setUp(self):
        self.d = Domain(4326, "-te 25 70 35 72 -ts 500 500")

    @patch.dict(os.environ,{'DIST2COAST':'/path/dos/not/exist'})
    def test_distance2coast_source_not_exists(self):
        with self.assertRaises(IOError) as err:
            distance2coast(self.d)
        self.assertEqual('Distance to the nearest coast product does not exist - '
                         'see Nansat documentation to get it (the path is '
                         '/path/dos/not/exist)', str(err.exception))
