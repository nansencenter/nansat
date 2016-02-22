# ------------------------------------------------------------------------------
# Name:         test_wkv.py
# Purpose:      Test the Yaml file with Well Knwon Variables
#
# Author:       Anton Korosov
#
# Created:      22.02.2016
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------------------
import unittest
import os
import yaml

import nansat



class WKVTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(os.path.dirname(nansat.__file__),
                  'wkv.yml')) as f:
            self.wkvList = yaml.load(f)

    def test_exists_standard_name(self):
        for wkvNode in self.wkvList:
            self.assertTrue('standard_name' in wkvNode)

    def test_exists_short_name(self):
        for wkvNode in self.wkvList:
            self.assertTrue('short_name' in wkvNode)

    def test_exists_long_name(self):
        for wkvNode in self.wkvList:
            self.assertTrue('long_name' in wkvNode)

    def test_exists_units(self):
        for wkvNode in self.wkvList:
            self.assertTrue('units' in wkvNode)

    def test_exists_colormap(self):
        for wkvNode in self.wkvList:
            self.assertTrue('colormap' in wkvNode)

if __name__ == "__main__":
    unittest.main()
