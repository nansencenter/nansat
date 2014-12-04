#------------------------------------------------------------------------------
# Name:         test_mosaic.py
# Purpose:      Test the Nansat class
#
# Author:       Anton Korosov
# Modified:     Anton Korosov
#
# Created:      04.12.2014
# Last modified:04.12.2014 09:00
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------
import unittest
import warnings
import os
import sys
import glob
from types import ModuleType, FloatType
import datetime
import matplotlib.pyplot as plt
import numpy as np

from nansat import Nansat, Domain, Mosaic
from nansat.tools import gdal

import nansat_test_data as ntd

class MosaicTest(unittest.TestCase):
    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_stere = os.path.join(ntd.test_data_path, 'stere.tif')
        self.test_file_complex = os.path.join(ntd.test_data_path, 'complex.nc')
        self.domain = Domain(4326, '-lle 27 70 31 72 -ts 700 650')
        plt.switch_backend('Agg')

        if not os.path.exists(self.test_file_gcps):
            raise ValueError('No test data available')

    def test_init(self):
        mo = Mosaic(domain=self.domain)

        self.assertEqual(type(mo), Mosaic)

    def test_average(self):
        mo = Mosaic(domain=self.domain)
        mo.average([self.test_file_gcps, self.test_file_stere],
                    bands=['L_645', 'L_555', 'L_469'])

        mask = mo['mask']
        L_645 = mo['L_645']
        L_555 = mo['L_555']
        L_469 = mo['L_469']

        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'mosaic_export.nc')
        bands = {
            'L_645' : {'type': '>i1'},
            'L_555' : {'type': '>i1'},
            'L_469' : {'type': '>i1'},
        }
        mo.export2thredds(tmpfilename, bands)


if __name__ == "__main__":
    unittest.main()
