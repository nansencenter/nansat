# ------------------------------------------------------------------------------
# Name:         nansat_test_base.py
# Purpose:      Basic class for Nansat tests
#
# Author:       Anton Korosov
#
# Created:      20.03.2018
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------------------
from __future__ import print_function, absolute_import, division
import os
import sys
import unittest
import tempfile
import pythesint

from mock import patch, PropertyMock, Mock, MagicMock, DEFAULT

from nansat.tests import nansat_test_data as ntd


class NansatTestBase(unittest.TestCase):

    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_stere = os.path.join(ntd.test_data_path, 'stere.tif')
        self.test_file_complex = os.path.join(ntd.test_data_path, 'complex.nc')
        self.test_file_arctic = os.path.join(ntd.test_data_path, 'arctic.nc')
        self.tmp_data_path = ntd.tmp_data_path
        self.default_mapper = 'generic'
        fd, self.tmp_filename = tempfile.mkstemp(suffix='.nc')
        os.close(fd)

        if not os.path.exists(self.test_file_gcps):
            raise ValueError('No test data available')
        # Mock several Pythesint functions to avoid network connection
        self.patcher = patch.multiple(pythesint, get_wkv_variable=DEFAULT,
                                                  get_gcmd_instrument=DEFAULT,
                                                  get_gcmd_platform=DEFAULT)
        self.mock_pti = self.patcher.start()
        self.mock_pti['get_gcmd_instrument'].return_value=dict(short_name='MODIS')
        self.mock_pti['get_gcmd_platform'].return_value=dict(short_name='AQUA')
        self.mock_pti['get_wkv_variable'].return_value=dict(short_name='swathmask')

    def tearDown(self):
        self.patcher.stop()
        try:
            os.unlink(self.tmp_filename)
        except OSError:
            pass

