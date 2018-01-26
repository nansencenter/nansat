#------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:16.03.2015 13:19
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------
from __future__ import absolute_import, unicode_literals

import unittest
import warnings
import os
import sys
import glob
import datetime
if sys.version_info.major == 2:
    from mock import patch, PropertyMock, Mock, MagicMock, DEFAULT
else:
    from unittest.mock import patch, PropertyMock, Mock, MagicMock, DEFAULT
try:
    import Image
    import ImageDraw
    import ImageFont
except:
    from PIL import Image, ImageDraw, ImageFont

import gdal
import numpy as np

from nansat import Figure, Nansat, Domain

from nansat.tests import nansat_test_data as ntd

IS_CONDA = 'conda' in os.environ['PATH']


class FigureTest(unittest.TestCase):
    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        if not os.path.exists(self.test_file_gcps):
            raise ValueError('No test data available')

    def test_init_array(self):
        f = Figure(np.zeros((10,10)))

        self.assertEqual(type(f), Figure)


    def test_get_auto_ticks_number(self):
        n = Nansat(self.test_file_gcps)
        lon, lat = n.get_geolocation_grids()
        f = Figure(lon)
        lonTicks = f._get_auto_ticks(5, lon)
        latTicks = f._get_auto_ticks(5, lat)

        self.assertEqual(len(lonTicks), 5)
        n.logger.error(str(lonTicks))
        n.logger.error(str(latTicks))

    def test_get_auto_ticks_vector(self):
        n = Nansat(self.test_file_gcps)
        lon, lat = n.get_geolocation_grids()
        f = Figure(lon)
        lonTicks = f._get_auto_ticks([28, 29, 30, 100], lon)

        self.assertEqual(len(lonTicks), 3)

    def test_add_latlon_grids_auto(self):
        ''' Should create figure with lon/lat gridlines spaced automatically '''
        tmpfilename = os.path.join(ntd.tmp_data_path, 'figure_latlon_grids_auto.png')
        n = Nansat(self.test_file_gcps)
        b = n[1]
        lon, lat = n.get_geolocation_grids()

        f = Figure(b)
        f.process(clim='hist', lonGrid=lon, latGrid=lat)
        f.save(tmpfilename)

        self.assertEqual(type(f), Figure)
        self.assertTrue(os.path.exists(tmpfilename))

    def test_add_latlon_grids_number(self):
        ''' Should create figure with lon/lat gridlines given manually '''
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'figure_latlon_grids_number.png')
        n = Nansat(self.test_file_gcps)
        n.resize(3)
        b = n[1]
        lon, lat = n.get_geolocation_grids()

        f = Figure(b)
        f.process(cmax=100, lonGrid=lon,
                               latGrid=lat,
                               lonTicks=7,
                               latTicks=7)
        f.save(tmpfilename)

        self.assertEqual(type(f), Figure)
        self.assertTrue(os.path.exists(tmpfilename))

    def test_add_latlon_grids_list(self):
        ''' Should create figure with lon/lat gridlines given manually '''
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'figure_latlon_grids_list.png')
        n = Nansat(self.test_file_gcps)
        b = n[1]
        lon, lat = n.get_geolocation_grids()

        f = Figure(b)
        f.process(clim='hist', lonGrid=lon,
                               latGrid=lat,
                               lonTicks=[28, 29, 30],
                               latTicks=[70.5, 71, 71.5, 73])
        f.save(tmpfilename)

        self.assertEqual(type(f), Figure)
        self.assertTrue(os.path.exists(tmpfilename))


    def test_get_tick_index_from_grid(self):
        ''' Should return indeces of pixel closest to ticks '''
        n = Nansat(self.test_file_gcps)
        lon, lat = n.get_geolocation_grids()

        f = Figure(lon)
        lonTicksIdx = f._get_tick_index_from_grid([28.5, 29], lon, 1, lon.shape[1])
        latTicksIdx = f._get_tick_index_from_grid([71, 71.5], lat, lat.shape[0], 1)
        n.logger.error(str(lonTicksIdx))
        n.logger.error(str(latTicksIdx))

    @patch.object(Figure, '__init__', return_value=None)
    def test_apply_logarithm(self, mock1):
        f = Figure()
        f.array = np.ones((1,2,2))*0.1
        f.apply_logarithm()
        self.assertTrue(np.allclose(np.ones((1,2,2))*0.31622777, f.array))

    @patch.object(Figure, '__init__', return_value=None)
    def test_make_transparent_color(self, mock1):
        f = Figure()

        img_array = np.array([[1.,2.],[3.,4.]])
        f.pilImg = Image.fromarray(img_array)
        f.reprojMask = np.zeros((2,2)) == 1
        f.transparency = [1,1,1]
        f._make_transparent_color()
        self.assertTrue((np.array(f.pilImg)[:,:,3] == np.array([[0,255],[255,255]])).all())


if __name__ == "__main__":
    unittest.main()
