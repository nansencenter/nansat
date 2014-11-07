#------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified:     Anton Korosov
#
# Created:      18.06.2014
# Last modified:29.09.2014
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

from nansat import Nansat, Domain
from nansat.tools import gdal

import nansat_test_data as ntd

IS_CONDA = 'conda' in os.environ['PATH']


class NansatTest(unittest.TestCase):
    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_stere = os.path.join(ntd.test_data_path, 'stere.tif')
        plt.switch_backend('Agg')

        if not os.path.exists(self.test_file_gcps):
            raise ValueError('No test data available')

    def test_init_filename(self):
        n = Nansat(self.test_file_gcps, logLevel=40)

        self.assertEqual(type(n), Nansat)

    def test_init_domain(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        n = Nansat(domain=d, logLevel=40)

        self.assertEqual(type(n), Nansat)
        self.assertEqual(n.shape(), (500, 500))

    def test_init_domain_array(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        n = Nansat(domain=d,
                   array=np.random.randn(500, 500),
                   parameters={'name': 'band1'},
                   logLevel=40)

        self.assertEqual(type(n), Nansat)
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertEqual(n.get_metadata('name', 1), 'band1')
        self.assertEqual(n[1].shape, (500, 500))

    def test_add_band(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)
        n = Nansat(domain=d, logLevel=40)
        n.add_band(arr, {'name': 'band1'})

        self.assertEqual(type(n), Nansat)
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertEqual(n.get_metadata('name', 1), 'band1')
        self.assertEqual(n[1].shape, (500, 500))

    def test_add_bands(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)

        n = Nansat(domain=d, logLevel=40)
        n.add_bands([arr, arr],
                    [{'name': 'band1'}, {'name': 'band2'}])

        self.assertEqual(type(n), Nansat)
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertEqual(type(n[2]), np.ndarray)
        self.assertEqual(n.get_metadata('name', 1), 'band1')
        self.assertEqual(n.get_metadata('name', 2), 'band2')

    def test_bands(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        bands = n.bands()

        self.assertEqual(type(bands), dict)
        self.assertTrue(1 in bands)
        self.assertTrue('name' in bands[1])
        self.assertEqual(bands[1]['name'], 'L_645')

    def test_has_band(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        hb = n.has_band('L_645')

        self.assertTrue(hb)

    def test_export(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansat_export.nc')
        n.export(tmpfilename, driver='GTiff')

        self.assertTrue(os.path.exists(tmpfilename))

    def test_export2thredds_stere(self):
        # skip the test if anaconda is used
        if IS_CONDA:
            return
        n = Nansat(self.test_file_stere, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds.nc')
        n.export(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_dont_export2thredds_gcps(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds.nc')

        self.assertRaises(RuntimeError, n.export2thredds, tmpfilename)

    def test_resize_by_pixelsize(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        n.resize(pixelsize=500, eResampleAlg=1)

        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_factor(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        n.resize(0.5, eResampleAlg=1)

        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_width(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        n.resize(width=100, eResampleAlg=1)

        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_height(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        n.resize(height=500, eResampleAlg=1)

        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_resize(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        n.resize(0.1)
        n.resize(10)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_resize_resize.png')
        n.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(type(n[1]), np.ndarray)

    def test_get_GDALRasterBand(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        b = n.get_GDALRasterBand(1)
        arr = b.ReadAsArray()

        self.assertEqual(type(b), gdal.Band)
        self.assertEqual(type(arr), np.ndarray)

    def test_list_bands_false(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        lb = n.list_bands(False)

        self.assertEqual(type(lb), str)

    def test_reproject_domain(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        d = Domain(4326, "-te 27 70 30 72 -ts 500 500")
        n.reproject(d)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_domain.png')
        n.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n.shape(), (500, 500))
        self.assertEqual(type(n[1]), np.ndarray)

    def test_reproject_stere(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        n2 = Nansat(self.test_file_stere, logLevel=40)
        n1.reproject(n2)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_stere.png')
        n1.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n1.shape(), n2.shape())
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_reproject_gcps(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        n2 = Nansat(self.test_file_gcps, logLevel=40)
        n1.reproject(n2)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_gcps.png')
        n1.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n1.shape(), n2.shape())
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_reproject_gcps_resize(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        n2 = Nansat(self.test_file_gcps, logLevel=40)
        n1.reproject(n2)
        n1.resize(2)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_gcps_resize.png')
        n1.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n1.shape()[0], n2.shape()[0]*2)
        self.assertEqual(n1.shape()[1], n2.shape()[1]*2)
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_undo(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        shape1 = n1.shape()
        n1.resize(10)
        n1.undo()
        shape2 = n1.shape()

        self.assertEqual(shape1, shape2)

    def test_write_figure(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_figure.png')
        n1.write_figure(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_figure_band(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_figure_band.png')
        n1.write_figure(tmpfilename, 2)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_figure_clim(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_figure_clim.png')
        n1.write_figure(tmpfilename, 3, clim='hist')

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_figure_clim(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_figure_legend.png')
        n1.write_figure(tmpfilename, 3, clim='hist', legend=True)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_geotiffimage(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_geotiffimage.tif')
        n1.write_geotiffimage(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_get_time(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        t = n1.get_time()

        self.assertEqual(len(t), len(n1.bands()))
        self.assertEqual(type(t[0]), datetime.datetime)

    def test_get_metadata(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        m = n1.get_metadata()

        self.assertEqual(type(m), dict)
        self.assertTrue('fileName' in m)

    def test_get_metadata_key(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        m = n1.get_metadata('fileName')

        self.assertEqual(type(m), str)

    def test_get_metadata_wrong_key(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        m = n1.get_metadata('some_crap')

        self.assertTrue(m is None)

    def test_get_metadata_bandid(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        m = n1.get_metadata(bandID=1)

        self.assertEqual(type(m), dict)
        self.assertTrue('name' in m)

    def test_set_metadata(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        n1.set_metadata('newKey', 'newVal')
        m = n1.get_metadata('newKey')

        self.assertEqual(m, 'newVal')

    def test_set_metadata_bandid(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        n1.set_metadata('newKey', 'newVal', 1)
        m = n1.get_metadata('newKey', 1)

        self.assertEqual(m, 'newVal')

    def test_export_band(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_geotiffimage.tif')
        n1.export_band(tmpfilename, driver='GTiff')

        self.assertTrue(os.path.exists(tmpfilename))

    def test_get_transect(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        v, xy, pl = n1.get_transect([[(28.31299128, 70.93709219),
                                      (28.93691525, 70.69646524)]])
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_get_transect.png')
        plt.plot(v['1:L_645']['shape0'], xy['shape0']['latitude'])
        plt.savefig(tmpfilename)
        plt.close('all')

        self.assertTrue(len(v['1:L_645']['shape0']) > 50)
        self.assertEqual(len(v['1:L_645']['shape0']),
                         len(xy['shape0']['latitude']))
        self.assertEqual(len(v['1:L_645']['shape0']),
                         len(pl['shape0'][0]))
        self.assertEqual(type(xy['shape0']['latitude']), np.ndarray)
        self.assertEqual(type(pl['shape0'][0]), np.ndarray)

    def test_get_transect_outside(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        v, xy, pl = n1.get_transect([[(28.31299128, 70.93709219),
                                      (0.0, 0.0)]])

        self.assertTrue(len(v['1:L_645']['shape0']) > 50)
        self.assertEqual(len(v['1:L_645']['shape0']),
                         len(xy['shape0']['latitude']))
        self.assertEqual(len(v['1:L_645']['shape0']),
                         len(pl['shape0'][0]))
        self.assertEqual(type(xy['shape0']['latitude']), np.ndarray)
        self.assertEqual(type(pl['shape0'][0]), np.ndarray)

    def test_get_transect_false(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        v, xy, pl = n1.get_transect([(28.31299128, 70.93709219),
                                     (28.93691525, 70.69646524)])

        self.assertEqual(len(v['1:L_645']), 2)
        self.assertEqual(len(v['1:L_645']), len(xy))
        self.assertEqual(len(v['1:L_645']), len(pl))
        self.assertEqual(type(xy['shape0']['latitude']), np.ndarray)
        self.assertEqual(type(pl['shape0'][0]), np.ndarray)

    def test_crop(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        st, ext = n1.crop(10, 20, 50, 60)

        self.assertEqual(st, 0)
        self.assertEqual(n1.shape(), (60, 50))
        self.assertEqual(ext, (10, 20, 50, 60))
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_crop_lonlat_lims(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        st, ext = n1.crop(lonlim=[28, 29], latlim=[70.5, 71])

        self.assertEqual(st, 0)
        self.assertEqual(n1.shape(), (111, 110))
        self.assertEqual(ext, (31, 89, 110, 111))
        self.assertEqual(type(n1[1]), np.ndarray)

if __name__ == "__main__":
    unittest.main()

