#------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Anton Korosov, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:03.03.2015 09:36
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
from scipy.io.netcdf import netcdf_file

from nansat import Nansat, Domain
from nansat.tools import gdal, OptionError

import nansat_test_data as ntd

IS_CONDA = 'conda' in os.environ['PATH']


class NansatTest(unittest.TestCase):
    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_stere = os.path.join(ntd.test_data_path, 'stere.tif')
        self.test_file_complex = os.path.join(ntd.test_data_path, 'complex.nc')
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

    def test_export_netcdf(self):
        ''' Test export and following import of data with bands containing
        np.nan values
        '''
        n = Nansat(self.test_file_gcps)
        arrNoNaN = np.random.randn(n.shape()[0], n.shape()[1])
        n.add_band(arrNoNaN, {'name': 'testBandNoNaN'})
        arrWithNaN = arrNoNaN.copy()
        arrWithNaN[n.shape()[0]/2-10:n.shape()[0]/2+10,
            n.shape()[1]/2-10:n.shape()[1]/2+10] = np.nan
        n.add_band(arrWithNaN, {'name': 'testBandWithNaN'})
        n.export('test.nc')
        exported = Nansat('test.nc')
        earrNoNaN = exported['testBandNoNaN']
        # Use allclose to allow some roundoff errors
        self.assertTrue(np.allclose(arrNoNaN, earrNoNaN))
        earrWithNaN = exported['testBandWithNaN']
        np.testing.assert_allclose(arrWithNaN, earrWithNaN)
        os.unlink('test.nc')

    def test_add_band_twice(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)
        n = Nansat(domain=d, logLevel=40)
        n.add_band(arr, {'name': 'band1'})
        n.add_band(arr, {'name': 'band2'})

        self.assertEqual(type(n), Nansat)
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertEqual(type(n[2]), np.ndarray)
        self.assertEqual(n.get_metadata('name', 1), 'band1')
        self.assertEqual(n.get_metadata('name', 2), 'band2')
        self.assertEqual(n[1].shape, (500, 500))
        self.assertEqual(n[2].shape, (500, 500))

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

    def test_add_subvrts_only_to_one_nansat(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)

        n1 = Nansat(domain=d, logLevel=40)
        n2 = Nansat(domain=d, logLevel=40)
        n1.add_band(arr, {'name': 'band1'})

        self.assertEqual(type(n1.vrt.bandVRTs), dict)
        self.assertTrue(len(n1.vrt.bandVRTs) > 0)
        self.assertEqual(n2.vrt.bandVRTs, {})

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
        n.export(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_export_gtiff(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansat_export.tif')
        n.export(tmpfilename, driver='GTiff')

        self.assertTrue(os.path.exists(tmpfilename))

    def test_export_band(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export_band.tif')
        n.export(tmpfilename, bands= [1], driver='GTiff')
        n = Nansat(tmpfilename, mapperName='generic')

        self.assertTrue(os.path.exists(tmpfilename))
        self.assertEqual(n.vrt.dataset.RasterCount, 1)

    def test_reproject_and_export_band(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        n2 = Nansat(self.test_file_stere, logLevel=40)
        n1.reproject(n2)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_export_band.nc')
        n1.export(tmpfilename, bands=[1])

        n = Nansat(tmpfilename, mapperName='generic')
        self.assertTrue(os.path.exists(tmpfilename))
        self.assertEqual(n.vrt.dataset.RasterCount, 1)

    def test_export2thredds_stere_one_band(self):
        # skip the test if anaconda is used
        if IS_CONDA:
            return
        n = Nansat(self.test_file_stere, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds_1b.nc')
        n.export2thredds(tmpfilename, ['L_469'])

        self.assertTrue(os.path.exists(tmpfilename))


    def test_export2thredds_stere_many_bands(self):
        # skip the test if anaconda is used
        if IS_CONDA:
            return
        n = Nansat(self.test_file_stere, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds_3b.nc')
        bands = {
            'L_645' : {'type': '>i1'},
            'L_555' : {'type': '>i1'},
            'L_469' : {'type': '>i1'},
        }
        n.export2thredds(tmpfilename, bands)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_dont_export2thredds_gcps(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds.nc')
        self.assertRaises(OptionError, n.export2thredds, tmpfilename, ['L_645'])

    def test_export2thredds_longlat(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   "-te 27 70 31 72 -ts 200 200")
        n.reproject(d)

        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds_longlat.nc')
        n.export2thredds(tmpfilename, ['L_469'])
        ncI = netcdf_file(tmpfilename, 'r')
        ncIVar = ncI.variables['L_469']
        self.assertTrue(ncIVar.grid_mapping in ncI.variables.keys())

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

    def test_resize_complex_algAverage(self):
        n = Nansat(self.test_file_complex, logLevel=40)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            n.resize(0.5, eResampleAlg=-1)

            self.assertTrue(len(w)==1)
            self.assertTrue(issubclass(w[-1].category, UserWarning))
            self.assertTrue(
                        'The imaginary parts of complex numbers ' \
                        'are lost when resampling by averaging '
                            in str(w[-1].message))

    def test_resize_complex_alg0(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=0)

        self.assertTrue(np.any(n[1].imag!=0))

    def test_resize_complex_alg1(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=1)

        self.assertTrue(np.any(n[1].imag!=0))

    def test_resize_complex_alg2(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=2)

        self.assertTrue(np.any(n[1].imag!=0))

    def test_resize_complex_alg3(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=3)

        self.assertTrue(np.any(n[1].imag!=0))

    def test_resize_complex_alg4(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=4)

        self.assertTrue(np.any(n[1].imag!=0))

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

    def test_get_no_transect_interactive(self):
        ''' Check that get_transect does returns None if interactive mode is
        '''
        import matplotlib.pyplot as plt
        plt.ion()
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        noneResult = n1.get_transect()

        self.assertEqual(noneResult, None)
        plt.ioff()

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

    def test_watermask(self):
        ''' if watermask data exists: should fetch array with watermask
            else:                     should raise an error'''
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        mod44path = os.getenv('MOD44WPATH')
        if mod44path is not None and os.path.exists(mod44path+ '/MOD44W.vrt'):
            wm = n1.watermask()[1]
            self.assertEqual(type(wm), np.ndarray)
            self.assertEqual(wm.shape[0], n1.shape()[0])
            self.assertEqual(wm.shape[1], n1.shape()[1])

    def test_watermask_fail(self):
        ''' Nansat.watermask should raise an IOError'''
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        os.environ['MOD44WPATH'] = '/fakepath'
        self.assertRaises(IOError, n1.watermask)

if __name__ == "__main__":
    unittest.main()
