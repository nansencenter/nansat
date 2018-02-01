# ------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa, Anton Korosov
#
# Created:      18.06.2014
# Last modified:24.08.2017 14:00
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------------------
from __future__ import unicode_literals, absolute_import

import os
import json
import logging
import unittest
import warnings
import datetime
import json
import sys
from xml.sax.saxutils import unescape
from mock import patch, PropertyMock, MagicMock, Mock, DEFAULT

import numpy as np

try:
    if 'DISPLAY' not in os.environ:
        import matplotlib; matplotlib.use('Agg')
    import matplotlib
    import matplotlib.pyplot as plt
except ImportError:
    MATPLOTLIB_IS_INSTALLED = False
else:
    MATPLOTLIB_IS_INSTALLED = True

from netCDF4 import Dataset

from nansat import Nansat, Domain, NSR
from nansat.tools import gdal

from nansat.warnings import NansatFutureWarning
from nansat.exceptions import NansatGDALError, WrongMapperError, NansatReadError

from nansat.tests import nansat_test_data as ntd

warnings.simplefilter("always", NansatFutureWarning)
warnings.simplefilter("always", UserWarning)


class NansatTest(unittest.TestCase):
    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_stere = os.path.join(ntd.test_data_path, 'stere.tif')
        self.test_file_complex = os.path.join(ntd.test_data_path, 'complex.nc')
        self.test_file_arctic = os.path.join(ntd.test_data_path, 'arctic.nc')
        self.tmpfilename = os.path.join(ntd.tmp_data_path, 'test.nc')

        if not os.path.exists(self.test_file_gcps):
            raise ValueError('No test data available')

    def tearDown(self):
        try:
            os.unlink(self.tmpfilename)
        except OSError:
            pass

    def test_open_gcps(self):
        with warnings.catch_warnings(record=True) as recorded_warnings:
            n = Nansat(self.test_file_gcps, log_level=40)

        nansat_warning_raised = False
        for rw in recorded_warnings:
            if rw.category == NansatFutureWarning:
                nansat_warning_raised = True
        self.assertFalse(nansat_warning_raised)

        self.assertEqual(type(n), Nansat)
        self.assertEqual(n.vrt.dataset.GetProjection(), '')
        self.assertTrue((n.vrt.dataset.GetGCPProjection().startswith('GEOGCS["WGS 84",')))
        self.assertEqual(n.vrt.dataset.RasterCount, 3)
        self.assertEqual(n.filename, self.test_file_gcps)
        self.assertIsInstance(n.logger, logging.Logger)
        self.assertEqual(n.name, os.path.split(self.test_file_gcps)[1])
        self.assertEqual(n.path, os.path.split(self.test_file_gcps)[0])

    def test_open_with_filename_warning(self):
        with warnings.catch_warnings(record=True) as w:
            n = Nansat(fileName=self.test_file_gcps)
            self.assertEqual(w[0].category, NansatFutureWarning)

    def test_open_with_mappername_warning(self):
        with warnings.catch_warnings(record=True) as w:
            n = Nansat(self.test_file_gcps, mapperName='generic')
            self.assertEqual(w[0].category, NansatFutureWarning)

    def test_open_with_loglevel_warning(self):
        with warnings.catch_warnings(record=True) as w:
            n = Nansat(self.test_file_gcps, logLevel=30)
            self.assertEqual(w[0].category, NansatFutureWarning)

    def test_open_with_domain_warning(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        with warnings.catch_warnings(record=True) as w:
            n = Nansat(domain=d)
            self.assertEqual(len(w), 1)
            self.assertIn('Nansat(domain', str(w[0].message))

    def test_get_time_coverage_start_end(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        n.set_metadata('time_coverage_start', '2016-01-20')
        n.set_metadata('time_coverage_end', '2016-01-21')

        self.assertEqual(type(n.time_coverage_start), datetime.datetime)
        self.assertEqual(type(n.time_coverage_end), datetime.datetime)

    def test_from_domain_array(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        n = Nansat.from_domain(d, np.random.randn(500, 500), {'name': 'band1'})

        self.assertEqual(type(n), Nansat)
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertEqual(n.get_metadata('name', 1), 'band1')
        self.assertEqual(n[1].shape, (500, 500))
        self.assertEqual(n.filename, '')
        self.assertIsInstance(n.logger, logging.Logger)
        self.assertEqual(n.name, '')
        self.assertEqual(n.path, '')

    def test_geolocation_of_exportedNC_vs_original(self):
        ''' Lon/lat in original and exported file should coincide '''
        orig = Nansat(self.test_file_gcps)
        orig.export(self.tmpfilename)

        copy = Nansat(self.tmpfilename)
        lon0, lat0 = orig.get_geolocation_grids()
        lon1, lat1 = copy.get_geolocation_grids()
        np.testing.assert_allclose(lon0, lon1)
        np.testing.assert_allclose(lat0, lat1)

    def test_special_characters_in_exported_metadata(self):
        orig = Nansat(self.test_file_gcps)
        orig.vrt.dataset.SetMetadataItem(str('jsonstring'), json.dumps({'meta1': 'hei',
                                                                        'meta2': 'derr'}))
        orig.export(self.tmpfilename)
        copy = Nansat(self.tmpfilename)
        dd = json.loads(unescape(copy.get_metadata(str('jsonstring')), {'&quot;': '"'}))
        self.assertIsInstance(dd, dict)

    def test_time_coverage_metadata_of_exported_equals_original(self):
        orig = Nansat(self.test_file_gcps)
        orig.set_metadata('time_coverage_start', '2010-01-02T08:49:02.347809')
        orig.set_metadata('time_coverage_end', '2010-01-02T08:50:03.599373')
        orig.export(self.tmpfilename)
        copy = Nansat(self.tmpfilename)

        self.assertEqual(orig.get_metadata('time_coverage_start'),
                copy.get_metadata('time_coverage_start'))
        self.assertEqual(orig.get_metadata('time_coverage_end'),
                copy.get_metadata('time_coverage_end'))

    def test_export_netcdf(self):
        ''' Test export and following import of data with bands containing
        np.nan values
        '''
        n = Nansat(self.test_file_gcps)
        arrNoNaN = np.random.randn(n.shape()[0], n.shape()[1])
        n.add_band(arrNoNaN, {'name': 'testBandNoNaN'})
        arrWithNaN = arrNoNaN.copy()
        arrWithNaN[int(n.shape()[0] / 2 - 10):int(n.shape()[0] / 2 + 10),
                   int(n.shape()[1] / 2 - 10):int(n.shape()[1] / 2 + 10)] = np.nan
        n.add_band(arrWithNaN, {'name': 'testBandWithNaN'})
        n.export(self.tmpfilename)
        exported = Nansat(self.tmpfilename)
        earrNoNaN = exported[str('testBandNoNaN')]
        # Use allclose to allow some roundoff errors
        self.assertTrue(np.allclose(arrNoNaN, earrNoNaN))
        earrWithNaN = exported[str('testBandWithNaN')]
        np.testing.assert_allclose(arrWithNaN, earrWithNaN)

    def test_add_band(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)
        n = Nansat.from_domain(d, log_level=40)
        n.add_band(arr, {'name': 'band1'})

        self.assertEqual(type(n), Nansat)
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertEqual(n.get_metadata('name', 1), 'band1')
        self.assertEqual(n[1].shape, (500, 500))

    def test_add_band_twice(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)
        n = Nansat.from_domain(d, log_level=40)
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

        n = Nansat.from_domain(d, log_level=40)
        n.add_bands([arr, arr],
                    [{'name': 'band1'}, {'name': 'band2'}])

        self.assertIsInstance(n, Nansat)
        self.assertEqual(n.vrt.vrt.vrt, None)
        self.assertIsInstance(n[1], np.ndarray)
        self.assertIsInstance(n[2], np.ndarray)
        self.assertEqual(n.get_metadata('name', 1), 'band1')
        self.assertEqual(n.get_metadata('name', 2), 'band2')

    def test_add_bands_no_parameter(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)

        n = Nansat.from_domain(d, log_level=40)
        n.add_bands([arr, arr])

        self.assertEqual(type(n), Nansat)
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertEqual(type(n[2]), np.ndarray)

    def test_add_subvrts_only_to_one_nansat(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)

        n1 = Nansat.from_domain(d, log_level=40)
        n2 = Nansat.from_domain(d, log_level=40)
        n1.add_band(arr, {'name': 'band1'})

        self.assertEqual(type(n1.vrt.band_vrts), dict)
        self.assertTrue(len(n1.vrt.band_vrts) > 0)
        self.assertEqual(n2.vrt.band_vrts, {})

    def test_bands(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        bands = n.bands()

        self.assertEqual(type(bands), dict)
        self.assertTrue(1 in bands)
        self.assertTrue('name' in bands[1])
        self.assertEqual(bands[1]['name'], 'L_645')

    def test_has_band_if_name_matches(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        hb = n.has_band('L_645')
        self.assertTrue(hb)

    def test_has_band_if_standard_name_matches(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        hb = n.has_band('surface_upwelling_spectral_radiance_in_air_emerging_from_sea_water')
        self.assertTrue(hb)

    def test_write_fig_tif(self):
        n = Nansat(self.test_file_arctic)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_fig_tif.tif')
        n.write_figure(tmpfilename)
        nn = Nansat(tmpfilename)
        # Asserts that the basic georeference (corners in this case) is still
        # present after opening the image
        self.assertTrue(np.allclose(n.get_corners(), nn.get_corners()))

    def test_resize_eResampleAlg_is_given(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        n.resize(pixelsize=500, eResampleAlg=5)

    def test_resize_by_pixelsize(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        n.resize(pixelsize=500, resample_alg=1)

        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_factor(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        n.resize(0.5, resample_alg=1)

        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_width(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        n.resize(width=100, resample_alg=1)

        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_by_height(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        n.resize(height=500, resample_alg=1)

        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_resize(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        n.resize(0.1)
        n.resize(10)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_resize_resize.png')
        n.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(type(n[1]), np.ndarray)

    def test_resize_complex_alg_average(self):
        n = Nansat(self.test_file_complex, log_level=40)
        with warnings.catch_warnings(record=True) as w:
            n.resize(0.5, resample_alg=-1)
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, UserWarning))
            self.assertIn('The imaginary parts of complex numbers '
                            'are lost when resampling by averaging ', str(w[-1].message))

    def test_resize_complex_alg0(self):
        n = Nansat(self.test_file_complex, log_level=40)
        n.resize(0.5, resample_alg=0)

        self.assertTrue(np.any(n[1].imag != 0))

    def test_resize_complex_alg1(self):
        n = Nansat(self.test_file_complex, log_level=40)
        n.resize(0.5, resample_alg=1)

        self.assertTrue(np.any(n[1].imag != 0))

    def test_resize_complex_alg2(self):
        n = Nansat(self.test_file_complex, log_level=40)
        n.resize(0.5, resample_alg=2)

        self.assertTrue(np.any(n[1].imag != 0))

    def test_resize_complex_alg3(self):
        n = Nansat(self.test_file_complex, log_level=40)
        n.resize(0.5, resample_alg=3)

        self.assertTrue(np.any(n[1].imag != 0))

    def test_resize_complex_alg4(self):
        n = Nansat(self.test_file_complex, log_level=40)
        n.resize(0.5, resample_alg=4)

        self.assertTrue(np.any(n[1].imag != 0))

    def test_get_GDALRasterBand(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        b = n.get_GDALRasterBand(1)
        arr = b.ReadAsArray()

        self.assertEqual(type(b), gdal.Band)
        self.assertEqual(type(arr), np.ndarray)

    def test_get_GDALRasterBand_if_bandID_is_given(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        b = n.get_GDALRasterBand(bandID=1)
        arr = b.ReadAsArray()

        self.assertEqual(type(b), gdal.Band)
        self.assertEqual(type(arr), np.ndarray)

    def test_list_bands_true(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        lb = n.list_bands(True)

        self.assertEqual(lb, None)

    def test_list_bands_false(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        lb = n.list_bands(False)

        self.assertEqual(type(lb), str)

    def test_reproject_domain(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        d = Domain(4326, "-te 27 70 30 72 -ts 500 500")
        n.reproject(d)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_domain.png')
        n.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n.shape(), (500, 500))
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertTrue(n.has_band('swathmask'))

    def test_reproject_domain_if_dstDomain_is_given(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        d = Domain(4326, "-te 27 70 30 72 -ts 500 500")
        n.reproject(dstDomain=d)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_domain.png')
        n.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n.shape(), (500, 500))
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertTrue(n.has_band('swathmask'))

    def test_reproject_domain_if_eRasampleAlg_is_given(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        d = Domain(4326, "-te 27 70 30 72 -ts 500 500")
        n.reproject(d,eResampleAlg=0)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_domain.png')
        n.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n.shape(), (500, 500))
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertTrue(n.has_band('swathmask'))

    @patch.multiple(Nansat, filename=DEFAULT, __init__ = Mock(return_value=None))
    def test_property_fileName(self, filename):
        n = Nansat()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            fn = n.fileName
            self.assertTrue( len(w) >= 1 )
            categories = [ww.category for ww in w]
            self.assertIn(NansatFutureWarning, categories)

    @patch.object(Nansat, 'get_corners',
                  return_value=(np.array([0, 0, 360, 360]), np.array([90,-90, 90, -90])))
    def test_reproject_domain_if_source_and_destination_domain_span_entire_lons(self, mock_Nansat):
        n = Nansat(self.test_file_arctic, log_level=40)
        d = Domain(4326, "-te -180 180 60 90 -ts 500 500")
        n.reproject(d)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_domain.png')
        n.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n.shape(), (500, 500))
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertTrue(n.has_band('swathmask'))

    def test_reproject_domain_if_tps_is_given(self):
        n = Nansat(self.test_file_gcps, log_level=40)
        d = Domain(4326, "-te 27 70 30 72 -ts 500 500")
        n.reproject(d, tps=False)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_domain.png')
        n.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n.shape(), (500, 500))
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertTrue(n.has_band('swathmask'))

        n = Nansat(self.test_file_gcps, log_level=40)
        d = Domain(4326, "-te 27 70 30 72 -ts 500 500")
        n.reproject(d, tps=True)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_domain.png')
        n.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n.shape(), (500, 500))
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertTrue(n.has_band('swathmask'))

    def test_reproject_of_complex(self):
        ''' Should return np.nan in areas out of swath '''
        n = Nansat(self.test_file_complex, log_level=40)
        d = Domain(4326, '-te -92.08 26.85 -92.00 26.91 -ts 200 200')
        n.reproject(d)
        b = n[1]

        self.assertTrue(n.has_band('swathmask'))
        self.assertTrue(np.isnan(b[0, 0]))
        self.assertTrue(np.isfinite(b[100, 100]))

    def test_add_band_and_reproject(self):
        ''' Should add band and swath mask
        and return 0 in areas out of swath '''
        n = Nansat(self.test_file_gcps, log_level=40)
        d = Domain(4326, "-te 27 70 30 72 -ts 500 500")
        n.add_band(np.ones(n.shape()))
        n.reproject(d)
        b1 = n[1]
        b4 = n[4]

        self.assertTrue(n.has_band('swathmask'))
        self.assertTrue(b1[0, 0] == 0)
        self.assertTrue(b1[300, 300] > 0)
        self.assertTrue(np.isnan(b4[0, 0]))
        self.assertTrue(b4[300, 300] == 1.)

    def test_reproject_no_addmask(self):
        ''' Should not add swath mask and return 0 in areas out of swath '''
        n = Nansat(self.test_file_gcps, log_level=40)
        d = Domain(4326, '-te -92.08 26.85 -92.00 26.91 -ts 200 200')
        n.reproject(d, addmask=False)
        b = n[1]

        self.assertTrue(not n.has_band('swathmask'))
        self.assertTrue(np.isfinite(b[0, 0]))
        self.assertTrue(np.isfinite(b[100, 100]))

    def test_reproject_stere(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        n2 = Nansat(self.test_file_stere, log_level=40)
        n1.reproject(n2)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_stere.png')
        n1.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n1.shape(), n2.shape())
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_reproject_gcps(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        n2 = Nansat(self.test_file_gcps, log_level=40)
        n1.reproject(n2)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_gcps.png')
        n1.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n1.shape(), n2.shape())
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_reproject_gcps_on_repro_gcps(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        n2 = Nansat(self.test_file_gcps, log_level=40)
        n2.reproject_gcps()
        n1.reproject(n2)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_gcps_on_repro_gcps.png')
        n1.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n1.shape(), n2.shape())
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_reproject_gcps_resize(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        n2 = Nansat(self.test_file_gcps, log_level=40)
        n1.reproject(n2)
        n1.resize(2)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_gcps_resize.png')
        n1.write_figure(tmpfilename, 2, clim='hist')

        self.assertEqual(n1.shape()[0], n2.shape()[0] * 2)
        self.assertEqual(n1.shape()[1], n2.shape()[1] * 2)
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_undo(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        shape1 = n1.shape()
        n1.resize(10)
        n1.undo()
        shape2 = n1.shape()

        self.assertEqual(shape1, shape2)

    def test_write_figure(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_figure.png')
        n1.write_figure(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_figure_band(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_figure_band.png')
        n1.write_figure(tmpfilename, 2)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_figure_clim(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_figure_clim.png')
        n1.write_figure(tmpfilename, 3, clim='hist')

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_figure_legend(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_figure_legend.png')
        n1.write_figure(tmpfilename, 3, clim='hist', legend=True)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_figure_logo(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_figure_logo.png')
        n1.write_figure(tmpfilename, 3, clim='hist',
                        logoFileName=self.test_file_gcps)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_geotiffimage(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_geotiffimage.tif')
        n1.write_geotiffimage(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_write_geotiffimage_if_bandID_is_given(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_geotiffimage.tif')
        n1.write_geotiffimage(tmpfilename, bandID=1)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_get_metadata(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        m = n1.get_metadata()

        self.assertEqual(type(m), dict)
        self.assertTrue('filename' in m)

    def test_get_metadata_key(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        m = n1.get_metadata('filename')

        self.assertEqual(type(m), str)

    def test_get_metadata_wrong_key(self):
        n1 = Nansat(self.test_file_stere, log_level=40)

        with self.assertRaises(ValueError):
            n1.get_metadata('some_crap')

    def test_get_metadata_band_id(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        m = n1.get_metadata(band_id=1)

        self.assertEqual(type(m), dict)
        self.assertTrue('name' in m)

    def test_get_metadata_bandID(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        m = n1.get_metadata(bandID=1)

        self.assertEqual(type(m), dict)
        self.assertTrue('name' in m)

    def test_set_metadata(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        n1.set_metadata('newKey', 'newVal')
        m = n1.get_metadata('newKey')

        self.assertEqual(m, 'newVal')

    def test_set_metadata_band_id(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        n1.set_metadata('newKey', 'newVal', band_id=1)
        m = n1.get_metadata('newKey', 1)

        self.assertEqual(m, 'newVal')

    def test_set_metadata_bandID(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        n1.set_metadata('newKey', 'newVal', bandID=1)
        m = n1.get_metadata('newKey', 1)

        self.assertEqual(m, 'newVal')

    def test_get_band_number(self):
        n1 = Nansat(self.test_file_stere, log_level=40)
        self.assertEqual(n1._get_band_number(1), 1)

    @unittest.skipUnless(MATPLOTLIB_IS_INSTALLED, 'Matplotlib is required')
    def test_get_transect(self):
        plt.switch_backend('agg')
        n1 = Nansat(self.test_file_gcps, log_level=40)
        t = n1.get_transect([[28.31299128, 28.93691525],
                             [70.93709219, 70.69646524]],
                            [str('L_645')])
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_get_transect.png')
        plt.plot(t['lat'], t['L_645'], '.-')
        plt.savefig(tmpfilename)
        plt.close('all')

        self.assertTrue('L_645' in t.dtype.fields)
        self.assertTrue('line' in t.dtype.fields)
        self.assertTrue('pixel' in t.dtype.fields)
        self.assertTrue('lat' in t.dtype.fields)
        self.assertTrue('lon' in t.dtype.fields)
        self.assertEqual(type(t['lat']), np.ndarray)
        self.assertEqual(type(t['lon']), np.ndarray)

    def test_get_transect_outside(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        t = n1.get_transect([[0, 28.31299128], [0, 70.93709219]], [1])

        self.assertTrue('L_645' in t.dtype.fields)
        self.assertTrue('line' in t.dtype.fields)
        self.assertTrue('pixel' in t.dtype.fields)
        self.assertTrue('lat' in t.dtype.fields)
        self.assertTrue('lon' in t.dtype.fields)
        self.assertEqual(type(t['lat']), np.ndarray)
        self.assertEqual(type(t['lon']), np.ndarray)

    def test_get_transect_wrong_points(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        self.assertRaises(ValueError, n1.get_transect, [1, 1], [1])

    def test_get_transect_wrong_band(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        t = n1.get_transect([[0, 28.31299128], [0, 70.93709219]], [10])

        self.assertTrue('line' in t.dtype.fields)
        self.assertTrue('pixel' in t.dtype.fields)
        self.assertTrue('lat' in t.dtype.fields)
        self.assertTrue('lon' in t.dtype.fields)
        self.assertEqual(type(t['lat']), np.ndarray)
        self.assertEqual(type(t['lon']), np.ndarray)

    def test_get_transect_pixlin(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        t = n1.get_transect([[10, 20],
                             [10, 10]],
                             [str('L_645')],
                            lonlat=False)

        self.assertTrue('L_645' in t.dtype.fields)
        self.assertTrue('line' in t.dtype.fields)
        self.assertTrue('pixel' in t.dtype.fields)
        self.assertTrue('lat' in t.dtype.fields)
        self.assertTrue('lon' in t.dtype.fields)
        self.assertEqual(type(t['lat']), np.ndarray)
        self.assertEqual(type(t['lon']), np.ndarray)
        self.assertEqual(len(t['lon']), 11)

    def test_get_transect_data(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        b1 = n1[1]
        t = n1.get_transect([[28.3], [70.9]], [], data=b1)

        self.assertTrue('input' in t.dtype.fields)
        self.assertTrue('L_645' not in t.dtype.fields)
        self.assertTrue('line' in t.dtype.fields)
        self.assertTrue('pixel' in t.dtype.fields)
        self.assertTrue('lat' in t.dtype.fields)
        self.assertTrue('lon' in t.dtype.fields)
        self.assertEqual(type(t['lat']), np.ndarray)
        self.assertEqual(type(t['lon']), np.ndarray)

    @unittest.skipUnless(MATPLOTLIB_IS_INSTALLED and 'DISPLAY' in os.environ, 'Matplotlib is required')
    def test_digitize_points(self):
        ''' shall return empty array in non interactive mode '''
        for backend in matplotlib.rcsetup.interactive_bk:
            # Find a supported interactive backend
            try:
                plt.switch_backend(backend)
                break;
            except:
                pass
        plt.ion()
        n1 = Nansat(self.test_file_gcps, log_level=40)
        points = n1.digitize_points(1)

        self.assertEqual(len(points), 0)
        plt.ioff()

    def test_crop(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        ext = n1.crop(10, 20, 50, 60)

        self.assertEqual(n1.shape(), (60, 50))
        self.assertEqual(ext, (10, 20, 50, 60))
        self.assertEqual(type(n1[1]), np.ndarray)

        n1 = Nansat(self.test_file_gcps, log_level=40)
        ext = n1.crop(0, 0, 200, 200)

        self.assertEqual(n1.shape(), (200, 200))
        self.assertEqual(ext, (0, 0, 200, 200))
        self.assertEqual(type(n1[1]), np.ndarray)


    def test_crop_gcpproj(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        n1.reproject_gcps()
        ext = n1.crop(10, 20, 50, 60)
        xmed = abs(np.median(np.array([gcp.GCPX
                                for gcp in n1.vrt.dataset.GetGCPs()])))
        gcpproj = NSR(n1.vrt.dataset.GetGCPProjection()
                                        ).ExportToProj4().split(' ')[0]

        self.assertTrue(xmed > 360)
        self.assertTrue(gcpproj=='+proj=stere')

    def test_crop_complex(self):
        n1 = Nansat(self.test_file_complex, log_level=40)
        ext = n1.crop(10, 20, 50, 60)

        self.assertEqual(n1.shape(), (60, 50))
        self.assertEqual(ext, (10, 20, 50, 60))
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_crop_no_gcps_arctic(self):
        n1 = Nansat(self.test_file_arctic, log_level=40)
        ext = n1.crop(10, 20, 50, 60)

        self.assertEqual(n1.shape(), (60, 50))
        self.assertEqual(ext, (10, 20, 50, 60))
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_crop_lonlat(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        ext = n1.crop_lonlat([28, 29], [70.5, 71])

        self.assertEqual(n1.shape(), (111, 110))
        self.assertEqual(ext, (31, 89, 110, 111))
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_crop_outside(self):
        n1 = Nansat(self.test_file_gcps, log_level=40)
        self.assertRaises(ValueError, n1.crop_lonlat, [-10, 10], [-10, 10])

    def test_watermask(self):
        ''' if watermask data exists: should fetch array with watermask
            else:                     should raise an error'''
        n1 = Nansat(self.test_file_gcps, log_level=40)
        mod44path = os.getenv('MOD44WPATH')
        if mod44path is not None and os.path.exists(mod44path + '/MOD44W.vrt'):
            wm = n1.watermask()[1]
            self.assertEqual(type(wm), np.ndarray)
            self.assertEqual(wm.shape[0], n1.shape()[0])
            self.assertEqual(wm.shape[1], n1.shape()[1])

    def test_watermask_fail_if_mod44path_is_wrong(self):
        ''' Nansat.watermask should raise an IOError'''
        n1 = Nansat(self.test_file_gcps, log_level=40)
        os.environ['MOD44WPATH'] = '/fakepath'
        self.assertRaises(IOError, n1.watermask)

    def test_watermask_fail_if_mod44path_not_exist(self):
        ''' Nansat.watermask should raise an IOError'''
        n1 = Nansat(self.test_file_gcps, log_level=40)
        del os.environ['MOD44WPATH']
        self.assertRaises(IOError, n1.watermask)

    def test_init_no_arguments(self):
        ''' No arguments should raise ValueError '''
        self.assertRaises(ValueError, Nansat)

    def test_get_item_basic_expressions(self):
        ''' Testing get_item with some basic expressions '''
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        n = Nansat.from_domain(d, np.empty((500, 500)), {'expression': 'np.ones((500, 500))'})
        self.assertIsInstance(n[1], np.ndarray)
        self.assertEqual(n[1].shape, (500, 500))
        self.assertTrue(np.allclose(n[1], np.ones((500, 500))))

    def test_get_item_inf_expressions(self):
        ''' inf should be replaced with nan '''
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        n = Nansat.from_domain(d, log_level=40)
        arr = np.empty((500, 500))
        n.add_band(arr, {'expression': 'np.array([0,1,2,3,np.inf,5,6,7])'})
        self.assertIsInstance(n[1], np.ndarray)
        self.assertTrue(np.isnan(n[1][4]))

    def test_repr_basic(self):
        ''' repr should include some basic elements '''
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        n = Nansat.from_domain(d, log_level=40)
        arr = np.empty((500, 500))
        exp = 'np.array([0,1,2,3,np.inf,5,6,7])'
        n.add_band(arr, {'expression': exp})
        n_repr = repr(n)
        self.assertIn(exp, n_repr, 'The expressions should be in repr')
        self.assertIn('SourceFilename', n_repr)
        self.assertIn('/vsimem/', n_repr)
        self.assertIn('500 x 500', n_repr)
        self.assertIn('Projection:', n_repr)
        self.assertIn('25', n_repr)
        self.assertIn('72', n_repr)
        self.assertIn('35', n_repr)
        self.assertIn('70', n_repr)

    @patch.object(Nansat, 'get_GDALRasterBand')
    def test_getitem(self, mock_Nansat):
        type(mock_Nansat()).GetMetadata = MagicMock(return_value={'a':1})
        type(mock_Nansat()).ReadAsArray = MagicMock(return_value=None)
        with self.assertRaises(NansatGDALError):
            Nansat(self.test_file_stere).__getitem__(1)

    @patch.object(Nansat, 'digitize_points')
    def test_crop_interactive(self, mock_digitize_points):
        mock_digitize_points.return_value=[np.array([[10, 20], [10, 30]])]
        n = Nansat(self.test_file_arctic, log_level=40)
        n.crop_interactive()
        self.assertEqual(n.shape(), (20, 10))

    def test_extend(self):
        n = Nansat(self.test_file_arctic, log_level=40)
        nshape1 = n.shape()
        n.extend(left=10, right=20, top=30, bottom=40)
        be = n[1]
        self.assertEqual(n.shape(), (nshape1[0]+70, nshape1[1]+30))
        self.assertIsInstance(be, np.ndarray)


if __name__ == "__main__":
    unittest.main()
