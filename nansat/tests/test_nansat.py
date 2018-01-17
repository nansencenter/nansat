# ------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified:     Morten Wergeland Hansen, Aleksander Vines
#
# Created:      18.06.2014
# Last modified:16.01.2018 11:00
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------------------
import unittest
import warnings
import os
import datetime
import json
import sys
from xml.sax.saxutils import unescape

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from nansat import Nansat, Domain, NSR
from nansat.tools import gdal, OptionError

import nansat_test_data as ntd
from __builtin__ import int


class NansatTest(unittest.TestCase):
    def setUp(self):
        self.test_file_gcps = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_stere = os.path.join(ntd.test_data_path, 'stere.tif')
        self.test_file_complex = os.path.join(ntd.test_data_path, 'complex.nc')
        self.test_file_arctic = os.path.join(ntd.test_data_path, 'arctic.nc')
        self.tmpfilename = os.path.join(ntd.tmp_data_path, 'test.nc')
        plt.switch_backend('Agg')

        if not os.path.exists(self.test_file_gcps):
            raise ValueError('No test data available')

    def tearDown(self):
        try:
            os.unlink(self.tmpfilename)
        except OSError:
            pass

    def test_open_gcps(self):
        n = Nansat(self.test_file_gcps, logLevel=40)

        self.assertEqual(type(n), Nansat)
        self.assertEqual(n.vrt.dataset.GetProjection(), '')
        self.assertTrue((n.vrt.dataset.GetGCPProjection()
                                            .startswith('GEOGCS["WGS 84",')))

    def test_get_time_coverage_start_end(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        n.set_metadata('time_coverage_start', '2016-01-20')
        n.set_metadata('time_coverage_end', '2016-01-21')

        self.assertEqual(type(n.time_coverage_start),
                         datetime.datetime)
        self.assertEqual(type(n.time_coverage_end),
                         datetime.datetime)

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
        orig.vrt.dataset.SetMetadataItem('jsonstring', json.dumps({'meta1':
                                         'hei', 'meta2': 'derr'}))
        orig.export(self.tmpfilename)
        copy = Nansat(self.tmpfilename)
        dd = json.loads(unescape(copy.get_metadata('jsonstring'), {'&quot;':
                                                                   '"'}))
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
        arrWithNaN[n.shape()[0] / 2 - 10:n.shape()[0] / 2 + 10,
                   n.shape()[1] / 2 - 10:n.shape()[1] / 2 + 10] = np.nan
        n.add_band(arrWithNaN, {'name': 'testBandWithNaN'})
        n.export(self.tmpfilename)
        exported = Nansat(self.tmpfilename)
        earrNoNaN = exported['testBandNoNaN']
        # Use allclose to allow some roundoff errors
        self.assertTrue(np.allclose(arrNoNaN, earrNoNaN))
        earrWithNaN = exported['testBandWithNaN']
        np.testing.assert_allclose(arrWithNaN, earrWithNaN)

    def test_add_band(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)
        n = Nansat(domain=d, logLevel=40)
        n.add_band(arr, {'name': 'band1'})

        self.assertEqual(type(n), Nansat)
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertEqual(n.get_metadata('name', 1), 'band1')
        self.assertEqual(n[1].shape, (500, 500))

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

    def test_add_bands_no_parameter(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        arr = np.random.randn(500, 500)

        n = Nansat(domain=d, logLevel=40)
        n.add_bands([arr, arr])

        self.assertEqual(type(n), Nansat)
        self.assertEqual(type(n[1]), np.ndarray)
        self.assertEqual(type(n[2]), np.ndarray)

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

    def test_export_gcps_to_netcdf(self):
        ''' Should export file with GCPs and write correct bands'''
        n0 = Nansat(self.test_file_gcps, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansat_export_gcps.nc')
        n0.export(tmpfilename)

        ncf = Dataset(tmpfilename)
        self.assertTrue(os.path.exists(tmpfilename))
        self.assertTrue('GCPX' in ncf.variables)
        self.assertTrue('GCPY' in ncf.variables)
        self.assertTrue('GCPPixel' in ncf.variables)
        self.assertTrue('GCPLine' in ncf.variables)

        n1 = Nansat(tmpfilename)
        b0 = n0['L_469']
        b1 = n1['L_469']
        np.testing.assert_allclose(b0, b1)

        lon0, lat0 = n0.get_geolocation_grids()
        lon1, lat1 = n1.get_geolocation_grids()
        np.testing.assert_allclose(lon0, lon1)
        np.testing.assert_allclose(lat0, lat1)

    def test_export_gcps_complex_to_netcdf(self):
        ''' Should export file with GCPs and write correct complex bands'''
        n0 = Nansat(self.test_file_gcps, logLevel=40)
        b0 = n0['L_469']

        n1 = Nansat(domain=n0)
        n1.add_band(b0.astype('complex64'),
                    parameters={'name': 'L_469'})

        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export_gcps_complex.nc')
        n1.export(tmpfilename)

        ncf = Dataset(tmpfilename)
        self.assertTrue(os.path.exists(tmpfilename))
        self.assertTrue('GCPX' in ncf.variables)
        self.assertTrue('GCPY' in ncf.variables)
        self.assertTrue('GCPPixel' in ncf.variables)
        self.assertTrue('GCPLine' in ncf.variables)

        n2 = Nansat(tmpfilename)
        b2 = n2['L_469']
        np.testing.assert_allclose(b0, b2)

        lon0, lat0 = n0.get_geolocation_grids()
        lon2, lat2 = n1.get_geolocation_grids()
        np.testing.assert_allclose(lon0, lon2)
        np.testing.assert_allclose(lat0, lat2)

    def test_export_gtiff(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansat_export.tif')
        n.export(tmpfilename, driver='GTiff')

        self.assertTrue(os.path.exists(tmpfilename))

    def test_export_band(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export_band.tif')
        n.export(tmpfilename, bands=[1], driver='GTiff')
        n = Nansat(tmpfilename, mapperName='generic')

        self.assertTrue(os.path.exists(tmpfilename))
        self.assertEqual(n.vrt.dataset.RasterCount, 1)

    def test_export_band_by_name(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export_band.tif')
        n.export(tmpfilename, bands=['L_645'], driver='GTiff')
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

    def test_export_selected_bands(self):
        n = Nansat(self.test_file_gcps)
        resfile = 'tmp.nc'
        new_band = np.random.randn(n.shape()[0], n.shape()[1])
        n.add_band(new_band, {'name': 'newBand'})
        # Test with band numbers
        n.export(resfile, bands=[4, 2])
        self.assertTrue(os.path.exists(resfile))
        nn = Nansat(resfile)
        self.assertTrue(nn.has_band('newBand'))
        self.assertTrue(nn.has_band('L_555'))
        os.unlink(resfile)
        # Test with band names - not yet implemented
#         n.export(resfile, bands=['newBand', 'L_555'])
#         nn = Nansat(resfile)
#         self.assertTrue(nn.has_band('newBand'))
#         self.assertTrue(nn.has_band('L_555'))
#         os.unlink(resfile)

    def test_export_option(self):
        n = Nansat(self.test_file_arctic)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export_option.nc')
        # Test with band numbers
        n.export(tmpfilename, options='WRITE_LONLAT=YES')
        n.export(tmpfilename + '2', options=['WRITE_LONLAT=NO'])
        nn = Nansat(tmpfilename, mapperName='generic')
        nn2 = Nansat(tmpfilename + '2', mapperName='generic')
        self.assertTrue(nn.has_band('lon'))
        self.assertTrue(nn.has_band('lat'))
        self.assertTrue(nn.has_band('Bristol'))
        self.assertFalse(nn2.has_band('lon'))
        self.assertFalse(nn2.has_band('lat'))
        self.assertTrue(nn2.has_band('Bristol'))

    def test_write_fig_wrong_type_filename(self):
        n = Nansat(self.test_file_arctic)
        with self.assertRaises(OptionError):
            n.write_figure(1.2)
        with self.assertRaises(OptionError):
            n.write_figure(['filename'])
        with self.assertRaises(OptionError):
            n.write_figure({'name': 'filename'})

    def test_write_fig_tif(self):
        n = Nansat(self.test_file_arctic)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_write_fig_tif.tif')
        n.write_figure(tmpfilename)
        nn = Nansat(tmpfilename)
        # Asserts that the basic georeference (corners in this case) is still
        # present after opening the image
        self.assertTrue(np.allclose(n.get_corners(), nn.get_corners()))

    def test_export2thredds_arctic_long_lat(self):
        n = Nansat(self.test_file_arctic, logLevel=40)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds_arctic.nc')
        bands = {
            'Bristol': {'type': '>i2'},
            'Bootstrap': {'type': '>i2'},
            'UMass_AES': {'type': '>i2'},
        }
        n.export2thredds(tmpfilename, bands, time=datetime.datetime(2016, 1, 20))

        self.assertTrue(os.path.exists(tmpfilename))
        g = gdal.Open(tmpfilename)
        metadata = g.GetMetadata_Dict()

        # Test that the long/lat values are set aproximately correct
        ncg = 'NC_GLOBAL#'
        if (metadata.get(ncg + 'easternmost_longitude') == None):
            ncg = ''
        easternmost_longitude = metadata.get(ncg + 'easternmost_longitude')
        self.assertTrue(float(easternmost_longitude) > 179,
                        'easternmost_longitude is wrong:' +
                        easternmost_longitude)
        westernmost_longitude = metadata.get(ncg + 'westernmost_longitude')
        self.assertTrue(float(westernmost_longitude) < -179,
                        'westernmost_longitude is wrong:' +
                        westernmost_longitude)
        northernmost_latitude = metadata.get(ncg + 'northernmost_latitude')
        self.assertTrue(float(northernmost_latitude) > 89.999,
                        'northernmost_latitude is wrong:' +
                        northernmost_latitude)
        southernmost_latitude = metadata.get(ncg + 'southernmost_latitude')
        self.assertTrue(float(southernmost_latitude) < 54,
                        'southernmost_latitude is wrong:' +
                        southernmost_latitude)
        self.assertTrue(float(southernmost_latitude) > 53,
                        'southernmost_latitude is wrong:' +
                        southernmost_latitude)

    def test_dont_export2thredds_gcps(self):
        n = Nansat(self.test_file_gcps, logLevel=40)
        n2 = Nansat(domain=n)
        n.add_band(np.ones(n2.shape(), np.float32))
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds.nc')
        self.assertRaises(OptionError, n2.export2thredds, tmpfilename,
                          ['L_645'])

    def test_export2thredds_longlat_list(self):
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   "-te 27 70 31 72 -ts 200 200")
        n = Nansat(domain=d)
        n.add_band(np.ones(d.shape(), np.float32),
                   parameters={'name': 'L_469'})
        n.set_metadata('time_coverage_start', '2016-01-19')

        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds_longlat.nc')
        n.export2thredds(tmpfilename, ['L_469'])
        ncI = Dataset(tmpfilename, 'r')
        ncIVar = ncI.variables['L_469']
        self.assertTrue(ncIVar.grid_mapping in ncI.variables.keys())

    def test_export2thredds_longlat_dict(self):
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   "-te 27 70 31 72 -ts 200 200")
        n = Nansat(domain=d)
        n.add_band(np.ones(d.shape(), np.float32),
                   parameters={'name': 'L_469'})
        n.set_metadata('time_coverage_start', '2016-01-19')

        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_export2thredds_longlat.nc')
        n.export2thredds(tmpfilename, {'L_469': {'type': '>i1'}})
        ncI = Dataset(tmpfilename, 'r')
        ncIVar = ncI.variables['L_469']
        self.assertTrue(ncIVar.grid_mapping in ncI.variables.keys())
        self.assertEqual(ncIVar[:].dtype, np.int8)

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

            self.assertTrue(len(w) == 1)
            self.assertTrue(issubclass(w[-1].category, UserWarning))
            self.assertTrue(
                        'The imaginary parts of complex numbers '
                        'are lost when resampling by averaging '
                        in str(w[-1].message)
                        )

    def test_resize_complex_alg0(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=0)

        self.assertTrue(np.any(n[1].imag != 0))

    def test_resize_complex_alg1(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=1)

        self.assertTrue(np.any(n[1].imag != 0))

    def test_resize_complex_alg2(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=2)

        self.assertTrue(np.any(n[1].imag != 0))

    def test_resize_complex_alg3(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=3)

        self.assertTrue(np.any(n[1].imag != 0))

    def test_resize_complex_alg4(self):
        n = Nansat(self.test_file_complex, logLevel=40)
        n.resize(0.5, eResampleAlg=4)

        self.assertTrue(np.any(n[1].imag != 0))

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
        self.assertTrue(n.has_band('swathmask'))

    def test_reproject_of_complex(self):
        ''' Should return np.nan in areas out of swath '''
        n = Nansat(self.test_file_complex, logLevel=40)
        d = Domain(4326, '-te -92.08 26.85 -92.00 26.91 -ts 200 200')
        n.reproject(d)
        b = n[1]

        self.assertTrue(n.has_band('swathmask'))
        self.assertTrue(np.isnan(b[0, 0]))
        self.assertTrue(np.isfinite(b[100, 100]))

    def test_add_band_and_reproject(self):
        ''' Should add band and swath mask
        and return 0 in areas out of swath '''
        n = Nansat(self.test_file_gcps, logLevel=40)
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
        n = Nansat(self.test_file_complex, logLevel=40)
        d = Domain(4326, '-te -92.08 26.85 -92.00 26.91 -ts 200 200')
        n.reproject(d, addmask=False)
        b = n[1]

        self.assertTrue(not n.has_band('swathmask'))
        self.assertTrue(np.isfinite(b[0, 0]))
        self.assertTrue(np.isfinite(b[100, 100]))

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

    def test_reproject_gcps_on_repro_gcps(self):
        n1 = Nansat(self.test_file_stere, logLevel=40)
        n2 = Nansat(self.test_file_gcps, logLevel=40)
        n2.reproject_GCPs()
        n1.reproject(n2)
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'nansat_reproject_gcps_on_repro_gcps.png')
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

        self.assertEqual(n1.shape()[0], n2.shape()[0] * 2)
        self.assertEqual(n1.shape()[1], n2.shape()[1] * 2)
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

    def test_write_figure_legend(self):
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

        with self.assertRaises(OptionError):
            n1.get_metadata('some_crap')

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
        t = n1.get_transect([[28.31299128, 28.93691525],
                             [70.93709219, 70.69646524]],
                            ['L_645'])
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
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        t = n1.get_transect([[0, 28.31299128], [0, 70.93709219]], [1])

        self.assertTrue('L_645' in t.dtype.fields)
        self.assertTrue('line' in t.dtype.fields)
        self.assertTrue('pixel' in t.dtype.fields)
        self.assertTrue('lat' in t.dtype.fields)
        self.assertTrue('lon' in t.dtype.fields)
        self.assertEqual(type(t['lat']), np.ndarray)
        self.assertEqual(type(t['lon']), np.ndarray)

    def test_get_transect_wrong_points(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        self.assertRaises(OptionError, n1.get_transect, [1, 1], [1])

    def test_get_transect_wrong_band(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        t = n1.get_transect([[0, 28.31299128], [0, 70.93709219]], [10])

        self.assertTrue('line' in t.dtype.fields)
        self.assertTrue('pixel' in t.dtype.fields)
        self.assertTrue('lat' in t.dtype.fields)
        self.assertTrue('lon' in t.dtype.fields)
        self.assertEqual(type(t['lat']), np.ndarray)
        self.assertEqual(type(t['lon']), np.ndarray)

    def test_get_transect_pixlin(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        t = n1.get_transect([[10, 20],
                             [10, 10]],
                            ['L_645'],
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
        n1 = Nansat(self.test_file_gcps, logLevel=40)
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
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        points = n1.digitize_points(1)

        self.assertEqual(len(points), 0)
        plt.ioff()

    def test_crop(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        ext = n1.crop(10, 20, 50, 60)

        self.assertEqual(n1.shape(), (60, 50))
        self.assertEqual(ext, (10, 20, 50, 60))
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_crop_gcpproj(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        n1.reproject_GCPs()
        ext = n1.crop(10, 20, 50, 60)
        xmed = abs(np.median(np.array([gcp.GCPX
                                for gcp in n1.vrt.dataset.GetGCPs()])))
        gcpproj = NSR(n1.vrt.dataset.GetGCPProjection()
                                        ).ExportToProj4().split(' ')[0]
        
        self.assertTrue(xmed > 360)
        self.assertTrue(gcpproj=='+proj=stere')

    def test_crop_complex(self):
        n1 = Nansat(self.test_file_complex, logLevel=40)
        ext = n1.crop(10, 20, 50, 60)

        self.assertEqual(n1.shape(), (60, 50))
        self.assertEqual(ext, (10, 20, 50, 60))
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_crop_no_gcps_arctic(self):
        n1 = Nansat(self.test_file_arctic, logLevel=40)
        ext = n1.crop(10, 20, 50, 60)

        self.assertEqual(n1.shape(), (60, 50))
        self.assertEqual(ext, (10, 20, 50, 60))
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_crop_lonlat(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        ext = n1.crop_lonlat([28, 29], [70.5, 71])

        self.assertEqual(n1.shape(), (111, 110))
        self.assertEqual(ext, (31, 89, 110, 111))
        self.assertEqual(type(n1[1]), np.ndarray)

    def test_crop_outside(self):
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        self.assertRaises(OptionError, n1.crop_lonlat, [-10, 10], [-10, 10])

    def test_watermask(self):
        ''' if watermask data exists: should fetch array with watermask
            else:                     should raise an error'''
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        mod44path = os.getenv('MOD44WPATH')
        if mod44path is not None and os.path.exists(mod44path + '/MOD44W.vrt'):
            wm = n1.watermask()[1]
            self.assertEqual(type(wm), np.ndarray)
            self.assertEqual(wm.shape[0], n1.shape()[0])
            self.assertEqual(wm.shape[1], n1.shape()[1])

    def test_watermask_fail(self):
        ''' Nansat.watermask should raise an IOError'''
        n1 = Nansat(self.test_file_gcps, logLevel=40)
        os.environ['MOD44WPATH'] = '/fakepath'
        self.assertRaises(IOError, n1.watermask)

    def test_init_no_arguments(self):
        ''' No arguments should raise OptionError '''
        self.assertRaises(OptionError, Nansat)

    def test_get_item_basic_expressions(self):
        ''' Testing get_item with some basic expressions '''
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        n = Nansat(domain=d, logLevel=40)
        arr = np.empty((500, 500))
        n.add_band(arr, {'expression': '1+1'})
        n.add_band(arr, {'expression': 'np.random.randn(500, 500)'})
        self.assertIsInstance(n[1], int)
        self.assertIsInstance(n[2], np.ndarray)
        self.assertEqual(n[1], 2)
        self.assertEqual(len(n[2]), 500)
        self.assertEqual(len(n[2][0]), 500)

    def test_get_item_inf_expressions(self):
        ''' inf should be replaced with nan '''
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        n = Nansat(domain=d, logLevel=40)
        arr = np.empty((500, 500))
        n.add_band(arr, {'expression': 'np.array([0,1,2,3,np.inf,5,6,7])'})
        self.assertIsInstance(n[1], np.ndarray)
        self.assertTrue(np.isnan(n[1][4]))

    def test_repr_basic(self):
        ''' repr should include some basic elements '''
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        n = Nansat(domain=d, logLevel=40)
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

    def test_export_netcdf_complex_remove_meta(self):
        ''' Test export of complex data with pixelfunctions
        '''
        n = Nansat(self.test_file_complex)
        self.assertEqual(n.get_metadata('PRODUCT_TYPE'), 'SLC')
        n.export(self.tmpfilename, rmMetadata=['PRODUCT_TYPE'])
        exported = Nansat(self.tmpfilename)
        with self.assertRaises(OptionError):
            exported.get_metadata('PRODUCT_TYPE')
        self.assertTrue((n[1] == exported[1]).any())

    def test_export_netcdf_arctic(self):
        ''' Test export of the arctic data without GCPS
        '''
        n = Nansat(self.test_file_arctic)
        n.export(self.tmpfilename)
        exported = Nansat(self.tmpfilename)
        self.assertTrue((n[1] == exported[1]).any())
        self.assertTrue((n[2] == exported[2]).any())
        self.assertTrue((n[3] == exported[3]).any())

if __name__ == "__main__":
    unittest.main()
