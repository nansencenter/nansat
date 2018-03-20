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
from __future__ import print_function, absolute_import, division
import os
import sys
import json
import unittest
import warnings
import datetime
from xml.sax.saxutils import unescape
from mock import patch, PropertyMock, Mock, MagicMock, DEFAULT

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

import gdal
from netCDF4 import Dataset

from nansat import Nansat, Domain, NSR
from nansat.warnings import NansatFutureWarning
from nansat.tests.nansat_test_base import NansatTestBase

warnings.simplefilter("always", NansatFutureWarning)
warnings.simplefilter("always", UserWarning)


class ExporterTest(NansatTestBase):

    def test_geolocation_of_exportedNC_vs_original(self):
        """ Lon/lat in original and exported file should coincide """
        orig = Nansat(self.test_file_gcps, mapper=self.default_mapper)
        orig.export(self.tmp_filename)

        copy = Nansat(self.tmp_filename, mapper=self.default_mapper)
        lon0, lat0 = orig.get_geolocation_grids()
        lon1, lat1 = copy.get_geolocation_grids()
        np.testing.assert_allclose(lon0, lon1)
        np.testing.assert_allclose(lat0, lat1)

    def test_special_characters_in_exported_metadata(self):
        orig = Nansat(self.test_file_gcps, mapper=self.default_mapper)
        orig.vrt.dataset.SetMetadataItem('jsonstring', json.dumps({'meta1':
                                         'hei', 'meta2': 'derr'}))
        orig.export(self.tmp_filename)
        copy = Nansat(self.tmp_filename, mapper=self.default_mapper)
        dd = json.loads(unescape(copy.get_metadata('jsonstring'), {'&quot;':
                                                                   '"'}))
        self.assertIsInstance(dd, dict)

    def test_time_coverage_metadata_of_exported_equals_original(self):
        orig = Nansat(self.test_file_gcps, mapper=self.default_mapper)
        orig.set_metadata('time_coverage_start', '2010-01-02T08:49:02.347809')
        orig.set_metadata('time_coverage_end', '2010-01-02T08:50:03.599373')
        orig.export(self.tmp_filename)
        copy = Nansat(self.tmp_filename, mapper=self.default_mapper)

        self.assertEqual(orig.get_metadata('time_coverage_start'),
                copy.get_metadata('time_coverage_start'))
        self.assertEqual(orig.get_metadata('time_coverage_end'),
                copy.get_metadata('time_coverage_end'))

    def test_export_netcdf(self):
        """ Test export and following import of data with bands containing np.nan values """
        n = Nansat(self.test_file_gcps, mapper=self.default_mapper)
        arrNoNaN = np.random.randn(n.shape()[0], n.shape()[1])
        n.add_band(arrNoNaN, {'name': 'testBandNoNaN'})
        arrWithNaN = arrNoNaN.copy()
        arrWithNaN[int(n.shape()[0] / 2.) - 10:int(n.shape()[0] / 2 + 10),
                   int(n.shape()[1] / 2.) - 10:int(n.shape()[1] / 2 + 10)] = np.nan
        n.add_band(arrWithNaN, {'name': 'testBandWithNaN'})
        n.export(self.tmp_filename)
        exported = Nansat(self.tmp_filename, mapper=self.default_mapper)
        earrNoNaN = exported['testBandNoNaN']
        # Use allclose to allow some roundoff errors
        self.assertTrue(np.allclose(arrNoNaN, earrNoNaN))
        earrWithNaN = exported['testBandWithNaN']
        np.testing.assert_allclose(arrWithNaN, earrWithNaN)

    def test_export_gcps_filename_warning(self):
        """ Should export file with GCPs and write correct bands"""
        n0 = Nansat(self.test_file_gcps, log_level=40, mapper=self.default_mapper)
        tmpfilename = os.path.join(self.tmp_data_path, 'temp.nc')
        with warnings.catch_warnings(record=True) as w:
            n0.export(fileName=tmpfilename)
            self.assertEqual(len(w), 1)
            self.assertIn('Nansat.export(fileName', str(w[0].message))

    def test_export_gcps_to_netcdf(self):
        """ Should export file with GCPs and write correct bands"""
        n0 = Nansat(self.test_file_gcps, log_level=40, mapper=self.default_mapper)
        tmpfilename = os.path.join(self.tmp_data_path, 'nansat_export_gcps.nc')
        n0.export(tmpfilename)

        ncf = Dataset(tmpfilename)
        self.assertTrue(os.path.exists(tmpfilename))
        self.assertTrue('GCPX' in ncf.variables)
        self.assertTrue('GCPY' in ncf.variables)
        self.assertTrue('GCPPixel' in ncf.variables)
        self.assertTrue('GCPLine' in ncf.variables)

        n1 = Nansat(tmpfilename, mapper=self.default_mapper)
        b0 = n0['L_469']
        b1 = n1['L_469']
        np.testing.assert_allclose(b0, b1)

        lon0, lat0 = n0.get_geolocation_grids()
        lon1, lat1 = n1.get_geolocation_grids()
        np.testing.assert_allclose(lon0, lon1)
        np.testing.assert_allclose(lat0, lat1)

    def test_export_gcps_complex_to_netcdf(self):
        """ Should export file with GCPs and write correct complex bands"""
        n0 = Nansat(self.test_file_gcps, log_level=40, mapper=self.default_mapper)
        b0 = n0['L_469']

        n1 = Nansat.from_domain(n0)
        n1.add_band(b0.astype('complex64'), parameters={'name': 'L_469'})

        tmpfilename = os.path.join(self.tmp_data_path, 'nansat_export_gcps_complex.nc')
        n1.export(tmpfilename)

        ncf = Dataset(tmpfilename)
        self.assertTrue(os.path.exists(tmpfilename))
        self.assertTrue('GCPX' in ncf.variables)
        self.assertTrue('GCPY' in ncf.variables)
        self.assertTrue('GCPPixel' in ncf.variables)
        self.assertTrue('GCPLine' in ncf.variables)

        n2 = Nansat(tmpfilename, mapper=self.default_mapper)
        b2 = n2['L_469']

        lon0, lat0 = n0.get_geolocation_grids()
        lon2, lat2 = n1.get_geolocation_grids()
        np.testing.assert_allclose(lon0, lon2)
        np.testing.assert_allclose(lat0, lat2)

    def test_export_gtiff(self):
        n = Nansat(self.test_file_gcps, log_level=40, mapper=self.default_mapper)
        tmpfilename = os.path.join(self.tmp_data_path, 'nansat_export.tif')
        n.export(tmpfilename, driver='GTiff')

        self.assertTrue(os.path.exists(tmpfilename))

    def test_export_band(self):
        n = Nansat(self.test_file_gcps, log_level=40, mapper=self.default_mapper)
        tmpfilename = os.path.join(self.tmp_data_path,
                                   'nansat_export_band.tif')
        n.export(tmpfilename, bands=[1], driver='GTiff')
        n = Nansat(tmpfilename, mapper=self.default_mapper)

        self.assertTrue(os.path.exists(tmpfilename))
        self.assertEqual(n.vrt.dataset.RasterCount, 1)

    def test_export_band_by_name(self):
        n = Nansat(self.test_file_gcps, log_level=40, mapper=self.default_mapper)
        tmpfilename = os.path.join(self.tmp_data_path,
                                   'nansat_export_band.tif')
        n.export(tmpfilename, bands=['L_645'], driver='GTiff')
        n = Nansat(tmpfilename, mapper=self.default_mapper)

        self.assertTrue(os.path.exists(tmpfilename))
        self.assertEqual(n.vrt.dataset.RasterCount, 1)

    def test_reproject_and_export_band(self):
        n1 = Nansat(self.test_file_gcps, log_level=40, mapper=self.default_mapper)
        n2 = Nansat(self.test_file_stere, log_level=40, mapper=self.default_mapper)
        n1.reproject(n2)
        tmpfilename = os.path.join(self.tmp_data_path,
                                   'nansat_reproject_export_band.nc')
        n1.export(tmpfilename, bands=[1])

        n = Nansat(tmpfilename, mapper=self.default_mapper)
        self.assertTrue(os.path.exists(tmpfilename))
        self.assertEqual(n.vrt.dataset.RasterCount, 1)

    def test_export_selected_bands(self):
        n = Nansat(self.test_file_gcps, mapper=self.default_mapper)
        resfile = 'tmp.nc'
        new_band = np.random.randn(n.shape()[0], n.shape()[1])
        n.add_band(new_band, {'name': 'newBand'})
        # Test with band numbers
        n.export(resfile, bands=[4, 2])
        self.assertTrue(os.path.exists(resfile))
        nn = Nansat(resfile, mapper=self.default_mapper)
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
        n = Nansat(self.test_file_arctic, mapper=self.default_mapper)
        tmpfilename = os.path.join(self.tmp_data_path,
                                   'nansat_export_option.nc')
        # Test with band numbers
        n.export(tmpfilename, options='WRITE_LONLAT=YES')
        n.export(tmpfilename + '2', options=['WRITE_LONLAT=NO'])
        nn = Nansat(tmpfilename, mapper=self.default_mapper)
        nn2 = Nansat(tmpfilename + '2', mapper=self.default_mapper)
        self.assertTrue(nn.has_band('lon'))
        self.assertTrue(nn.has_band('lat'))
        self.assertTrue(nn.has_band('Bristol'))
        self.assertFalse(nn2.has_band('lon'))
        self.assertFalse(nn2.has_band('lat'))
        self.assertTrue(nn2.has_band('Bristol'))

    def test_export2thredds_arctic_long_lat(self):
        n = Nansat(self.test_file_arctic, mapper=self.default_mapper, log_level=40)
        tmpfilename = os.path.join(self.tmp_data_path,
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

		# GDAL behaves differently:
		# Windows: nc-attributes are accessible without 'NC_GLOBAL#' prefix
		# Linux: nc-attributes are accessible only with 'NC_GLOBAL#' prefix
        # OSX: ?
        # Therefore we have to add NC_GLOBAL# and test if such metadata exists
        nc_prefix = 'NC_GLOBAL#'
        if not nc_prefix + 'easternmost_longitude' in metadata:
            nc_prefix = ''
        self.assertIn(nc_prefix + 'easternmost_longitude', metadata)

        # Test that the long/lat values are set correctly
        test_metadata_keys = ['easternmost_longitude', 'westernmost_longitude',
                              'northernmost_latitude', 'southernmost_latitude']
        test_metadata_min = [179, -180, 89.9, 53]
        test_metadata_max = [180, -179, 90, 54]
        for i, test_metadata_key in enumerate(test_metadata_keys):
            medata_value = float(metadata[nc_prefix + test_metadata_key])
            self.assertTrue(medata_value >= test_metadata_min[i],
                            '%s is wrong: %f'%(test_metadata_key, medata_value))
            self.assertTrue(medata_value <= test_metadata_max[i],
                            '%s is wrong: %f'%(test_metadata_key, medata_value))

    def test_dont_export2thredds_gcps(self):
        n = Nansat(self.test_file_gcps, log_level=40, mapper=self.default_mapper)
        n2 = Nansat.from_domain(n)
        n.add_band(np.ones(n2.shape(), np.float32))
        tmpfilename = os.path.join(self.tmp_data_path,
                                   'nansat_export2thredds.nc')
        self.assertRaises(ValueError, n2.export2thredds, tmpfilename,
                          ['L_645'])

    def test_export2thredds_longlat_list(self):
        n = Nansat(self.test_file_gcps, log_level=40, mapper=self.default_mapper)
        with self.assertRaises(ValueError):
            n.export2thredds('aa', ['L_469'])

    def test_export2thredds_longlat_dict(self):
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   "-te 27 70 31 72 -ts 200 200")
        n = Nansat.from_domain(d)
        n.add_band(np.ones(d.shape(), np.float32),
                   parameters={'name': 'L_469'})
        n.set_metadata('time_coverage_start', '2016-01-19')

        tmpfilename = os.path.join(self.tmp_data_path,
                                   'nansat_export2thredds_longlat.nc')
        n.export2thredds(tmpfilename, {'L_469': {'type': '>i1'}})
        ncI = Dataset(tmpfilename, 'r')
        ncIVar = ncI.variables['L_469']
        self.assertTrue(ncIVar.grid_mapping in ncI.variables.keys())
        self.assertEqual(ncIVar[:].dtype, np.int8)


    def test_export_netcdf_complex_remove_meta(self):
        n = Nansat(self.test_file_complex, mapper=self.default_mapper)
        self.assertEqual(n.get_metadata('PRODUCT_TYPE'), 'SLC')
        with warnings.catch_warnings(record=True) as recorded_warnings:
            n.export(self.tmp_filename, rmMetadata=['PRODUCT_TYPE'])
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)
        exported = Nansat(self.tmp_filename, mapper=self.default_mapper)
        with self.assertRaises(ValueError):
            exported.get_metadata('PRODUCT_TYPE')
        self.assertTrue((n[1] == exported[1]).any())

    def test_export_netcdf_arctic(self):
        n = Nansat(self.test_file_arctic, mapper=self.default_mapper)
        n.export(self.tmp_filename)
        exported = Nansat(self.tmp_filename, mapper=self.default_mapper)
        self.assertTrue((n[1] == exported[1]).any())
        self.assertTrue((n[2] == exported[2]).any())
        self.assertTrue((n[3] == exported[3]).any())

    def test_export_netcdf_arctic_hardcopy(self):
        n = Nansat(self.test_file_arctic, mapper=self.default_mapper)
        n.export(self.tmp_filename, hardcopy=True)
        exported = Nansat(self.tmp_filename, mapper=self.default_mapper)
        self.assertTrue((n[1] == exported[1]).any())
        self.assertTrue((n[2] == exported[2]).any())
        self.assertTrue((n[3] == exported[3]).any())

    @patch('nansat.exporter.VRT._add_geolocation')
    def test_export_add_geoloc(self, mock_add_geolocation):
        n = Nansat(self.test_file_arctic, mapper=self.default_mapper)
        with warnings.catch_warnings(record=True) as recorded_warnings:
            n.export(self.tmp_filename, addGeoloc=True)
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)
        self.assertTrue(mock_add_geolocation.called)

    def test_export_add_gcps(self):
        n = Nansat(self.test_file_arctic, mapper=self.default_mapper)
        with warnings.catch_warnings(record=True) as recorded_warnings:
            n.export(self.tmp_filename, addGCPs=True, bottomup=True)
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)
            self.assertEqual(recorded_warnings[1].category, NansatFutureWarning)

    def test_export2thredds_rmmetadata(self):
        n = Nansat(self.test_file_arctic, mapper=self.default_mapper, log_level=40)
        with warnings.catch_warnings(record=True) as recorded_warnings:
            n.export2thredds(self.tmp_filename, {'Bristol': {'type': '>i2'}},
                            time=datetime.datetime(2016, 1, 20),
                            rmMetadata=['description'])
            self.assertEqual(recorded_warnings[0].category, NansatFutureWarning)

if __name__ == "__main__":
    unittest.main()
