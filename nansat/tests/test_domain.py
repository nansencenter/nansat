#------------------------------------------------------------------------------
# Name:         test_domain.py
# Purpose:      Test the Domain class
#
# Author:       Anton Korosov
#
# Created:      29.09.2014
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
import numpy as np

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
except ImportError:
    BASEMAP_LIB_EXISTS = False
else:
    BASEMAP_LIB_EXISTS = True

from nansat.nsr import NSR
from nansat.domain import Domain
from nansat.tools import OptionError, gdal, ogr, ProjectionError
from nansat.figure import Image
import sys
import nansat_test_data as ntd


class DomainTest(unittest.TestCase):
    def setUp(self):
        self.test_file_raw_proj = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_projected = os.path.join(ntd.test_data_path, 'stere.tif')
        if BASEMAP_LIB_EXISTS:
            plt.switch_backend('Agg')

        if not os.path.exists(self.test_file):
            raise ValueError('No test data available')

    def test_init_from_strings(self):
        d = Domain("+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   "-te 25 70 35 72 -ts 2000 2000")

        self.assertEqual(type(d), Domain)

    def test_init_from_epsg_and_te_string(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")

        self.assertEqual(type(d), Domain)

    def test_init_from_epsg_and_lle_string(self):
        d = Domain(4326, "-lle 25 70 35 72 -ts 500 500")

        self.assertEqual(type(d), Domain)

    def test_init_from_lonlat(self):
        lat, lon = np.mgrid[-90:90:0.5, -180:180:0.5]
        d = Domain(lon=lon, lat=lat)

        self.assertEqual(type(d), Domain)
        self.assertEqual(d.shape(), lat.shape)

    def test_init_from_GDALDataset(self):
        ds = gdal.Open(self.test_file)
        d = Domain(ds=ds)

        self.assertEqual(type(d), Domain)

    def test_dont_init_from_invalid(self):
        self.assertRaises(OptionError, Domain)
        self.assertRaises(OptionError, Domain, None)
        with self.assertRaises(OptionError):
            Domain(ds=gdal.Open(self.test_file),
                   srs="+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   ext="-te 25 70 35 72 -ts 2000 2000")
        with self.assertRaises(ProjectionError):
            Domain(ds=gdal.Open(self.test_file),
                   srs="unmatched srs")

    def test_init_use_AutoCreateWarpedVRT_to_determine_bounds(self):
        d = Domain(ds=gdal.Open(self.test_file),
                   srs="+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs")
        self.assertEqual(type(d), Domain)

    def test_write_kml(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(ntd.tmp_data_path, 'domain_write_kml.kml')
        d.write_kml(kmlFileName=tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    def test_get_geolocation_grids(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.get_geolocation_grids()

        self.assertEqual(type(lon), np.ndarray)
        self.assertEqual(type(lat), np.ndarray)
        self.assertEqual(lat.shape, (500, 500))

    def test_get_border_wkt(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        bwkt = d.get_border_wkt()

        self.assertEqual(type(bwkt), str)
        self.assertTrue('POLYGON' in bwkt)

    def test_get_border_geometry(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        geom = d.get_border_geometry()

        self.assertEqual(type(geom), ogr.Geometry)
        # the below test doesn't work in Travis
        # probably some libs are missing in default anaconda install
        # self.assertTrue(geom.IsValid())

    def test_get_corners(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.get_corners()

        self.assertTrue(all(lon - [25., 25., 35., 35.] < 0.01))
        self.assertTrue(all(lat - [72., 70., 72., 70.] < 0.01))

    def test_get_pixelsize_meters(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        x, y = d.get_pixelsize_meters()

        self.assertTrue(x - 444 < 1)
        self.assertTrue(y - 723 < 1)

    def test_transform_points(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.transform_points([1, 2, 3], [1, 2, 3])

        self.assertEqual(type(lon), np.ndarray)
        self.assertEqual(type(lat), np.ndarray)

    def test_transform_points_inverse(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        x, y = d.transform_points([25, 26, 27], [70, 71, 72], 1)

        self.assertTrue(all(np.round(x) == [0, 50, 100]))
        self.assertTrue(all(np.round(y) == [500, 250, 0]))

    def test_transform_points_dstsrs(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.transform_points([1, 2, 3], [1, 2, 3],
        dstSRS=NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=75 +lon_0=10 +no_defs'))

        self.assertEqual(type(lon), np.ndarray)
        self.assertEqual(type(lat), np.ndarray)

    def test_azimuth_y(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        au = d.azimuth_y()

        self.assertEqual(np.round(au[0, 0]), 0)
        self.assertEqual(np.round(au[10, 10]), 0)

    def test_shape(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        shape = d.shape()

        self.assertEqual(shape, (500, 500))

    @unittest.skipUnless(BASEMAP_LIB_EXISTS, 'Basemap is required')
    def test_write_map(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(ntd.tmp_data_path, 'domain_write_map.png')
        d.write_map(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (50, 50))

    @unittest.skipUnless(BASEMAP_LIB_EXISTS, 'Basemap is required')
    def test_write_map_dpi100(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'domain_write_map_dpi100.png')
        d.write_map(tmpfilename, dpi=100)

        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (100, 100))

    @unittest.skipUnless(BASEMAP_LIB_EXISTS, 'Basemap is required')
    def test_write_map_labels(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'domain_write_map_labels.png')
        d.write_map(tmpfilename,
                    merLabels=[False, False, False, True],
                    parLabels=[True, False, False, False])

        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()

    def test_reproject_GCPs(self):
        ds = gdal.Open(self.test_file)
        d = Domain(ds=ds)
        d.reproject_GCPs('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=75 +lon_0=10 +no_defs')
        gcp = d.vrt.dataset.GetGCPs()[0]

        self.assertTrue(gcp.GCPX > 636161)
        self.assertTrue(gcp.GCPY < -288344)

    def test_reproject_GCPs_auto(self):
        ds = gdal.Open(self.test_file)
        d = Domain(ds=ds)
        d.reproject_GCPs()
        
        gcpproj = NSR(d.vrt.dataset.GetGCPProjection())
        self.assertEqual(gcpproj.GetAttrValue('PROJECTION'),
                        'Stereographic')

    def test_overlaps_contains(self):
        Bergen = Domain(4326, "-te 5 60 6 61 -ts 500 500")
        WestCoast = Domain(4326, "-te 1 58 6 64 -ts 500 500")
        Norway = Domain(4326, "-te 3 55 30 72 -ts 500 500")
        Paris = Domain(4326, "-te 2 48 3 49 -ts 500 500")
        self.assertTrue(Bergen.overlaps(Norway))
        self.assertTrue(Norway.contains(Bergen))
        self.assertFalse(Bergen.contains(Norway))
        self.assertTrue(Norway.overlaps(WestCoast))
        self.assertFalse(Norway.contains(WestCoast))
        self.assertFalse(Paris.overlaps(Norway))
        self.assertFalse(Paris.contains(Norway))

    def test_add_to_dict(self):
        input_1 = ['-te', '5', '60', '6', '61.1']
        output_1 = {'te': [5., 60., 6., 61.1]}
        key_1, extent_1 = Domain._add_to_dict(dict(), input_1)
        self.assertIsInstance(extent_1, dict)
        self.assertIsInstance(key_1, str)
        self.assertEqual(key_1, input_1[0].replace('-', ''))
        self.assertEqual(len(extent_1), 1)
        self.assertIsInstance(extent_1.values(), list)
        map(lambda el: self.assertIsInstance(el, float), *extent_1.values())
        self.assertEqual(extent_1, output_1)
        input_2 = ['-te', '5', 'str', '6', '61']
        try:
            key_2, extent_2 = Domain._add_to_dict(dict(), input_2)
        except OptionError as opt_err:
            self.assertEqual(opt_err.message, 'Input values must be int or float')

    def test_validate_te_lle(self):
        input_1 = [5., 60., 6., 61.]
        input_2 = ([5., 60., 5., 61.], [5., 60., 4., 61.], [5., 60., 6., 59.])
        input_3 = [60., 5., 61.]
        self.assertEqual(Domain._validate_te_lle(input_1), None)
        for inp in input_2:
            try:
                Domain._validate_te_lle(inp)
            except OptionError as opt_err:
                self.assertEqual(opt_err.message, 'Min cannot be bigger than max: '
                                                  '<-te x_min y_min x_max y_max> or '
                                                  '<-lle min_lon min_lat max_lon max_lat>')

        try:
            Domain._validate_te_lle(input_3)
        except OptionError as opt_err:
            self.assertEqual(opt_err.message, '-te and -lle requires exactly 4 parameters '
                                              '(3 given): <-te x_min y_min x_max y_max> or <-lle'
                                              ' min_lon min_lat max_lon max_lat>')

    def test_validate_ts_tr(self):
        input_1 = [100, 200]
        input_2 = ([0, 0], [0, 50], [10, 0], [-1, 10], [10, -100], [-100, -100])
        input_3 = [10]
        self.assertEqual(Domain._validate_ts_tr(input_1), None)

        for inp in input_2:
            try:
                Domain._validate_ts_tr(inp)
            except OptionError as opt_err:
                self.assertEqual(opt_err.message, 'Resolution or width and height must be bigger '
                                                  'than 0: <-tr x_resolution y_resolution> or '
                                                  '<-ts width height>')
        try:
            Domain._validate_ts_tr(input_3)
        except OptionError as opt_err:
            self.assertEqual(opt_err.message, '-ts and -tr requires exactly 2 parameters '
                                              '(1 given): <-tr x_resolution y_resolution> or '
                                              '<-ts width height>')

    def test_check_size(self):
        te_lle_example = '<-te x_min y_min x_max y_max> or <-lle min_lon min_lat max_lon max_lat>'
        tr_ts_example = '<-tr x_resolution y_resolution> or <-ts width height>'
        self.assertEqual(Domain._check_size(2, 2, ('-te', '-lle'), te_lle_example), None)

        try:
            Domain._check_size(1, 2, ('-ts', '-tr'), tr_ts_example)
        except OptionError as opt_err:
            self.assertEqual(opt_err.message, '-ts and -tr requires exactly 2 parameters '
                                              '(1 given): <-tr x_resolution y_resolution> or '
                                              '<-ts width height>')
        try:
            Domain._check_size(2, 4, ('-te', '-lle'), te_lle_example)
        except OptionError as opt_err:
            self.assertEqual(opt_err.message, '-te and -lle requires exactly 4 parameters '
                                              '(2 given): <-te x_min y_min x_max y_max> or <-lle'
                                              ' min_lon min_lat max_lon max_lat>')

    def test_gen_regexp(self):
        test_1 = '(-te|-lle)(\\s+[-+]?\\d*[.\\d*]*)(\\s+[-+]?\\d*[.\\d*]*)(\\s+[-+]?' \
                 '\\d*[.\\d*]*)(\\s+[-+]?\\d*[.\\d*]*)\\s?'
        result_1 = Domain._gen_regexp('te', 'lle', 4)
        self.assertIsInstance(result_1, str)
        self.assertEqual(result_1, test_1)
        test_2 = '(-ts|-tr)(\\s+[-+]?\\d*[.\\d*]*)(\\s+[-+]?\\d*[.\\d*]*)\\s?'
        result_2 = Domain._gen_regexp('ts', 'tr', 2)
        self.assertEqual(result_2, test_2)

    def test_create_extent_dict(self):
        test = ('-te 5 60 6 61.1 -ts 500 500',
                '-te -92.08 26.85 -92.00 26.91 -ts 200 200',
                '-te 5 60 6 61.1',
                '-te 5 60 6 61.1 -te 5 60 6 61.1')

        output_1 = {'te': [5., 60., 6., 61.1], 'ts': [500, 500]}
        output_2 = {'te': [-92.08, 26.85, -92.00, 26.91], 'ts': [200, 200]}
        result_1 = Domain._create_extent_dict(test[0])
        self.assertIsInstance(result_1, dict)
        self.assertEquals(len(result_1.keys()), 2)
        self.assertEquals(result_1, output_1)
        result_2 = Domain._create_extent_dict(test[1])
        self.assertEquals(result_2, output_2)

        try:
            test = Domain._create_extent_dict(test[2])
        except OptionError as opt_err:
            self.assertEquals(opt_err.message, '<extent_dict> must contains exactly 2 parameters '
                                               '("-te" or "-lle") and ("-ts" or "-tr")')

        try:
            test = Domain._create_extent_dict(test[3])
        except OptionError as opt_err:
            self.assertEquals(opt_err.message, '<extent_dict> must contains exactly 2 parameters '
                                               '("-te" or "-lle") and ("-ts" or "-tr")')

    def test_get_min_max_lat_lon(self):
        dom = Domain(4326, "-te 5 60 6 61 -ts 500 500")
        result = dom.get_min_max_lat_lon()
        self.assertIsInstance(result, tuple)
        self.assertLess(result[0], result[1])
        self.assertLess(result[2], result[3])
        self.assertEqual(dom.get_min_max_lat_lon(), (60.002, 61.0, 5.0, 5.998))

    def test_get_row_col_vector(self):
        test_1 = Domain._get_row_col_vector(250, 500)
        self.assertIsInstance(test_1, list)
        self.assertEquals(test_1, range(251))
        self.assertEquals(len(test_1), 251)
        test_2 = Domain._get_row_col_vector(500, 10)
        self.assertEquals(test_2, range(0, 550, 50))
        self.assertEquals(len(test_2), 10 + 1)

    def test_compound_row_col_vectors(self):

        result = Domain._compound_row_col_vectors(30, 40, range(0, 33, 3), range(0, 44, 4))
        output_col, output_row = result
        self.assertIsInstance(result, tuple)
        self.assertEquals(len(result), 2)
        test_col = [0,  3,  6,  9,  12,  15,  18,  21,  24,  27,  30,  30,  30,  30,  30,  30,
                    30,  30,  30,  30,  30,  30,  30,  27,  24,  21,  18,  15,  12,  9,  6,  3,
                    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]
        test_row = [0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  4,  8,  12,  16,  20,  24,
                    28,  32,  36,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,
                    40,  36,  32,  28,  24,  20,  16,  12,  8,  4,  0]

        self.assertEquals(output_col, test_col)
        self.assertEquals(output_row, test_row)

    def test_get_border(self):
        dom = Domain(4326, "-te 4.5 60 6 61 -ts 750 500")
        result = dom.get_border(nPoints=10)
        lat, lon = result
        self.assertEqual(type(lat), np.ndarray)
        self.assertEqual(type(lon), np.ndarray)
        self.assertIsInstance(result, tuple)
        self.assertEquals(len(result), 2)
        test_x = [4.5, 4.65, 4.8, 4.95, 5.1, 5.25, 5.4, 5.55, 5.7, 5.85, 6., 6., 6.,
                  6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 5.85, 5.7, 5.55, 5.4, 5.25,
                  5.1, 4.95, 4.8, 4.65, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5]
        test_y = [61., 61., 61., 61., 61., 61., 61., 61., 61., 61., 61., 61., 60.9, 60.8, 60.7, 60.6,
                  60.5, 60.4, 60.3, 60.2, 60.1, 60., 60., 60., 60., 60., 60., 60., 60., 60., 60.,
                  60., 60., 60., 60.1, 60.2, 60.3, 60.4, 60.5, 60.6, 60.7, 60.8, 60.9, 61.]

        self.assertEquals(list(lat), test_x)
        self.assertEquals(list(lon), test_y)

    def test_transform_ts(self):
        result = Domain._transform_ts(1.5, 1.0, [750.0, 500.0])
        self.assertIsInstance(result, tuple)
        self.assertEquals(len(result), 4)
        map(lambda el: self.assertIsInstance(el, float), result)

    def test_transform_tr(self):
        result = Domain._transform_ts(4.0, 1.3, [0.015, 0.005])
        self.assertIsInstance(result, tuple)
        self.assertEquals(len(result), 4)
        map(lambda el: self.assertIsInstance(el, float), result)

        try:
            result = Domain._transform_ts(4.0, 1.3, [5.0, 0.005])
        except OptionError as param_err:
            self.assertEqual(param_err.message,
                             '"-tr" is too large. width is 4.0, height is 1.3 ')

    def test_get_geotransform(self):
        input_1 = {'te': [27.0, 70.3, 31.0, 71.6], 'tr': [0.015, 0.005]}
        test_1 = ([27.0, 0.015, 0.0, 71.6, 0.0, -0.005], 266, 259)
        result = Domain._get_geotransform(input_1)
        self.assertIsInstance(result, tuple)
        self.assertEquals(len(result), 3)
        self.assertIsInstance(result[0], list)
        self.assertEquals(len(result[0]), 6)
        map(lambda el: self.assertIsInstance(el, float), result[0])
        self.assertIsInstance(result[1], int)
        self.assertIsInstance(result[2], int)

        self.assertEquals(result, test_1)

        input_2 = {'te': [4.5, 60.0, 6.0, 61.0], 'ts': [750.0, 500.0]}
        test_2 = ([4.5, 0.002, 0.0, 61.0, 0.0, -0.002], 750, 500)
        result = Domain._get_geotransform(input_2)
        self.assertEquals(result, test_2)

    def test_get_border_postgis(self):
        d = Domain(4326, '-te 25 70 35 72 -ts 500 500')
        result = d.get_border_postgis()
        self.assertIsInstance(result, str)
        self.assertEquals(result, "PolygonFromText('POLYGON((25.0 72.0,26.0 72.0,27.0 72.0,28.0 "
                                  "72.0,29.0 72.0,30.0 72.0,31.0 72.0,32.0 72.0,33.0 72.0,34.0 "
                                  "72.0,35.0 72.0,35.0 72.0,35.0 71.8,35.0 71.6,35.0 71.4,35.0 "
                                  "71.2,35.0 71.0,35.0 70.8,35.0 70.6,35.0 70.4,35.0 70.2,35.0 "
                                  "70.0,35.0 70.0,34.0 70.0,33.0 70.0,32.0 70.0,31.0 70.0,30.0 "
                                  "70.0,29.0 70.0,28.0 70.0,27.0 70.0,26.0 70.0,25.0 70.0,25.0 "
                                  "70.0,25.0 70.2,25.0 70.4,25.0 70.6,25.0 70.8,25.0 71.0,25.0 "
                                  "71.2,25.0 71.4,25.0 71.6,25.0 71.8,25.0 72.0))')")

    def test_repr(self):
        dom = Domain(4326, "-te 4.5 60 6 61 -ts 750 500")
        result = dom.__repr__()
        test = 'Domain:[750 x 500]\n----------------------------------------\nProjection:\nGEOGC' \
               'S["WGS 84",\n    DATUM["WGS_1984",\n        SPHEROID["WGS 84",6378137,298.257223' \
               '563]],\n    PRIMEM["Greenwich",0],\n    UNIT["degree",0.0174532925199433]]\n-----' \
               '-----------------------------------\nCorners (lon, lat):\n\t (  4.50,  61.00)  ' \
               '(  6.00,  61.00)\n\t (  4.50,  60.00)  (  6.00,  60.00)\n'

        self.assertIsInstance(result, str)
        self.assertEquals(result, test)

if __name__ == "__main__":
    unittest.main()
