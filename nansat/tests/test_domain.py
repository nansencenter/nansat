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
import os
import numpy as np

try:
    if 'DISPLAY' not in os.environ:
        import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
except ImportError:
    BASEMAP_LIB_IS_INSTALLED = False
else:
    BASEMAP_LIB_IS_INSTALLED = True

from nansat.nsr import NSR
from nansat.vrt import VRT
from nansat.domain import Domain
from nansat.tools import gdal, ogr
from nansat.figure import Image
import sys
from nansat.tests import nansat_test_data as ntd
try:
    from mock import patch, Mock, MagicMock, PropertyMock, DEFAULT
except:
    from unittest.mock import patch, Mock, MagicMock, PropertyMock, DEFAULT

from nansat.exceptions import NansatProjectionError


EXTENT_TE_TS = "-te 25 70 35 72 -ts 500 500"
EXTENT_DICT_TE_TS = {'te': [25.0, 70.0, 35.0, 72.0], 'ts': [500.0, 500.0]}
EXTENT_LLE_TS = "-lle 25 70 35 72 -ts 500 500"
EXTENT_DICT_LLE_TS = {'lle': [25.0, 70.0, 35.0, 72.0], 'ts': [500.0, 500.0]}
EXTENT_DICT_LLE_TE_TS = {'lle': [25.0, 70.0, 35.0, 72.0],
                         'te': [25.0, 70.0, 35.0, 72.0], 'ts': [500.0, 500.0]}
GEO_TRANSFORM = [25.0, 0.02, 0.0, 72.0, 0.0, -0.004]
RASTER_X_SIZE = 500
RASTER_Y_SIZE = 500
SRS_PROJ4 = "+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs"
SRS_EPSG = 4326
EXTENT_BERGEN = "-te 5 60 6 61 -ts 500 500"
EXTENT_WESTCOAST = "-te 1 58 6 64 -ts 500 500"
EXTENT_NORWAY = "-te 3 55 30 72 -ts 500 500"
EXTENT_PARIS = "-te 2 48 3 49 -ts 500 500"


class DomainTest(unittest.TestCase):

    def setUp(self):
        self.test_file = os.path.join(ntd.test_data_path, 'gcps.tif')
        self.test_file_projected = os.path.join(ntd.test_data_path, 'stere.tif')
        if BASEMAP_LIB_IS_INSTALLED:
            plt.switch_backend('Agg')
        if (    not os.path.exists(self.test_file)
             or not os.path.exists(self.test_file_projected) ):
            raise ValueError('No test data available')

    def test_dont_init_from_invalid_combination(self):
        self.assertRaises(ValueError, Domain)
        self.assertRaises(ValueError, Domain, None)
        with self.assertRaises(ValueError):
            Domain(ds=gdal.Open(self.test_file),
                   srs="+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs",
                   ext="-te 25 70 35 72 -ts 500 500")

    def test_init_from_GDALDataset(self):
        d = Domain(ds=gdal.Open(self.test_file))
        self.assertEqual(type(d), Domain)

    def test_init_from_GDALDataset_and_srs(self):
        d = Domain(ds=gdal.Open(self.test_file),
                   srs="+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs")
        self.assertEqual(type(d), Domain)

    @patch('nansat.domain.gdal')
    def test_dont_init_if_gdal_AutoCreateWarpedVRT_fails(self, mock_gdal):
        mock_gdal.AutoCreateWarpedVRT.return_value = None
        with self.assertRaises(NansatProjectionError):
            Domain(ds=gdal.Open(self.test_file),
                   srs="+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs")

    @patch.object(Domain, '_create_extent_dict',
        return_value={'te': [25.0, 70.0, 35.0, 72.0], 'ts': [500.0, 500.0]})
    @patch.object(Domain, '_get_geotransform',
        return_value=([25.0, 0.02, 0.0, 72.0, 0.0, -0.004], 500, 500))
    def test_init_from_srs_and_ext_te(self, mock__get_geotransform, mock__create_extent_dict):
        d = Domain(srs=4326,
                   ext="-te 25 70 35 72 -ts 500 500")
        self.assertEqual(type(d), Domain)

    @patch.object(Domain, '_create_extent_dict',
        return_value={'lle': [25.0, 70.0, 35.0, 72.0], 'ts': [500.0, 500.0]})
    @patch.object(Domain, '_convert_extentDic',
        return_value={'lle': [25.0, 70.0, 35.0, 72.0], 'te': [25.0, 70.0, 35.0, 72.0], 'ts': [500.0, 500.0]})
    @patch.object(Domain, '_get_geotransform',
        return_value=([25.0, 0.02, 0.0, 72.0, 0.0, -0.004], 500, 500))
    def test_init_from_srs_and_ext_lle(self, mock__get_geotransform, mock__convert_extentDic,
                                       mock__create_extent_dict):
        d = Domain(srs=4326,
                   ext="-lle 25 70 35 72 -ts 500 500")
        self.assertEqual(type(d), Domain)

    def test_init_from_lonlat(self):
        lat, lon = np.mgrid[-90:90:0.5, -180:180:0.5]
        d = Domain(lon=lon, lat=lat)
        self.assertEqual(type(d), Domain)
        self.assertEqual(d.shape(), lat.shape)

    @patch.object(Domain, 'get_corners',
        return_value=(np.array([ 25.,  25.,  35.,  35.]), np.array([ 72.,  70.,  72.,  70.])))
    def test_repr(self, mock_get_corners):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        result = d.__repr__()
        test = ('Domain:[500 x 500]\n'
                '----------------------------------------\n'
                'Projection:\nGEOGCS["WGS 84",\n'
                '    DATUM["WGS_1984",\n'
                '        SPHEROID["WGS 84",6378137,298.257223563]],\n'
                '    PRIMEM["Greenwich",0],\n'
                '    UNIT["degree",0.0174532925199433]]\n'
                '----------------------------------------\n'
                'Corners (lon, lat):\n'
                '\t ( 25.00,  72.00)  ( 35.00,  72.00)\n'
                '\t ( 25.00,  70.00)  ( 35.00,  70.00)\n' )
        self.assertIsInstance(result, str)
        self.assertEquals(result, test)

    def test_write_kml(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(ntd.tmp_data_path, 'domain_write_kml.kml')
        d.write_kml(kmlFileName=tmpfilename)
        self.assertTrue(os.path.exists(tmpfilename))

    #def test__get_border_kml(self):


    #def test_write_kml_image(self):


    @patch.object(Domain, 'transform_points',
        return_value=(np.meshgrid(range(0,500),range(0,500))[0].flatten()*(35-25)/500.+25,
                      np.meshgrid(range(0,500),range(0,500))[1].flatten()*(70-72)/500.+72))
    def test_get_geolocation_grids_from_GDAL_transformer(self, mock_transform_points):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.get_geolocation_grids()
        self.assertEqual(type(lon), np.ndarray)
        self.assertEqual(type(lat), np.ndarray)
        self.assertEqual(lat.shape, (500, 500))

    def test_get_geolocation_grids_from_geolocationArray(self):
        lat, lon = np.mgrid[25:35:0.02, 70:72:0.004]
        d = Domain(lon=lon, lat=lat)
        lon, lat = d.get_geolocation_grids()
        self.assertEqual(type(lon), np.ndarray)
        self.assertEqual(type(lat), np.ndarray)
        self.assertEqual(lat.shape, (500, 500))




    def test_convert_extentDic(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        result = d._convert_extentDic(NSR(4326), {'lle': [25.0, 70.0, 35.0, 72.0], 'ts': [500.0, 500.0]})
        self.assertEqual(result, {'lle': [25.0, 70.0, 35.0, 72.0], 'te': [25.0, 70.0, 35.0, 72.0], 'ts': [500.0, 500.0]})

    def test_add_to_dict(self):
        input_1 = ['-te', '5', '60', '6', '61.1']
        output_1 = {'te': [5., 60., 6., 61.1]}
        key_1, extent_1 = Domain._add_to_dict(dict(), input_1)
        self.assertIsInstance(extent_1, dict)
        self.assertIsInstance(key_1, str)
        self.assertEqual(key_1, input_1[0].replace('-', ''))
        self.assertEqual(len(extent_1), 1)
        self.assertIsInstance(list(extent_1.values()), list)

        for el in list(extent_1.values())[0]:
            self.assertIsInstance(el, float)

        self.assertEqual(extent_1, output_1)
        input_2 = ['-te', '5', 'str', '6', '61']
        with self.assertRaises(ValueError) as opt_err:
            key_2, extent_2 = Domain._add_to_dict(dict(), input_2)
            self.assertEqual(opt_err.args[0], 'Input values must be int or float')

    def test_validate_ts_tr(self):
        input_1 = [100, 200]
        input_2 = ([0, 0], [0, 50], [10, 0], [-1, 10], [10, -100], [-100, -100])
        input_3 = [10]
        self.assertEqual(Domain._validate_ts_tr(input_1), None)
        for inp in input_2:
            with self.assertRaises(ValueError) as opt_err:
                Domain._validate_ts_tr(inp)
                self.assertEqual(opt_err.args[0], 'Resolution or width and height must be bigger '
                                                  'than 0: <-tr x_resolution y_resolution> or '
                                                  '<-ts width height>')
        with self.assertRaises(ValueError) as opt_err:
            Domain._validate_ts_tr(input_3)
            self.assertEqual(opt_err.args[0], '-ts and -tr requires exactly 2 parameters '
                                              '(1 given): <-tr x_resolution y_resolution> or '
                                              '<-ts width height>')
    def test_validate_te_lle(self):
        input_1 = [5., 60., 6., 61.]
        input_2 = ([5., 60., 5., 61.], [5., 60., 4., 61.], [5., 60., 6., 59.])
        input_3 = [60., 5., 61.]
        self.assertEqual(Domain._validate_te_lle(input_1), None)
        for inp in input_2:
            with self.assertRaises(ValueError) as opt_err:
                Domain._validate_te_lle(inp)
                self.assertEqual(opt_err.message, 'Min cannot be bigger than max: '
                                                  '<-te x_min y_min x_max y_max> or '
                                                  '<-lle min_lon min_lat max_lon max_lat>')
        with self.assertRaises(ValueError) as opt_err:
            Domain._validate_te_lle(input_3)
            self.assertEqual(opt_err.args[0], '-te and -lle requires exactly 4 parameters '
                                              '(3 given): <-te x_min y_min x_max y_max> or <-lle'
                                              ' min_lon min_lat max_lon max_lat>')

    def test_check_size(self):
        te_lle_example = '<-te x_min y_min x_max y_max> or <-lle min_lon min_lat max_lon max_lat>'
        tr_ts_example = '<-tr x_resolution y_resolution> or <-ts width height>'
        self.assertEqual(Domain._check_size(2, 2, ('-te', '-lle'), te_lle_example), None)
        with self.assertRaises(ValueError) as opt_err:
            Domain._check_size(1, 2, ('-ts', '-tr'), tr_ts_example)
            self.assertEqual(opt_err.args[0], '-ts and -tr requires exactly 2 parameters '
                                              '(1 given): <-tr x_resolution y_resolution> or '
                                              '<-ts width height>')
        with self.assertRaises(ValueError) as opt_err:
            Domain._check_size(2, 4, ('-te', '-lle'), te_lle_example)
            self.assertEqual(opt_err.args[0], '-te and -lle requires exactly 4 parameters '
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
        self.assertEqual(len(list(result_1.keys())), 2)
        self.assertEqual(result_1, output_1)
        result_2 = Domain._create_extent_dict(test[1])
        self.assertEquals(result_2, output_2)
        with self.assertRaises(ValueError) as opt_err:
            Domain._create_extent_dict(test[2])
            self.assertEquals(opt_err.args[0], '<extent_dict> must contains exactly 2 parameters ')
        self.assertEqual(result_2, output_2)

        try:
            test = Domain._create_extent_dict(test[2])
        except ValueError as opt_err:
            self.assertEqual(opt_err.args[0], '<extent_dict> must contains exactly 2 parameters '
                                               '("-te" or "-lle") and ("-ts" or "-tr")')
        with self.assertRaises(ValueError) as opt_err:
            Domain._create_extent_dict(test[3])
            self.assertEquals(opt_err.args[0], '<extent_dict> must contains exactly 2 parameters '
                                               '("-te" or "-lle") and ("-ts" or "-tr")')

    def test_get_border(self):
        dom = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        result = dom.get_border(nPoints=10)
        lat, lon = result
        self.assertEqual(type(lat), np.ndarray)
        self.assertEqual(type(lon), np.ndarray)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        test_x = [25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35.,
                  35., 35., 35., 35., 35., 35., 35., 35., 35., 35., 35.,
                  35., 34., 33., 32., 31., 30., 29., 28., 27., 26., 25.,
                  25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25.]
        test_y = [72., 72., 72., 72., 72., 72., 72., 72.,  72.,
                   72. ,  72. ,  72. ,  71.8,  71.6,  71.4,  71.2,  71. ,  70.8,
                   70.6,  70.4,  70.2,  70. ,  70. ,  70. ,  70. ,  70. ,  70. ,
                   70. ,  70. ,  70. ,  70. ,  70. ,  70. ,  70. ,  70.2,  70.4,
                   70.6,  70.8,  71. ,  71.2,  71.4,  71.6,  71.8,  72.]
        self.assertEqual(list(lat), test_x)
        self.assertEqual(list(lon), test_y)

    def test_compound_row_col_vectors(self):
        result = Domain._compound_row_col_vectors(30, 40, list(range(0, 33, 3)), list(range(0, 44, 4)))
        output_col, output_row = result
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        test_col = [0,  3,  6,  9,  12,  15,  18,  21,  24,  27,  30,  30,  30,  30,  30,  30,
                    30,  30,  30,  30,  30,  30,  30,  27,  24,  21,  18,  15,  12,  9,  6,  3,
                    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]
        test_row = [0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  4,  8,  12,  16,  20,  24,
                    28,  32,  36,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,  40,
                    40,  36,  32,  28,  24,  20,  16,  12,  8,  4,  0]

        self.assertEqual(output_col, test_col)
        self.assertEqual(output_row, test_row)

    def test_get_row_col_vector(self):
        test_1 = Domain._get_row_col_vector(250, 500)
        self.assertIsInstance(test_1, list)
        self.assertEqual(test_1, list(range(251)))
        self.assertEqual(len(test_1), 251)
        test_2 = Domain._get_row_col_vector(500, 10)
        self.assertEqual(test_2, list(range(0, 550, 50)))
        self.assertEqual(len(test_2), 10 + 1)

    def test_get_border_wkt(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        bwkt = d.get_border_wkt()
        self.assertEqual(type(bwkt), str)
        self.assertTrue('POLYGON' in bwkt)

    def test_get_border_geometry(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        geom = d.get_border_geometry()
        self.assertEqual(type(geom), ogr.Geometry)

    def test_overlaps_intersects_and_contains(self):
        Bergen = Domain(4326, "-te 5 60 6 61 -ts 500 500")
        WestCoast = Domain(4326, "-te 1 58 6 64 -ts 500 500")
        Norway = Domain(4326, "-te 3 55 30 72 -ts 500 500")
        Paris = Domain(4326, "-te 2 48 3 49 -ts 500 500")

        self.assertFalse(Bergen.overlaps(Norway))
        self.assertFalse(Norway.overlaps(Bergen))
        self.assertTrue(Norway.overlaps(WestCoast))
        self.assertFalse(Paris.overlaps(Norway))

        self.assertTrue(Bergen.intersects(Norway))
        self.assertTrue(Norway.intersects(Bergen))
        self.assertTrue(Norway.intersects(WestCoast))
        self.assertFalse(Paris.intersects(Norway))

        self.assertTrue(Norway.contains(Bergen))
        self.assertFalse(Bergen.contains(Norway))
        self.assertFalse(Norway.contains(WestCoast)) # why false?
        self.assertFalse(Paris.contains(Norway))


    def test_contains(self):
        Bergen = Domain(4326, EXTENT_BERGEN)
        WestCoast = Domain(4326, EXTENT_WESTCOAST)
        Norway = Domain(4326, EXTENT_NORWAY)
        Paris = Domain(4326, EXTENT_PARIS)
        self.assertTrue(Norway.contains(Bergen))
        self.assertFalse(Bergen.contains(Norway))
        self.assertFalse(Norway.contains(WestCoast))
        self.assertFalse(Paris.contains(Norway))

    def test_transform_ts(self):
        result = Domain._transform_ts(1.5, 1.0, [750.0, 500.0])
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 4)
        for el in result:
            self.assertIsInstance(el, float)
    def test_get_border_postgis(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
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

    def test_get_corners(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.get_corners()
        self.assertTrue(all(lon - [25., 25., 35., 35.] < 0.01))
        self.assertTrue(all(lat - [72., 70., 72., 70.] < 0.01))

    def test_get_min_max_lon_lat(self):
        dom = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        result = dom.get_min_max_lon_lat()
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 4)
        for el in result:
            self.assertIsInstance(el, float)
        self.assertLess(result[0], result[1])
        self.assertLess(result[2], result[3])
        self.assertEqual(result, (25.0, 34.980000000000004, 70.004000000000005, 72.0))

    def test_get_pixelsize_meters(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        x, y = d.get_pixelsize_meters()
        self.assertEqual(int(x), 444)
        self.assertEqual(int(y), 723)
        d = Domain(ds=gdal.Open(self.test_file_projected))
        x, y = d.get_pixelsize_meters()
        self.assertEqual(int(x), 500)
        self.assertEqual(int(y), 500)
    def test_get_geotransform(self):
        input_1 = {'te': [25.0, 70.0, 35.0, 72.0], 'ts': [500.0, 500.0]}
        test_1 = ([25.0, 0.02, 0.0, 72.0, 0.0, -0.004], 500, 500)
        result = Domain._get_geotransform(input_1)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 3)
        self.assertIsInstance(result[0], list)
        self.assertEqual(len(result[0]), 6)
        for el in result[0]:
            self.assertIsInstance(el, float)
        self.assertIsInstance(result[1], int)
        self.assertIsInstance(result[2], int)
        self.assertEqual(result, test_1)

    def test_transform_tr(self):
        result = Domain._transform_tr(4.0, 1.3, [0.015, 0.005])
        self.assertIsInstance(result, tuple)
        self.assertEquals(len(result), 4)
        map(lambda el: self.assertIsInstance(el, float), result)
        with self.assertRaises(ValueError) as param_err:
            Domain._transform_tr(4.0, 1.3, [5.0, 0.005])
            self.assertEqual(param_err.message,
                             '"-tr" is too large. width is 4.0, height is 1.3 ')

    def test_transform_ts(self):
        result = Domain._transform_ts(1.5, 1.0, [750.0, 500.0])
        self.assertIsInstance(result, tuple)
        self.assertEquals(len(result), 4)
        map(lambda el: self.assertIsInstance(el, float), result)

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

    @patch.object(Domain, 'get_geolocation_grids',
                  return_value=(np.array([[-4.        , -0.66666667,  2.66666667],
                                          [-4.        , -0.66666667,  2.66666667]]),
                                np.array([[ 7.,  7.,  7.],
                                          [ 1.,  1.,  1.]])))
    def test_azimuth_y(self, mock_get_geolocation_grids):
        d = Domain(4326, "-te -4 -5 +6 +7 -ts 3 2")
        self.assertTrue((d.azimuth_y()==np.array([[ 0.,  0.,  0.], [ 0.,  0.,  0.]])).all())

    def test_shape(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        self.assertEqual(d.shape(), (500, 500))

    @unittest.skipUnless(BASEMAP_LIB_IS_INSTALLED, 'Basemap is required')
    def test_write_map(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(ntd.tmp_data_path, 'domain_write_map.png')
        d.write_map(tmpfilename)
        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (50, 50))

    @unittest.skipUnless(BASEMAP_LIB_IS_INSTALLED, 'Basemap is required')
    def test_write_map_dpi100(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'domain_write_map_dpi100.png')
        d.write_map(tmpfilename, dpi=100)
        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (100, 100))

    @unittest.skipUnless(BASEMAP_LIB_IS_INSTALLED, 'Basemap is required')
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

    def test_reproject_gcps(self):
        ds = gdal.Open(self.test_file)
        d = Domain(ds=ds)
        d.reproject_gcps('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=75 +lon_0=10 +no_defs')
        gcp = d.vrt.dataset.GetGCPs()[0]
        self.assertTrue(gcp.GCPX > 636161)
        self.assertTrue(gcp.GCPY < -288344)

    def test_reproject_gcps_auto(self):
        ds = gdal.Open(self.test_file)
        d = Domain(ds=ds)
        d.reproject_gcps()
        gcpproj = NSR(d.vrt.dataset.GetGCPProjection())
        self.assertEqual(gcpproj.GetAttrValue('PROJECTION'),
                        'Stereographic')

    def test_reproject_GCPs(self):
        ds = gdal.Open(self.test_file)
        d = Domain(ds=ds)
        d.reproject_GCPs('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=75 +lon_0=10 +no_defs')
        gcp = d.vrt.dataset.GetGCPs()[0]
        self.assertTrue(gcp.GCPX > 636161)
        self.assertTrue(gcp.GCPY < -288344)


    def test_repr(self):
        dom = Domain(4326, "-te 4.5 60 6 61 -ts 750 500")
        result = dom.__repr__()
        test = 'Domain:[750 x 500]\n----------------------------------------\nProjection:\nGEOGC' \
               'S["WGS 84",\n    DATUM["WGS_1984",\n        SPHEROID["WGS 84",6378137,298.257223' \
               '563]],\n    PRIMEM["Greenwich",0],\n    UNIT["degree",0.0174532925199433]]\n-----' \
               '-----------------------------------\nCorners (lon, lat):\n\t (  4.50,  61.00)  ' \
               '(  6.00,  61.00)\n\t (  4.50,  60.00)  (  6.00,  60.00)\n'

        self.assertIsInstance(result, str)
        self.assertEqual(result, test)

    @patch.multiple(Domain, get_border_geometry=DEFAULT, __init__ = Mock(return_value=None))
    def test_intersects(self, get_border_geometry):
        other_domain = MagicMock()
        d = Domain()
        d.intersects(other_domain)
        d.get_border_geometry.assert_called_once()
        d.get_border_geometry().Intersects.assert_called_once_with(other_domain.get_border_geometry())

    @patch.multiple(Domain, get_border_geometry=DEFAULT, __init__ = Mock(return_value=None))
    def test_overlaps(self, get_border_geometry):
        other_domain = MagicMock()
        d = Domain()
        d.overlaps(other_domain)
        d.get_border_geometry.assert_called_once()
        d.get_border_geometry().Overlaps.assert_called_once_with(other_domain.get_border_geometry())

if __name__ == "__main__":
    unittest.main()
