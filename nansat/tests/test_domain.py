#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the nansat module
#
# Author:       Morten Wergeland Hansen, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:  18.06.2014
# Last modified:27.08.2014 10:58
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import unittest, warnings
import os, sys, glob
from types import ModuleType, FloatType
import numpy as np

from nansat import Domain
from nansat.tools import OptionError, gdal, ogr
from nansat.figure import Image

tmp_data_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'data', 'test_data')

if not os.path.exists(tmp_data_path):
    os.mkdir(tmp_data_path)


class DomainTest(unittest.TestCase):
    def setUp(self):
        self.data_path = os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        'data')

        self.test_data = os.path.join(self.data_path, 'gcps.tif')

        if not os.path.exists(self.test_data):
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
        ds = gdal.Open(self.test_data)
        d = Domain(ds=ds)

        self.assertEqual(type(d), Domain)

    def test_dont_init_from_invalid(self):
        self.assertRaises(OptionError, Domain)
        self.assertRaises(OptionError, Domain, None)

    def test_write_kml(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(tmp_data_path, 'domain_write_kml.kml')
        d.write_kml(kmlFileName=tmpfilename)

        self.assertEqual(os.path.exists(tmpfilename), True)

    def test_get_geolocation_grids(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.get_geolocation_grids()

        self.assertEqual(type(lon), np.ndarray)
        self.assertEqual(type(lat), np.ndarray)
        self.assertEqual(lat.shape, (500,500))

    def test_get_border(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.get_border()

        self.assertEqual(type(lon), np.ndarray)
        self.assertEqual(type(lat), np.ndarray)
        self.assertEqual(len(lat), 44)
        self.assertEqual(lat[0], lat[-1])
        self.assertEqual(lon[0], lon[-1])

    def test_get_border_wkt(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        bwkt = d.get_border_wkt()

        self.assertEqual(type(bwkt), str)
        self.assertEqual('POLYGON' in bwkt, True)

    def test_get_border_geometry(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        geom = d.get_border_geometry()

        self.assertEqual(type(geom), ogr.Geometry)
        self.assertEqual(geom.IsValid(), True)

    def test_get_corners(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.get_corners()

        self.assertEqual(all(lon - [ 25.,  25.,  35.,  35.] < 0.01), True)
        self.assertEqual(all(lat - [ 72.,  70.,  72.,  70.] < 0.01), True)

    def test_get_pixelsize_meters(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        x, y = d.get_pixelsize_meters()

        self.assertEqual(x - 444 < 1, True)
        self.assertEqual(y - 723 < 1, True)

    def test_transform_points(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        lon, lat = d.transform_points([1,2,3], [1,2,3])

        self.assertEqual(type(lon), np.ndarray)
        self.assertEqual(type(lat), np.ndarray)

    def test_transform_points_inverse(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        x, y = d.transform_points([25, 26, 27], [70, 71, 72], 1)

        self.assertEqual(all(np.round(x) == [0, 50, 100]), True)
        self.assertEqual(all(np.round(y) == [500, 250, 0]), True)

    def test_upwards_azimuth_direction(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        uad = d.upwards_azimuth_direction()

        self.assertEqual(np.round(uad), 180)

    def test_azimuth_up(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        au = d.azimuth_up()

        self.assertEqual(np.round(au[0,0]), 0)
        self.assertEqual(np.round(au[10,10]), 0)

    def test_shape(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        shape = d.shape()

        self.assertEqual(shape, (500, 500))

    def test_write_map(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(tmp_data_path, 'domain_write_map.png')
        d.write_map(tmpfilename)

        self.assertEqual(os.path.exists(tmpfilename), True)
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (50, 50))

    def test_write_map_dpi100(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(tmp_data_path, 'domain_write_map_dpi100.png')
        d.write_map(tmpfilename, dpi=100)

        self.assertEqual(os.path.exists(tmpfilename), True)
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (100, 100))

    def tearDown(self):
        # if any test plots are created, they could be deleted here
        pass
        #tmpfiles = glob.glob(os.path.join(self.data_path, 'temp_domain*'))
        #for tmpfile in tmpfiles:
        #    os.remove(tmpfile)
