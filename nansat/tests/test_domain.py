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
import matplotlib.pyplot as plt

from nansat.nsr import NSR
from nansat.domain import Domain
from nansat.tools import OptionError, gdal, ogr
from nansat.figure import Image

import nansat_test_data as ntd


class DomainTest(unittest.TestCase):
    def setUp(self):
        self.test_file = os.path.join(ntd.test_data_path, 'gcps.tif')
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

    def test_write_map(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(ntd.tmp_data_path, 'domain_write_map.png')
        d.write_map(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (50, 50))

    def test_write_map_dpi100(self):
        d = Domain(4326, "-te 25 70 35 72 -ts 500 500")
        tmpfilename = os.path.join(ntd.tmp_data_path,
                                   'domain_write_map_dpi100.png')
        d.write_map(tmpfilename, dpi=100)

        self.assertTrue(os.path.exists(tmpfilename))
        i = Image.open(tmpfilename)
        i.verify()
        self.assertEqual(i.info['dpi'], (100, 100))

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


if __name__ == "__main__":
    unittest.main()
