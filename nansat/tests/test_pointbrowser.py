# ------------------------------------------------------------------------------
# Name:         test_pointbrowser.py
# Purpose:      Test the PointBrowser class
#
# Author:       Aleksander Vines, Anton Korosov
#
# Created:      2015-10-22
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------------------
import os
import unittest
from mock import patch, PropertyMock, Mock, MagicMock, DEFAULT

import numpy as np
from nansat.pointbrowser import PointBrowser

try:
    import matplotlib
    import matplotlib.pyplot as plt
    plt.switch_backend('Agg')
except ImportError:
    MATPLOTLIB_IS_INSTALLED = False
else:
    MATPLOTLIB_IS_INSTALLED = True


class PointBrowserTest(unittest.TestCase):
    @unittest.skipUnless(MATPLOTLIB_IS_INSTALLED, 'Matplotlib is required')
    def setUp(self):
        self.data = np.zeros((4, 4))
        self.event = MagicMock()

    def test_init(self):
        """ Create Pointbrowser """
        pb = PointBrowser(self.data, force_interactive=False)
        self.assertIsInstance(pb.fig, plt.Figure)
        self.assertTrue(np.alltrue(pb.data == self.data))
        self.assertTrue(np.alltrue(pb.ax.get_images()[0].get_array() == self.data))
        self.assertEqual(pb.fmt, 'x-k')
        self.assertEqual(pb.points, [])
        self.assertEqual(pb.coordinates, [[]])

    def test_onclick(self):
        """ Mimic click """
        self.event = MagicMock()
        self.event.xdata = 10
        self.event.ydata = 10
        self.event.key = None
        pb = PointBrowser(self.data, force_interactive=False)

        pb.onclick(self.event)
        self.assertIsInstance(pb.points[0][0], matplotlib.lines.Line2D)
        self.assertEqual(pb.coordinates, [[(self.event.xdata, self.event.ydata)]])

    def test_onclick_none(self):
        """ Mimic click outside figure  """
        self.event.xdata = None
        self.event.ydata = None
        self.event.key = None
        pb = PointBrowser(self.data, force_interactive=False)

        pb.onclick(self.event)
        self.assertEqual(pb.points, [])
        self.assertEqual(pb.coordinates, [[]])

    def test_onclick_key_z(self):
        """ Mimic click with 'z' pressed """
        self.event.xdata = 10
        self.event.ydata = 10
        self.event.key = 'z'
        pb = PointBrowser(self.data, force_interactive=False)

        pb.onclick(self.event)
        self.assertEqual(pb.points, [])
        self.assertEqual(pb.coordinates, [[]])

    def test_onclick_key(self):
        """ Mimic click with 'anykey' pressed """
        self.event = MagicMock()
        self.event.xdata = 10
        self.event.ydata = 10
        self.event.key = 'newkey'
        pb = PointBrowser(self.data, force_interactive=False)

        pb.onclick(self.event)
        self.assertIsInstance(pb.points[0][0], matplotlib.lines.Line2D)
        self.assertEqual(pb.coordinates, [[],[(self.event.xdata, self.event.ydata)]])

    def test_convert_coordinates(self):
        """ Mimic click with 'anykey' pressed """
        pb = PointBrowser(self.data, force_interactive=False)
        pb.coordinates = [[[1,2,3],[4,5,6]]]
        new_coordinates = pb._convert_coordinates()
        self.assertTrue(np.all(new_coordinates[0] == np.array([[1,2,3], [4,5,6]]).T))

    @patch('nansat.pointbrowser.plt')
    def test_get_points(self, plt_mock):
        plt_mock.show.return_value = None
        pb = PointBrowser(self.data, force_interactive=False)
        points = pb.get_points()
        self.assertTrue(pb.fig.canvas.mpl_connect.called)
        self.assertTrue(plt_mock.show.called)
        self.assertEqual(points, [])


if __name__ == "__main__":
    unittest.main()
