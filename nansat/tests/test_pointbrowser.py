# ------------------------------------------------------------------------------
# Name:         test_pointbrowser.py
# Purpose:      Test the PointBrowser class
#
# Author:       Aleksander Vines
#
# Created:      2015-10-22
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------------------
import unittest
import numpy as np
from nansat.pointbrowser import PointBrowser


class PointBrowserTest(unittest.TestCase):

    def test_onclick(self):
        data = np.ndarray(shape=(4, 4), dtype=float, order='F')
        point = PointBrowser(data)
        event = Event(xdata=0, ydata=0, key=None)

        point.onclick(event)
        t = point._convert_coordinates()[0]
        self.assertIsInstance(t, np.ndarray)
        xPoints = t[0]
        self.assertIsInstance(xPoints, np.ndarray)
        yPoints = t[1]
        self.assertIsInstance(yPoints, np.ndarray)
        self.assertEqual(xPoints[0], event.xdata, "x coordinates is set wrong")
        self.assertEqual(yPoints[0], event.ydata, "y coordinates is set wrong")

    def test_onclick_multilines(self):
        data = np.ndarray(shape=(4, 4), dtype=float, order='F')
        point = PointBrowser(data)
        events = []
        events.append(Event(xdata=0, ydata=0, key=None))
        events.append(Event(xdata=1, ydata=0, key=None))
        events.append(Event(xdata=2, ydata=2, key='AnyKeyButZorAltZ'))
        events.append(Event(xdata=2, ydata=3, key=None))
        for event in events:
            point.onclick(event)
        points = point._convert_coordinates()
        self.assertEqual(len(points), 2, 'There should be two transects')
        self.assertTrue(np.alltrue(points[0] == np.array([[0, 1], [0, 0]])),
                        't1 is not correct')
        self.assertTrue(np.alltrue(points[1] == np.array([[2, 2], [2, 3]])),
                        't2 is not correct')


class Event:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


if __name__ == "__main__":
    unittest.main()
