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


class NansatTest(unittest.TestCase):

    def test_onclick(self):
        data = np.ndarray(shape=(4, 4), dtype=float, order='F')
        point = PointBrowser(data)
        event = Event(xdata=0, ydata=0, key=None)

        point.onclick(event)
        x, y = point.coordinates[0]
        self.assertEqual(x, event.xdata, "x coordinates is set wrong")
        self.assertEqual(y, event.ydata, "y coordinates is set wrong")
        self.assertEqual(point.connect[0], 1, "connect is set wrong")

    def test_onclick_multilines(self):
        data = np.ndarray(shape=(4, 4), dtype=float, order='F')
        point = PointBrowser(data)
        events = []
        events.append(Event(xdata=0, ydata=0, key=None))
        events.append(Event(xdata=1, ydata=0, key=None))
        events.append(Event(xdata=2, ydata=2, key='AnyKeyButZorAltZ'))
        events.append(Event(xdata=2, ydata=3, key=None))
        connect = [1, 1, 0, 1]
        for event in events:
            point.onclick(event)
        for i in range(0, 4):
            x, y = point.coordinates[i]
            self.assertEqual(x, events[i].xdata,
                             "%d-x coordinates is set wrong" % i)
            self.assertEqual(y, events[i].ydata,
                             "%d-y coordinates is set wrong" % i)
            self.assertEqual(point.connect[i], connect[i],
                             "%d-connect is set wrong" % i)


class Event:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


if __name__ == "__main__":
    unittest.main()
