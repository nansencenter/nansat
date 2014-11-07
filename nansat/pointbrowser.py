# Name:    pointbrowser.py
# Purpose: contains PointBrowser class
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2013
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

import matplotlib.pyplot as plt


class PointBrowser():
    '''
    Click on raster images shown by plt.imshow and get the X-Y coordinates.

    '''
    def __init__(self, data, **kwargs):
        ''' Open figure with imshow and colorbar

        Parameters
        -----------
        data : ndarray
            image to imshow
        **kwargs : dict
            optional parameters for imshow

        Creates
        --------
        self.fig        : pyplot Figure
        self.data       : ndarray with data
        self.ax         : axes
        self.points     : plot with points
        self.line       : plot with points

        '''
        self.fig = plt.figure()
        self.data = data
        self.ax = self.fig.add_subplot(111)
        img = self.ax.imshow(self.data, extent=(0, self.data.shape[1],
                                                0, self.data.shape[0]),
                             origin='lower', **kwargs)
        self.fig.colorbar(img)
        self.points, = self.ax.plot([], [], '+', ms=12, color='b')
        self.line, = self.ax.plot([], [])
        self.coordinates = []

    def onclick(self, event):
        ''' Append onclick event '''
        if event.xdata is not None and event.ydata is not None:
            self.coordinates.append((event.xdata, event.ydata))
            tCoordinates = map(tuple, zip(*self.coordinates))
            self.points.set_data(tCoordinates)
            self.points.figure.canvas.draw()
            self.line.set_data(tCoordinates)
            self.line.figure.canvas.draw()

    def get_points(self):
        ''' Process click event '''
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.axes[0].set_xlim([0, self.data.shape[1]])
        self.fig.axes[0].set_ylim([0, self.data.shape[0]])
        text = ('1. Please click on the figure and mark a point or '
                'draw a line.\n2. Then close the figure.')
        plt.text(0, int(self.data.shape[0]*1.05), text, fontsize=13,
                 verticalalignment='top', horizontalalignment='left')
        plt.gca().invert_yaxis()
        plt.show()
