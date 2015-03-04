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

import matplotlib
import matplotlib.pyplot as plt


class PointBrowser():
    '''
    Click on raster images shown by plt.imshow and get the X-Y coordinates.

    '''
    def __init__(self, data, transect=True, **kwargs):
        ''' Open figure with imshow and colorbar

        Parameters
        -----------
        data : ndarray
            image to imshow
        transect : bool
            if True, get transects / points
            if False, get only points
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
        self.text_ax = plt.axes([0.0, 0.85, 1.0, 0.15])
        self.ax = plt.axes([0.0, 0.0, 1.0, 0.85])
        img = self.ax.imshow(self.data, extent=(0, self.data.shape[1],
                                                0, self.data.shape[0]),
                             origin='lower', **kwargs)

        self.fig.colorbar(img)
        self.points, = self.ax.plot([], [], '+', ms=12, color='b')
        self.lines = []
        self.coordinates = []
        self.connect = []
        self.drawLine = transect

    def onclick(self, event):
        ''' Append onclick event '''
        if event.xdata is not None and event.ydata is not None:
            if str(event.key)=='alt+z' or str(event.key)=='z':
                pass
            else:
                # ignore clicked point if "z" key is held down
                # - holding down any other key (NOT cmd (mac),shift,alt,ctrl)
                #   means a new line is started at the clicked point
                # - holding down no key means current line is extended to include the clicked point
                self.coordinates.append((event.xdata, event.ydata))
                # press (any) key (NOT 'cmd','ctrl','alt','shift', or 'z' - see above) means to start new line.
                # if pressed, then set 0 to self.connect. otherwise set 1.
                if event.key is None and self.drawLine:
                   self.connect.append(1)
                else:
                   self.connect.append(0)

                # get coordinate of clicked point
                tCoordinates = map(tuple, zip(*self.coordinates))
                self.points.set_data(tCoordinates)
                self.points.figure.canvas.draw()

                # separate points by each line
                linesCoords = []
                for i, iLine in enumerate(self.coordinates):
                   if i == 0:
                       oneLine = [self.coordinates[0]]
                   elif self.connect[i] == 0:
                       linesCoords.append(oneLine)
                       oneLine = [self.coordinates[i]]
                   else:
                       oneLine.append(self.coordinates[i])
                linesCoords.append(oneLine)

                # draw lines
                if self.drawLine:
                    line, = self.ax.plot([], [])
                    for iLinePoints in linesCoords:
                       tCoordinates = map(tuple, zip(*iLinePoints))
                       self.lines.append(line.set_data(tCoordinates))
                       line.figure.canvas.draw()

    def get_points(self):
        ''' Process click event '''
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.ax.set_xlim([0, self.data.shape[1]])
        self.ax.set_ylim([0, self.data.shape[0]])
        self.ax.invert_yaxis()
        self.ax.tick_params(direction='in', pad=-20, labelsize=10)

        self.text_ax.set_xlim([0, 1])
        self.text_ax.set_ylim([0, 1])
        self.text_ax.set_xticks([])
        self.text_ax.set_yticks([])

        text = ('To draw a transect line: click in several locations\n'
                'To start drawing a new line: press "space", click in the next\n'
                '        location, release "space" and continue clicking\n'
                'To zoom: press "z" and use pan/zoom tools, then release "z"')
        self.text_ax.text(0.01, 0.9, text, fontsize=13,
                 verticalalignment='top', horizontalalignment='left')

        plt.show()
