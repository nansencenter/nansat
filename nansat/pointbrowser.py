# Name:    pointbrowser.py
# Purpose: contains PointBrowser class
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov,
#               Aleksander Vines
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2015
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
import os
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


class PointBrowser():
    """
    Click on raster images shown by plt.imshow and get the X-Y coordinates.

    Parameters
    ----------
    data : ndarray
        image to imshow
    transect : bool
        if True, get transects / points
        if False, get only points
    **kwargs : dict
        optional parameters for imshow

    Creates
    -------
    self.fig        : pyplot Figure
    self.data       : ndarray with data
    self.ax         : axes
    self.points     : plot with points
    self.line       : plot with points
    self.coordinates: container for recorded coordinates

    """
    # instance attributes
    fig = None
    data = None
    fmt = None
    text_ax = None
    ax = None
    points = None
    lines = None
    coordinates = None

    def __init__(self, data, fmt='x-k', **kwargs):
        """Open figure with imshow and colorbar"""
        if not MATPLOTLIB_IS_INSTALLED:
            raise ImportError(' Matplotlib is not installed ')
        if not matplotlib.is_interactive():
            raise SystemError('''
        Python is started with -pylab option, transect will not work.
        Please restart python without -pylab.''')

        self.fig = plt.figure()
        self.data = data
        self.fmt = fmt
        self.text_ax = plt.axes([0.0, 0.85, 1.0, 0.15])
        self.ax = plt.axes([0.0, 0.0, 1.0, 0.85])
        img = self.ax.imshow(self.data, extent=(0, self.data.shape[1],
                                                0, self.data.shape[0]),
                             origin='lower', **kwargs)

        self.fig.colorbar(img)
        self.points = []
        self.lines = [self.ax.plot([], [], self.fmt)[0]]
        self.coordinates = [[]]

    def onclick(self, event):
        """Append onclick event"""
        # ignore click outside image
        if event.xdata is None or event.ydata is None:
            return

        # ignore clicked point if "z" key is held down
        if str(event.key) == 'alt+z' or str(event.key) == 'z':
            return

        if event.key is not None:
            # - holding down any other key (NOT cmd (mac),shift,alt,ctrl)
            #   means a new line is started at the clicked point
            self.coordinates.append([])
            self.lines.append(self.ax.plot([], [], self.fmt)[0])

        self.coordinates[-1].append((event.xdata, event.ydata))
        self.points.append(self.ax.plot(event.xdata, event.ydata, self.fmt))
        self.lines[-1].set_data(np.array(self.coordinates[-1]).T)
        self.ax.figure.canvas.draw()

    def _convert_coordinates(self):
        ''' Converts the coordinates array to points array for return in
        get_points, so that the internal structure can be tested

        The format of the returned array:
        [array([[x1,...,xn],[y1,...,yn]]),array([[xn+1,...],[yn+1,...]]),...]
        Each 'array' element is a numpy.ndarray and represents one transect,
        where x1,y1 is the first point in the first transect,
        and xn,yn the last point in the first transect.
        The inner x/y-arrays are also numpy.ndarrays
        '''
        return [np.array(p).T for p in self.coordinates if len(p) > 0]

    def get_points(self):
        ''' Enables the onclick events and returns the points.

        The format of the returned array:
        [array([[x1,...,xn],[y1,...,yn]]),array([[xn+1,...],[yn+1,...]]),...]
        Each 'array' element is a numpy.ndarray and represents one transect,
        where x1,y1 is the first point in the first transect,
        and xn,yn the last point in the first transect.
        The inner x/y-arrays are also numpy.ndarrays
        '''
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
                'To start drawing a new line: press "space",'
                ' click in the next\n'
                '        location, release "space" and continue clicking\n'
                'To zoom: press "z" and use pan/zoom tools, then release "z"')
        self.text_ax.text(0.01, 0.9, text, fontsize=13,
                          verticalalignment='top', horizontalalignment='left')

        # collect data
        plt.show()

        # convert list of lists of coordinates to list of arrays
        points = self._convert_coordinates()

        return points
