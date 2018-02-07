# Name:    nansat_map.py
# Purpose: Container of Nansatap class
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
from __future__ import absolute_import

import os
import re
import warnings

import numpy as np

try:
    if 'DISPLAY' not in os.environ:
        import matplotlib; matplotlib.use('Agg')
    from mpl_toolkits.basemap import Basemap
    import matplotlib as mpl
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from scipy import ndimage
except ImportError:
    BASEMAP_LIB_IS_INSTALLED = False
    MATPLOTLIB_IS_INSTALLED = False
else:
    BASEMAP_LIB_IS_INSTALLED = True
    MATPLOTLIB_IS_INSTALLED = True


from nansat.nsr import NSR
from nansat.tools import get_random_color
from nansat.warnings import NansatFutureWarning

NANSATMAP_WARNING = ('Nansatmap() will be disabled in Nansat 1.1. and moved to a separate package'
                     'It is not covered by unittests intentionally.')
warnings.warn(NANSATMAP_WARNING, NansatFutureWarning)

class Nansatmap(Basemap):
    '''Perform opeartions with graphical files: create,
    add legend and geolocation_grids, save.

    NansatMap instance is created in the Nansat.write_map method.
    The methods below are applied consequently in order to get projection,
    generate a basemap from array(s), add legend and geolocation grids,
    save to a file.

    '''
    # general attributes
    cmap = cm.jet
    colorbar = None
    mpl = []
    lon, lat, x, y = None, None, None, None
    # parameters for smoothing
    # convolve
    convolve_weightSize = 7
    convolve_weights = None
    convolve_mode = 'reflect'
    convolve_cval = 0.0
    convolve_origin = 0
    # fourier_gaussian
    fourier_sigma = 1.0
    fourier_n = -1
    fourier_axis = -1
    # spline
    spline_order = 3
    spline_axis = -1
    # gaussian filter
    gaussian_sigma = 2.5
    gaussian_order = 0
    gaussian_mode = 'reflect'
    gaussian_cval = 0.0
    # saving parameters
    DEFAULT_EXTENSION = '.png'


    def __init__(self, domain, **kwargs):
        ''' Set attributes
        Get proj4 from the given domain and convert the proj4 projection to
        the basemap projection.

        Parameters
        -----------
        domain : domain object
        kwargs : dictionary
            parameters that are used for all operations.

        Modifies
        ---------
        self.fig : figure
            matplotlib.pyplot.figure
        self.colorbar : boolean
            if colorbar is True, it is possible to put colorbar.
            e.g. contour_plots(contour_style='fill'), put_color()
        self.mpl : list
            elements are matplotlib.contour.QuadContourSet instance,
                         matplotlib.quiver.Quiver instance or
                         matplotlib.collections.QuadMesh object

        See also
        ----------
        http://matplotlib.org/basemap/api/basemap_api.html

        '''
        if not MATPLOTLIB_IS_INSTALLED:
            raise ImportError('Matplotlib is not installed')
        if not BASEMAP_LIB_IS_INSTALLED:
            raise ImportError('Basemap is not installed')
        warnings.warn(NANSATMAP_WARNING, NansatFutureWarning)

        self.domain = domain

        # get proj4
        proj4 = NSR(domain.vrt.get_projection()).ExportToProj4()

        # convert proj4 to basemap projection
        projStr = proj4.split(' ')[0][6:]
        projection = {'aea': 'aea', 'ocea': 'aea',
                      'aeqd': 'aeqd', 'xxx1': 'spaeqd', 'xxx2': 'npaeqd',
                      'cass': 'cass',
                      'cea': 'cea',
                      'eqc': 'cyl', 'longlat': 'cyl',
                      'eck4': 'eck4',
                      'eqdc': 'eqdc',
                      'gall': 'gall',
                      'geos': 'geos',
                      'gnom': 'gnom',
                      'hammer': 'hammer', 'nell_h': 'hammer',
                      'kav7': 'kav7',
                      'laea': 'laea', 'xxx3': 'splaea', 'xxx4': 'nplaea',
                      'lcc': 'lcc', 'lcca': 'lcc',
                      'mbtfpq': 'mbtfpq',
                      'somerc': 'merc', 'merc': 'merc', 'omerc': 'merc',
                      'mill': 'mill',
                      'moll': 'moll',
                      'nsper': 'nsper',
                      'omerc': 'omerc',
                      'ortho': 'ortho',
                      'poly': 'poly', 'rpoly': 'poly', 'imw_p': 'poly',
                      'robin': 'robin',
                      'sinu': 'sinu', 'fouc_s': 'sinu', 'gn_sinu': 'sinu',
                      'mbtfps': 'sinu', 'urmfps': 'sinu',
                      'stere': 'stere', 'sterea': 'stere', 'lee_os': 'stere',
                      'mil_os': 'stere', 'rouss': 'stere',
                      'ups': 'npstere', 'ups': 'spstere',  # CHECK!!
                      'tmerc': 'tmerc', 'gstmerc': 'tmerc', 'utm': 'tmerc',
                      'vandg': 'vandg', 'vandg2': 'vandg',
                      'vandg3': 'vandg', 'vandg4': 'vandg',
                      }.get(projStr, 'cyl')

        if projection in ['stere']:
            lon_0 = float(re.findall('lon_0=+[-+]?\d*[.\d*]*',
                                     proj4)[0].split('=')[1])
            lat_0 = float(re.findall('lat_0=+[-+]?\d*[.\d*]*',
                                     proj4)[0].split('=')[1])
            kwargs['lon_0'] = lon_0
            kwargs['lat_0'] = lat_0

        if projStr == 'utm':
            kwargs['lon_0'] = -180 + NSR(proj4).GetUTMZone()*6 - 3
            kwargs['lat_0'] = 0

        self.extensionList = ['png', 'emf', 'eps', 'pdf', 'rgba',
                              'ps', 'raw', 'svg', 'svgz']

        # set llcrnrlat, urcrnrlat, llcrnrlon and urcrnrlon to kwargs.
        # if required, modify them from -90. to 90.
        # get min/max lon/lat
        lonCrn, latCrn = domain.get_corners()
        self.lonMin = min(lonCrn)
        self.lonMax = max(lonCrn)
        self.latMin = max(min(latCrn), -90.)
        self.latMax = min(max(latCrn), 90.)

        if not('llcrnrlat' in kwargs.keys()):
            kwargs['llcrnrlat'] = latCrn[1]
        if not('urcrnrlat' in kwargs.keys()):
            kwargs['urcrnrlat'] = latCrn[2]
        if not('llcrnrlon' in kwargs.keys()):
            kwargs['llcrnrlon'] = lonCrn[1]
        if not('urcrnrlon' in kwargs.keys()):
            kwargs['urcrnrlon'] = lonCrn[2]

        # separate kwarge of plt.figure() from kwargs
        figArgs = ['num', 'figsize', 'dpi', 'facecolor', 'edgecolor',
                   'frameon']
        figKwargs = {}
        for iArg in figArgs:
            if iArg in kwargs.keys():
                figKwargs[iArg] = kwargs.pop(iArg)

        Basemap.__init__(self, projection=projection, **kwargs)

        # create figure and set it as an attribute
        plt.close()
        self.fig = plt.figure(**figKwargs)

    def smooth(self, idata, mode, **kwargs):
        '''Smooth data for contour() and contourf()

        idata is smoothed by convolve, fourier_gaussian, spline or
        gaussian (default). If contour_mode is 'convolve' and weight is None,
        the weight matrix is created automatically.

        Parameters
        -----------
        idata : numpy 2D array
            Input data
        mode : string
            'convolve','fourier','spline' or 'gaussian'

        Returns
        ---------
        odata : numpy 2D array

        See also
        ----------
        http://docs.scipy.org/doc/scipy/reference/ndimage.html

        '''
        # modify default parameter
        self._set_defaults(kwargs)

        if mode == 'convolve':
            # if weight is None, create a weight matrix
            if self.convolve_weights is None:
                weights = np.ones((self.convolve_weightSize,
                                   self.convolve_weightSize))
                center = (self.convolve_weightSize - 1) / 2
                for i in range(- (center), center + 1, 1):
                    for j in range(- (center), center + 1, 1):
                        weights[i][j] /= pow(2.0, max(abs(i), abs(j)))
                self.convolve_weights = weights
            odata = ndimage.convolve(idata,
                                     weights=self.convolve_weights,
                                     mode=self.convolve_mode,
                                     cval=self.convolve_cval,
                                     origin=self.convolve_origin)
        elif mode == 'fourier':
            odata = ndimage.fourier_gaussian(idata,
                                             sigma=self.fourier_sigma,
                                             n=self.fourier_n,
                                             axis=self.fourier_axis)
        elif mode == 'spline':
            odata = ndimage.spline_filter1d(idata,
                                            order=self.spline_order,
                                            axis=self.spline_axis)
        else:
            if mode != 'gaussian':
                print 'apply Gaussian filter in image_process()'
            odata = ndimage.gaussian_filter(idata,
                                            sigma=self.gaussian_sigma,
                                            order=self.gaussian_order,
                                            mode=self.gaussian_mode,
                                            cval=self.gaussian_cval)
        return odata

    def _do_contour(self, bmfunc, data, v, smooth, mode, **kwargs):
        ''' Prepare data and make contour or contourf plots

        1. Smooth data
        1. Add colormap
        1. Append contour or contourf plot to self.mpl

        bmfunc : Basemap function
            Basemap.contour, Basemap.contourf
        data : numpy 2D array
            Input data
        v : list with values
            draw contour lines at the values specified in sequence v
        smooth : Boolean
            Apply smoothing?
        mode : string
            'gaussian', 'spline', 'fourier', 'convolve'
            mname of smoothing algorithm to apply

        '''
        self._create_xy_grids()

        # if cmap is given, set to self.cmap
        if 'cmap' in kwargs.keys():
            self.cmap = kwargs.pop('cmap')

        # smooth data
        if smooth:
            data = self.smooth(data, mode, **kwargs)

        # draw contour lines
        if v is None:
            self.mpl.append(bmfunc(self, self.x, self.y, data, **kwargs))
        else:
            self.mpl.append(bmfunc(self, self.x, self.y, data, v, **kwargs))

    def contour(self, data, v=None, smooth=False, mode='gaussian',
                label=True, **kwargs):
        '''Draw lined contour plots

        If smooth is True, data is smoothed. Then draw lined contour.

        Parameters
        ----------
        data : numpy 2D array
            Input data
        v : list with values
            draw contour lines at the values specified in sequence v
        smooth : Boolean
            Apply smoothing?
        mode : string
            'gaussian', 'spline', 'fourier', 'convolve'
            mname of smoothing algorithm to apply
        label : boolean
            Add lables?
        **kwargs:
            Optional parameters for Nansatmap.smooth()
            Optional parameters for pyplot.contour().
            Optional parameters for pyplot.clabel()

        Modifies
        ---------
        self.mpl : list
            append QuadContourSet instance
        '''

        self._do_contour(Basemap.contour, data, v, smooth, mode, **kwargs)

        # add lables to the contour lines
        if label:
            plt.clabel(self.mpl[-1], **kwargs)

    def contourf(self, data, v=None,
                 smooth=False, mode='gaussian', **kwargs):
        '''Draw filled contour plots

        If smooth is True, data is smoothed. Then draw filled contour.

        Parameters
        ----------
        data : numpy 2D array
            Input data
        v : list with values
            draw contour lines at the values specified in sequence v
        smooth : Boolean
            Apply smoothing?
        mode : string
            'gaussian', 'spline', 'fourier', 'convolve'
            mname of smoothing algorithm to apply
        **kwargs:
            cmap : colormap (e.g. cm.jet)
            Optional parameters for Nansatmap.smooth()
            Optional parameters for pyplot.contourf().

        Modifies
        ---------
        self.mpl : list
            append QuadContourSet instance

        '''
        self._do_contour(Basemap.contourf, data, v, smooth, mode, **kwargs)
        self.colorbar = len(self.mpl) - 1

    def imshow(self, data, low=0, high=255, **kwargs):
        ''' Make RGB plot over the map

        data : numpy array
            RGB or RGBA input data
        **kwargs:
            Parameters for Basemap.imshow

        Modifies
        ---------
        self.mpl : list
            append AxesImage object with imshow

        '''
        # Create X/Y axes
        self._create_xy_grids()

        # add random colormap
        if 'cmap' in kwargs and kwargs['cmap'] == 'random':
            values = np.unique(data[np.isfinite(data)])
            cmap, norm = self._create_random_colormap(values,
                                                      low=low, high=high)
            kwargs['cmap'] = cmap
            kwargs['norm'] = norm

        # Plot data using imshow
        self.mpl.append(Basemap.imshow(self, data,
                                       extent=[self.x.min(), self.x.max(),
                                               self.y.min(), self.y.max()],
                                       origin='upper', **kwargs))
        self.colorbar = len(self.mpl) - 1

    def pcolormesh(self, data, **kwargs):
        '''Make a pseudo-color plot over the map

        Parameters
        ----------
        data : numpy 2D array
            Input data
        **kwargs:
            Parameters for Basemap.pcolormesh (e.g. vmin, vmax)

        Modifies
        ---------
        self.mpl : list
            append matplotlib.collections.QuadMesh object

        '''
        # mask nan data
        data = np.ma.array(data, mask=np.isnan(data))
        # Plot a quadrilateral mesh.
        self._create_xy_grids()
        self.mpl.append(Basemap.pcolormesh(self, self.x, self.y, data,
                                           **kwargs))
        self.colorbar = len(self.mpl) - 1

    def quiver(self, dataX, dataY, step=None, quivectors=None, **kwargs):
        '''Draw quiver plots

        Parameters
        ----------
        dataX :  numpy array
            Input data with X-component
        dataY :  numpy array
            Input data with Y-component
        step : int or (int, int)
            Skip <step> pixels along both dimentions(alternative to quivectors)
        quivectors : int or (int,int)
            Number of vectors along both dimentions
        Parameters for Basemap.quiver()

        Modifies
        ---------
        self.mpl : list
            append matplotlib.quiver.Quiver instance

        '''
        # if Nan is included, apply mask
        dataX = np.ma.array(dataX, mask=np.isnan(dataX))
        dataY = np.ma.array(dataY, mask=np.isnan(dataY))

        # get subsetting parameters
        if type(step) is int:
            step0 = step1 = step
        elif type(step) in [list, tuple]:
            step0 = step[0]
            step1 = step[1]
        elif quivectors is not None:
            if type(quivectors) is int:
                quivectors0 = quivectors
                quivectors1 = quivectors
            if type(quivectors) in [list, tuple]:
                quivectors0 = quivectors[0]
                quivectors1 = quivectors[1]
            step0 = dataX.shape[0] / quivectors0
            step1 = dataX.shape[1] / quivectors1
        else:
            step0 = step1 = 5

        dataX2 = dataX[::step0, ::step1]
        dataY2 = dataY[::step0, ::step1]
        self._create_lonlat_grids()
        lon2 = self.lon[::step0, ::step1]
        lat2 = self.lat[::step0, ::step1]
        x2, y2 = self(lon2, lat2)

        qKwargs = {}
        for iKey in ['width', 'scale', 'units', 'angles', 'scale_units']:
            if iKey in kwargs.keys():
                qKwargs[iKey] = kwargs.pop(iKey)
        Q = Basemap.quiver(self, x2, y2, dataX2, dataY2, **qKwargs)

        qkargs = {}
        for iKey in ['X', 'Y', 'U', 'label']:
            if iKey in kwargs.keys():
                qkargs[iKey] = kwargs.pop(iKey)

        if all(iKey in qkargs.keys() for iKey in ('X', 'Y', 'U', 'label')):
            self.mpl.append(plt.quiverkey(Q, qkargs['X'], qkargs['Y'],
                                          qkargs['U'], qkargs['label'],
                                          **kwargs))
        else:
            self.mpl.append(Q)

    def add_colorbar(self, fontsize=6, **kwargs):
        '''Add color bar

        Parameters
        ----------
        fontsize : int
        Parameters for matplotlib.pyplot.colorbar

        Modifies
        ---------
        Adds colorbar to self.fig

        '''
        if kwargs is None:
            kwargs = {}
        if not ('orientation' in kwargs.keys()):
            kwargs['orientation'] = 'horizontal'
        if not ('pad' in kwargs.keys()):
            kwargs['pad'] = 0.01

        # add colorbar and set font size
        if self.colorbar is not None:
            origin = self.mpl[self.colorbar]

            # if colormap is ListedColormap
            # add integer ticks
            ticks = None
            listedColormap = False
            if (hasattr(origin, 'cmap') and
                (type(origin.cmap) == mpl.colors.ListedColormap) and
                hasattr(origin.norm, 'boundaries')) :
                ticks = (origin.norm.boundaries[:-1] +
                         np.diff(origin.norm.boundaries) / 2.)
                listedColormap = True

            cbar = self.fig.colorbar(origin, ticks=ticks, **kwargs)
            if listedColormap:
                labels = origin.norm.boundaries[:-1]
                if np.all(labels == np.floor(labels)):
                    labels = labels.astype('int32')
                cbar.ax.set_xticklabels(labels)
            imaxes = plt.gca()
            plt.axes(cbar.ax)
            plt.xticks(fontsize=fontsize)
            plt.axes(imaxes)

    def drawgrid(self, lat_num=5, lon_num=5,
                 lat_labels=[True, False, False, False],
                 lon_labels=[False, False, True, False],
                 **kwargs):
        '''Draw and label parallels (lat and lon lines) for values (in degrees)

        Parameters
        -----------
        fontsize : int
        lat_num : int
            Number of latitude lables
        lon_num :
            Number of longitude lables
        lat_labels : list of Bool
            Location of latitude labels
        lon_labels : list of Bool
            Location of longitude labels

        See also: Basemap.drawparallels(), Basemap.drawmeridians()

        '''
        self.drawparallels(np.arange(self.latMin, self.latMax,
                           (self.latMax - self.latMin) / lat_num),
                           labels=lat_labels, **kwargs)
        self.drawmeridians(np.arange(self.lonMin, self.lonMax,
                           (self.lonMax - self.lonMin) / lon_num),
                           labels=lon_labels, **kwargs)

    def draw_continents(self, **kwargs):
        ''' Draw continents

        Parameters
        ----------
        Parameters for basemap.fillcontinents

        '''

        if kwargs is None:
            kwargs = {}
        if not ('color' in kwargs.keys()):
            kwargs['color'] = '#999999'
        if not ('lake_color' in kwargs.keys()):
            kwargs['lake_color'] = '#99ffff'

        # draw continets
        self.fillcontinents(**kwargs)

    def save(self, fileName, landmask=True, dpi=75,
             pad_inches=0, bbox_inches='tight', **kwargs):
        '''Draw continents and save

        Parameters
        -----------
        fileName : string
            name of outputfile
        landmask : Boolean
            Draw landmask?
        Parameters for basemap.fillcontinents

        '''
        if landmask:
            self.draw_continents(**kwargs)

        # set default extension
        if not((fileName.split('.')[-1] in self.extensionList)):
            fileName = fileName + self.DEFAULT_EXTENSION
        self.fig.savefig(fileName, dpi=dpi,
                         pad_inches=pad_inches,
                         bbox_inches=bbox_inches)

    def _set_defaults(self, idict):
        '''Check input params and set defaut values

        Look throught default parameters (self.d) and given parameters (dict)
        and paste value from input if the key matches

        Parameters
        ----------
        idict : dictionary
            parameter names and values

        Modifies
        ---------
            default self attributes

        '''
        for key in idict:
            if hasattr(self, key):
                setattr(self, key, idict[key])

    def _create_lonlat_grids(self):
        '''Generate grids with lon/lat coordinates in each cell

        Modifies
        ---------
        self.lon : numpy array with lon coordinates
        self.lat : numpy array with lat coordinates
        '''
        if self.lon is None or self.lat is None:
            self.lon, self.lat = self.domain.get_geolocation_grids()

    def _create_xy_grids(self):
        '''Generate grids with x/y coordinates in each cell

        Modifies
        ---------
        self.x : numpy array with X coordinates
        self.y : numpy array with Y coordinates
        '''
        self._create_lonlat_grids()
        if self.x is None or self.y is None:
            self.x, self.y = self(self.lon, self.lat)

    def _create_random_colormap(self, values, low=0, high=255):
        ''' Generate colormap and colorbar with random discrete colors

        Parameters
        ----------
            values : list or 1D array
                values for which the random colors are to be generated
        Returns
        -------
            cmap : matplotlib.color.Colormap
            norm : matplotlib.color.BoundaryNorm
        '''
        # create first random color
        randomColors = [get_random_color(low=low, high=high)]
        # add more random colors
        for v in values[1:]:
            randomColors.append(get_random_color(randomColors[-1],
                                                 low=low, high=high))

        # create colormap and norm
        cmap = mpl.colors.ListedColormap(randomColors)
        bounds = sorted(list(values))
        bounds += [max(bounds) + 1]  # bounds should be longer than values by 1
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        return cmap, norm

    def add_zone_labels(self, zones, fontsize=5):
        ''' Finds best place of labels for a zone map, adds labels to the map

        Parameters
        ----------
            zones : numpy array with integer zones
                the same array as usied in Nansatmap.imshow
        '''
        zoneIndices = np.unique(zones[np.isfinite(zones)])
        for zi in zoneIndices:
            zrows, zcols = np.nonzero(zones == zi)
            zrc = np.median(zrows)
            zcc = np.median(zcols)
            lon, lat = self.domain.transform_points([zcc], [zrc], 0)
            x, y = self(lon[0], lat[0])
            plt.text(x, y, '%d' % zi, fontsize=fontsize)
