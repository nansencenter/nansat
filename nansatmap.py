# Name:    nansat_map.py
# Purpose: Container of NansatMap class
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

from nansat_tools import *

class Nansatmap(Basemap):
    '''Perform opeartions with graphical files: create,
    add legend and geolocation_grids, save.

    NansatMap instance is created in the Nansat.write_map method.
    The methods below are applied consequently in order to get projection,
    generate a basemap from array(s), add legend and geolocation grids,
    save to a file.

    '''
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
        self.lon, self.lat : numpy arrays
            lat and lon of the domain in degrees
        self.x, self.y :numpy arrays
            map projection coordinates
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
        # get proj4
        spatialRef = osr.SpatialReference()
        projection = domain._get_projection(domain.vrt.dataset)
        spatialRef.ImportFromWkt(projection)
        proj4 = spatialRef.ExportToProj4()

        # convert proj4 to basemap projection
        projStr = proj4.split(' ')[0][6:]
        projection = {  'aea':'aea', 'ocea':'aea',
                        'aeqd':'aeqd', 'xxx1':'spaeqd', 'xxx2':'npaeqd',
                        'cass':'cass',
                        'cea':'cea',
                        'eqc':'cyl', 'longlat':'cyl',
                        'eck4':'eck4',
                        'eqdc':'eqdc',
                        'gall':'gall',
                        'geos':'geos',
                        'gnom':'gnom',
                        'hammer':'hammer', 'nell_h':'hammer',
                        'kav7':'kav7',
                        'laea':'laea', 'xxx3':'splaea', 'xxx4':'nplaea',
                        'lcc':'lcc', 'lcca':'lcc',
                        'mbtfpq':'mbtfpq',
                        'somerc':'merc', 'merc':'merc', 'omerc':'merc',
                        'mill':'mill',
                        'moll':'moll',
                        'nsper':'nsper',
                        'omerc':'omerc',
                        'ortho':'ortho',
                        'poly':'poly', 'rpoly':'poly', 'imw_p':'poly',
                        'robin':'robin',
                        'sinu':'sinu', 'fouc_s':'sinu', 'gn_sinu':'sinu',
                        'mbtfps':'sinu','urmfps':'sinu',
                        'stere':'stere', 'sterea':'stere', 'lee_os':'stere',
                        'mil_os':'stere', 'rouss':'stere',
                        'ups':'npstere', 'ups':'spstere', # CHECK!! #
                        'tmerc':'tmerc', 'gstmerc':'tmerc', 'utm':'tmerc',
                        'vandg':'vandg', 'vandg2':'vandg',
                        'vandg3':'vandg', 'vandg4':'vandg',
                     }.get(projStr, 'cyl')

        # set default values of ALL params of NansatMap
        self.d = {}
        self.d['llcrnrlon'] = None
        self.d['llcrnrlat'] = None
        self.d['urcrnrlon'] = None
        self.d['urcrnrlat'] = None
        self.d['projection'] = projection
        # fillcontinents
        self.d['continentColor'] = '#999999'
        self.d['lakeColor'] = '#99ffff'
        self.d['continent'] = True
        # figure
        self.d['fignum'] = 1
        # pcolormesh
        self.d['color_data'] = None
        # contour_plots
        self.d['contour_data'] = None
        self.d['contour_style'] = 'line' # 'fill' or 'line'
        self.d['contour_smoothing'] = False
        self.d['contour_mode'] = 'gaussian'
        self.d['contour_alpha'] = None
        self.d['contour_label'] = False
         # contourf & contour
        self.d['contour_colors'] = None
        self.d['contour_alpha'] = None
        self.d['contour_cmap'] = None
        self.d['contour_norm'] = None
        self.d['contour_vmin'] = None
        self.d['contour_vmax'] = None
        self.d['contour_levels'] = None
        self.d['contour_origin'] = None
        self.d['contour_extent'] = None
        self.d['contour_locator'] = None
        self.d['contour_extend'] = 'neither'
        self.d['contour_xunits'] = None
        self.d['contour_yunits'] = None
        self.d['contour_antialiased'] = True
        # contour
        self.d['contour_linewidths'] = None
        self.d['contour_linestyles'] = None
        # contour line label
        self.d['contour_linesfontsize'] = 3
        self.d['contour_inline'] = 1
        # contourf
        self.d['contour_nchunk'] = 0
        self.d['contour_hatches'] = None
        # convolve
        self.d['convolve_weightSize'] = 7
        self.d['convolve_weights'] = None
        self.d['convolve_mode'] = 'reflect'
        self.d['convolve_cval'] = 0.0
        self.d['convolve_origin'] = 0
        # fourier_gaussian
        self.d['fourier_sigma'] = 1.0
        self.d['fourier_n'] = -1
        self.d['fourier_axis'] = -1
        # spline
        self.d['spline_order'] = 3
        self.d['spline_axis'] = -1
        # gaussian filter
        self.d['gaussian_sigma'] = 2.5
        self.d['gaussian_order'] = 0
        self.d['gaussian_mode'] = 'reflect'
        self.d['gaussian_cval'] = 0.0
        # quiver_plots
        self.d['quiver_data'] = None
        self.d['quivectors'] = 30
        # legend (title)
        self.d['title'] = ''
        self.d['title_fontsize'] = 7
        # legend (color bar)
        self.d['colorbar'] = True
        self.d['colorbar_orientation'] = 'horisontal'
        self.d['colorbar_pad'] = 0.01
        self.d['colorbar_fontsize'] = 6
        # draw_geoCoordinates
        self.d['geocoordinates'] = False
        self.d['lat_num'] = 5
        self.d['lat_fontsize'] = 4
        self.d['lat_labels'] = [True, False, False, False]
        self.d['lon_num'] = 5
        self.d['lon_fontsize'] = 4
        self.d['lon_labels'] = [False, False, True, False]

        # set lon and lat attributes from nansat
        self.lon, self.lat = domain.get_geolocation_grids()

        # set llcrnrlat and urcrnrlat and
        # if required, modify them from -90. to 90.
        self.d['llcrnrlat'] = max(self.lat.min(), -90.)
        self.d['urcrnrlat'] = min(self.lat.max(), 90.)
        # set llcrnrlon and urcrnrlon from self.lon
        self.d['llcrnrlon'] = self.lon.min()
        self.d['urcrnrlon'] = self.lon.max()

        # modify the default values using input values
        self._set_defaults(kwargs)

        Basemap.__init__(self,
                         llcrnrlon=self.d['llcrnrlon'],
                         llcrnrlat=self.d['llcrnrlat'],
                         urcrnrlon=self.d['urcrnrlon'],
                         urcrnrlat=self.d['urcrnrlat'],
                         projection=self.d['projection'],
                         **kwargs)

        # convert to map projection coords and set them to x and y
        self.x, self.y = self(self.lon, self.lat)

        # create figure and set it as an attribute
        plt.close()
        self.fig = plt.figure(num=self.d['fignum'], **kwargs)

        # set attributes
        self.colorbar = None
        self.mpl = []

    def contour_plots(self, data=None, **kwargs):
        '''Draw contour plots
        If contour_smoothing is True, data is smoothed by convolve,
        fourier_gaussian, spline or gaussian (default).
        If contour_mode is 'convolve' and weight is None,
        the weight matrix is created automatically.
        If contour_style is 'fill', then filled contour.
        Otherwise, lined contour.

        Parameters
        ----------
        data : numpy 2D array
        Any of NansatMap.__init__() parameters

        Modifies
        ---------
        self.mpl : list
            append QuadContourSet instance

        see also
        ---------
        http://docs.scipy.org/doc/scipy/reference/ndimage.html

        '''
        # modify default parameter
        self._set_defaults(kwargs)
        # if direct argument 'data' is None, set from kwargs.
        if data is None:
            data =self.d['contour_data']

        # smooth data
        if self.d['contour_smoothing']:
            if self.d['contour_mode']=='convolve':
                # if weight is None, create a weight matrix
                if self.d['convolve_weights'] is None:
                    weights = np.ones((self.d['convolve_weightSize'],
                                       self.d['convolve_weightSize']))
                    center = (self.d['convolve_weightSize'] - 1) / 2
                    for i in range(-(center), center+1, 1):
                        for j in range(-(center), center+1, 1):
                            weights[i][j] /= pow(2.0, max(abs(i),abs(j)))
                    self.d['convolve_weights'] = weights
                data = ndimage.convolve(data,
                                        weights=self.d['convolve_weights'],
                                        mode=self.d['convolve_mode'],
                                        cval=self.d['convolve_cval'],
                                        origin=self.d['convolve_origin'])
            elif self.d['contour_mode']=='fourier_gaussian':
                data = ndimage.fourier_gaussian(data,
                                                sigma=self.d['fourier_sigma'],
                                                n=self.d['fourier_n'],
                                                axis=self.d['fourier_axis'])
            elif self.d['contour_mode']=='spline':
                data = ndimage.spline_filter1d(data,
                                               order=self.d['spline_order'],
                                               axis=self.d['spline_axis'])
            else:
                if self.d['contour_mode']!='gaussian':
                    print 'apply Gaussian filter in image_process()'
                data = ndimage.gaussian_filter(data,
                                               sigma=self.d['gaussian_sigma'],
                                               order=self.d['gaussian_order'],
                                               mode=self.d['gaussian_mode'],
                                               cval=self.d['gaussian_cval'])
        if self.d['contour_style']=='fill':
            # draw filled contour
            self.mpl.append(self.contourf(self.x, self.y, data,
                                   colors=self.d['contour_colors'],
                                   alpha=self.d['contour_alpha'],
                                   cmap=self.d['contour_cmap'],
                                   norm=self.d['contour_norm'],
                                   vmin=self.d['contour_vmin'],
                                   vmax=self.d['contour_vmax'],
                                   levels=self.d['contour_levels'],
                                   origin=self.d['contour_origin'],
                                   extent=self.d['contour_extent'],
                                   locator=self.d['contour_locator'],
                                   extend=self.d['contour_extend'],
                                   xunits=self.d['contour_xunits'],
                                   yunits=self.d['contour_yunits'],
                                   antialiased=self.d['contour_antialiased'],
                                   nchunk=self.d['contour_nchunk'],
                                   hatches=self.d['contour_hatches']))
            self.colorbar = len(self.mpl)-1
        else:
            # draw contour lines
            self.mpl.append(self.contour(self.x, self.y, data,
                                   colors=self.d['contour_colors'],
                                   alpha=self.d['contour_alpha'],
                                   cmap=self.d['contour_cmap'],
                                   norm=self.d['contour_norm'],
                                   vmin=self.d['contour_vmin'],
                                   vmax=self.d['contour_vmax'],
                                   levels=self.d['contour_levels'],
                                   origin=self.d['contour_origin'],
                                   extent=self.d['contour_extent'],
                                   locator=self.d['contour_locator'],
                                   extend=self.d['contour_extend'],
                                   xunits=self.d['contour_xunits'],
                                   yunits=self.d['contour_yunits'],
                                   antialiased=self.d['contour_antialiased'],
                                   linewidths=self.d['contour_linewidths'],
                                   linestyles=self.d['contour_linestyles']))
            # add values for the contour lines
            if self.d['contour_label']:
                plt.clabel(self.mpl[-1],
                           inline=self.d['contour_inline'],
                           fontsize=self.d['contour_linesfontsize'])

    def quiver_plots(self, data=None, **kwargs):
        '''Draw quiver plots

        Parameters
        ----------
        data : list
            elements are numpy array
        Any of NansatMap.__init__() parameters

        Modifies
        ---------
        self.mpl : list
            append matplotlib.quiver.Quiver instance

        Raises
        --------
        TypeError : occurs if data is not a list or len(data) is not two

        '''
        # modify default parameter
        self._set_defaults(kwargs)
        # if direct argument 'data' is None, set from kwargs.
        if data is None:
            data = self.d['quiver_data']
        if type(data)!=list or len(data)!=2:
            raise TypeError('data for quiver_plots() must be '
                            'a list as: [np.array, np.array]')

        # subsample for quiver plot
        step0 = data[0].shape[0]/self.d['quivectors']
        step1 = data[0].shape[1]/self.d['quivectors']
        dataX2 = data[0][::step0, ::step1]
        dataY2 = data[1][::step0, ::step1]
        lon2 = self.lon[::step0, ::step1]
        lat2 = self.lat[::step0, ::step1]
        x2, y2 = self(lon2, lat2)
        self.mpl.append(self.quiver(x2, y2, dataX2, dataY2))

    def put_color(self, data=None, **kwargs):
        '''Make a pseudo-color plot over the map

        Parameters
        ----------
        data : numpy 2D array
        Any of NansatMap.__init__() parameters

        Modifies
        ---------
        self.mpl : list
            append matplotlib.collections.QuadMesh object

        '''
        # modify default parameter
        self._set_defaults(kwargs)
        # if direct argument 'data' is None, set from kwargs.
        if data is None:
            data = self.d['color_data']
        # Plot a quadrilateral mesh.
        self.mpl.append(self.pcolormesh(self.x, self.y, data, **kwargs))
        self.colorbar = len(self.mpl)-1

    def add_legend(self, **kwargs):
        '''Add color bar and title

        Parameters
        ----------
        Any of NansatMap.__init__() parameters

        '''
        # modify default parameters
        self._set_defaults(kwargs)
        # add colorbar and reduce font size
        if self.colorbar is not None and self.d['colorbar']:
            cbar = self.fig.colorbar(self.mpl[self.colorbar],
                                     orientation=self.d['colorbar_orientation'],
                                     pad=self.d['colorbar_pad'])
            imaxes = plt.gca()
            plt.axes(cbar.ax)
            plt.xticks(fontsize=self.d['colorbar_fontsize'])
            plt.axes(imaxes)

        # add title
        if self.d['title'] != '':
            plt.title(self.d['title'], fontsize=self.d['title_fontsize'])

    def process(self, **kwargs):
        '''Do all common operations for preparation of a map for saving

        #. Draw contour plots (line/fill) if requied
        #. Darw quiver plots if required
        #. Put colors on the map if required
        #. Draw continents if required
        #. Add geo-coordinates if required
        #. Add colorBar / title if required

        Parameters
        -----------
        Any of NansatMap.__init__() parameters

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        # if d['color_data'] is not None, put colors on the map
        if self.d['color_data'] is not None:
            self.put_color()

        # if d['contour_data'] is not None, draw contour plots.
        if self.d['contour_data'] is not None:
            self.contour_plots()

        #if d['quiver_data'] is not None, draw quiver plots.
        if self.d['quiver_data'] is not None:
            self.quiver_plots()

        # if d['continent'] is True, then draw continets
        if self.d['continent']:
            self.fillcontinents(color=self.d['continentColor'],
                                lake_color=self.d['lakeColor'])

        # if d['geocoordinates'] is True, then draw geoCoordinates
        if self.d['geocoordinates']:
            self.drawparallels(np.arange(self.lat.min(), self.lat.max(),
                               (self.lat.max()-self.lat.min())/self.d['lat_num']),
                               labels=self.d['lat_labels'],
                               fontsize=self.d['lat_fontsize'])
            self.drawmeridians(np.arange(self.lon.min(), self.lon.max(),
                               (self.lon.max()-self.lon.min())/self.d['lon_num']),
                               labels=self.d['lon_labels'],
                               fontsize=self.d['lon_fontsize'])

        # if d['colorbar'] is True or d['title'] is not empty, then add legend
        if self.d['colorbar'] or self.d['title'] != '':
            self.add_legend()

    def save(self, fileName):
        ''' Save self.fig to a physical file

        Parameters
        ----------
        fileName : string
            name of outputfile

        '''
        self.fig.savefig(fileName)

    def _set_defaults(self, kwargs):
        '''Check input params and set defaut values

        Look throught default parameters (self.d) and given parameters (kwargs)
        and paste value from input if the key matches

        Parameters
        ----------
        kwargs : dictionary
            parameter names and values

        Modifies
        ---------
        self.d

        '''
        for key in kwargs:
            if key in self.d:
                self.d[key] = kwargs[key]





