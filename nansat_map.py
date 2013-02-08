#-------------------------------------------------------------------------------
# Name:        nansatmap
# Purpose:
#
# Author:      asumak
#
# Created:     06.02.2013
# Copyright:   (c) asumak 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from nansat_tools import *
import time

class NansatMap():

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
        self.map : basemap object
        self.x, self.y :numpy arrays
            map projection coordinates
        self.fig : figure
            matplotlib.pyplot.figure
        self.colorbar : boolean
            if colorbar is True, it is possible to put colorbar.
            e.g. contour_plots(contour_style='fill'), put_color()
        self.mpl : matplotlib

        '''
        # get proj4
        spatialRef = osr.SpatialReference()
        projection = domain._get_projection(domain.vrt.dataset)
        spatialRef.ImportFromWkt(projection)
        proj4 = spatialRef.ExportToProj4()
        print proj4

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
        self.d['llcrnrx'] = None
        self.d['llcrnry'] = None
        self.d['urcrnrx'] = None
        self.d['urcrnry'] = None
        self.d['width'] = None
        self.d['height'] = None
        self.d['projection'] = projection
        self.d['resolution'] = 'c'
        self.d['area_thresh'] = None
        self.d['rsphere'] = 6370997.0
        self.d['lat_ts'] = None
        self.d['lat_0'] = None
        self.d['lat_1'] = None
        self.d['lat_2'] = None
        self.d['lon_0'] = None
        self.d['lon_1'] = None
        self.d['lon_2'] = None
        self.d['k_0'] = None
        self.d['no_rot'] = None
        self.d['suppress_ticks'] = True
        self.d['satellite_height'] = 35786000
        self.d['boundinglat'] = None
        self.d['fix_aspect'] = True
        self.d['anchor'] = 'C'
        self.d['celestial'] = False
        self.d['round'] = False
        self.d['ax'] = None
        self.d['continentColor'] = '#999999'
        self.d['lakeColor'] = '#99ffff'
        self.d['continentStyle'] = None
        # figure
        self.d['fignum'] = 1
        self.d['figsize'] = None
        self.d['dpi'] = None
        self.d['facecolor'] = None
        self.d['edgecolor'] = None
        self.d['linewidth'] = 0.0
        self.d['frameon'] = True
        self.d['subplotpars'] = None
        # pcolormesh
        self.d['color_data'] = None
        self.d['cmap'] = plt.cm.jet
        self.d['norm'] = None
        self.d['vmin'] = None
        self.d['vmax'] = None
        self.d['shading'] = 'flat'
        self.d['edgecolors'] = None
        self.d['alpha'] = None
        # contour_plots
        self.d['contour_data'] = None
        self.d['contour_style'] = 'line'
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
        self.d['quiver_data'] = [None, None]
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
        print "getting geolocation grids..."
        start_time = time.time()
        self.lon, self.lat = domain.get_geolocation_grids()
        print "Completed get_geolocation_grids", time.time() - start_time

        # set llcrnrlat and urcrnrlat and
        # if required, modify them from -90. to 90.
        self.d['llcrnrlat'] = max(self.lat.min(), -90.)
        self.d['urcrnrlat'] = min(self.lat.max(), 90.)
        # set llcrnrlon and urcrnrlon from self.lon
        self.d['llcrnrlon'] = self.lon.min()
        self.d['urcrnrlon'] = self.lon.max()

        # modify the default values using input values
        self._set_defaults(kwargs)

        print "creating basemap..."
        start_time = time.time()
        # create basemap object
        self.map = Basemap(llcrnrlon=self.d['llcrnrlon'],
                           llcrnrlat=self.d['llcrnrlat'],
                           urcrnrlon=self.d['urcrnrlon'],
                           urcrnrlat=self.d['urcrnrlat'],
                           llcrnrx=self.d['llcrnrx'],
                           llcrnry=self.d['llcrnry'],
                           urcrnrx=self.d['urcrnrx'],
                           urcrnry=self.d['urcrnry'],
                           width=self.d['width'],
                           height=self.d['height'],
                           projection=self.d['projection'],
                           resolution=self.d['resolution'],
                           area_thresh=self.d['area_thresh'],
                           rsphere=self.d['rsphere'],
                           lat_ts=self.d['lat_ts'],
                           lat_0=self.d['lat_0'],
                           lat_1=self.d['lat_1'],
                           lat_2=self.d['lat_2'],
                           lon_0=self.d['lon_0'],
                           lon_1=self.d['lon_1'],
                           lon_2=self.d['lon_2'],
                           k_0=self.d['k_0'],
                           no_rot=self.d['no_rot'],
                           suppress_ticks=self.d['suppress_ticks'],
                           satellite_height=self.d['satellite_height'],
                           boundinglat=self.d['boundinglat'],
                           fix_aspect=self.d['fix_aspect'],
                           anchor=self.d['anchor'],
                           celestial=self.d['celestial'],
                           round=self.d['round'],
                           ax=self.d['ax'])
        print "Completed basemap", time.time() - start_time
        # convert to map projection coords and set them to x and y
        self.x, self.y = self.map(self.lon, self.lat)

        # create figure
        plt.close()
        print "creating figure..."
        start_time = time.time()
        self.fig = plt.figure( num=self.d['fignum'],
                               figsize=self.d['figsize'],
                               dpi=self.d['dpi'],
                               facecolor=self.d['facecolor'],
                               edgecolor=self.d['edgecolor'],
                               linewidth=self.d['linewidth'],
                               frameon=self.d['frameon'],
                               subplotpars=self.d['subplotpars'])
        print "Completed figure", time.time() - start_time
        # set attributes
        self.colorbar = False
        self.mpl = None

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
        self.mpl : matplotlib

        see also:
        ---------
        http://docs.scipy.org/doc/scipy/reference/ndimage.html
        http://matplotlib.org/basemap/api/basemap_api.html

        '''
        # modify default parameter
        self._set_defaults(kwargs)
        # if direct argument 'data' is None, set from kwargs.
        if data is None:
            data =self.d.clear['contour_data']

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
                            weights[i][j] /= math.pow(2.0, max(abs(i),abs(j)))
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
            self.mpl = self.map.contourf(self.x, self.y, data,
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
                                   hatches=self.d['contour_hatches'])

            self.colorbar = True
        else:
            # draw contour lines
            self.mpl = self.map.contour(self.x, self.y, data,
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
                                   linestyles=self.d['contour_linestyles'])

            # add values for the contour lines
            if self.d['contour_label']:
                plt.clabel(self.mpl,
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
        self.mpl : matplotlib

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
        x2, y2 = self.map(lon2, lat2)
        self.mpl = self.map.quiver(x2, y2, dataX2, dataY2)

    def put_color(self, data=None, **kwargs):
        '''Make a pseudo-color plot over the map

        Parameters
        ----------
        data : numpy 2D array
        Any of NansatMap.__init__() parameters

        Modifies
        ---------
        self.mpl : matplotlib

        '''
        # modify default parameter
        self._set_defaults(kwargs)
        # if direct argument 'data' is None, set from kwargs.
        if data is None:
            data = self.d['color_data']
        # Plot a quadrilateral mesh.
        self.mpl = self.map.pcolormesh(self.x, self.y, data,
                            cmap=self.d['cmap'], norm=self.d['norm'],
                            vmin=self.d['vmin'], vmax=self.d['vmax'],
                            shading=self.d['shading'],
                            edgecolors=self.d['edgecolors'],
                            alpha=self.d['alpha'])
        self.colorbar = True

    def draw_continent(self, **kwargs):
        '''Draw continent

        Parameters
        ----------
        Any of NansatMap.__init__() parameters
        self.d['continentStyle'] takes 'mask', 'bluemarble' or 'fill'.
        If it is not either None, 'mask' or 'bluemarble', then 'fill'.

        Modifies
        ---------
        self.map : basemap

        NB: Unfortunately, the 'fill' method does not always do the right thing.
            Matplotlib always tries to fill the inside of a polygon.
            Under certain situations, what is the inside of a coastline
            polygon can be ambiguous, and the outside may be filled
            instead of the inside. In these situations, the recommended
            workaround is to use the 'mask'.
            But Landmask and Bluemarble image cannot be overlaid on top of
            other images due to limitations in matplotlib image handling.

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        # draw continent
        if self.d['continentStyle']=='mask':
            self.map.drawlsmask(self.d['continentColor'], self.d['lakeColor'],
                                lakes=True)
        elif self.d['continentStyle']=='bluemarble':
            self.map.bluemarble()
        else:
            self.map.fillcontinents(color=self.d['continentColor'],
                                    lake_color=self.d['lakeColor'])

    def draw_geoCoordinates(self, **kwargs):
        '''Add longitude and latitude

        Parameters
        ----------
        Any of NansatMap.__init__() parameters

        Modifies
        ---------
        self.map : basemap

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        # draw lat and lon
        self.map.drawparallels(np.arange(self.lat.min(), self.lat.max(),
                               (self.lat.max()-self.lat.min())/self.d['lat_num']),
                               labels=self.d['lat_labels'],
                               fontsize=self.d['lat_fontsize'])
        self.map.drawmeridians(np.arange(self.lon.min(), self.lon.max(),
                               (self.lon.max()-self.lon.min())/self.d['lon_num']),
                               labels=self.d['lon_labels'],
                               fontsize=self.d['lon_fontsize'])

    def add_legend(self, **kwargs):
        '''Add color bar and title

        Parameters
        ----------
        Any of NansatMap.__init__() parameters

        '''
        # modify default parameters
        self._set_defaults(kwargs)
        # add colorbar and reduce font size
        if self.colorbar and self.d['colorbar']:
            cbar = self.fig.colorbar(self.mpl,
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

        # if d['contour_data'] is not None, draw contour plots.
        if self.d['contour_data'] is not None:
            self.contour_plots()

        # if d['quiver_data'] is not None, draw quiver plots.
        if self.d['quiver_data'] is not None:
            self.quiver_plots()

        # if d['color_data'] is not None, put colors on the map
        if self.d['color_data'] is not None:
            self.put_color()

        # if d['continentStyle'] is 'mask', 'bluemarble' or 'fill',
        # then draw continets
        if self.d['continentStyle'] is not None:
            self.draw_continent()

        # if d['geocoordinates'] is True, then draw geoCoordinates
        if self.d['geocoordinates']:
            self.draw_geoCoordinates()

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
        #plt.show()

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





