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
        # figure
        self.d['fignum'] = 1
        self.d['figsize'] = None
        self.d['dpi'] = None
        self.d['facecolor'] = None
        self.d['edgecolor'] = None
        self.d['frameon'] = True
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
        # save
        self.d['DEFAULT_EXTENSION'] = '.png'

        # set lon and lat attributes from nansat
        self.lon, self.lat = domain.get_geolocation_grids()

        # set llcrnrlat and urcrnrlat and
        # if required, modify them from -90. to 90.
        self.d['llcrnrlat'] = max(self.lat.min(), -90.)
        self.d['urcrnrlat'] = min(self.lat.max(), 90.)
        # set llcrnrlon and urcrnrlon from self.lon
        self.d['llcrnrlon'] = self.lon.min()
        self.d['urcrnrlon'] = self.lon.max()

        self.extensionList = ['png', 'emf', 'eps', 'pdf', 'rgba',
                              'ps', 'raw', 'svg', 'svgz']

        # modify the default values using input values
        kwargs = self._set_defaults(kwargs)
        if kwargs is None:
            kwargs = {}

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
        self.fig = plt.figure(num=self.d['fignum'],
                              figsize=self.d['figsize'],
                              dpi=self.d['dpi'],
                              facecolor=self.d['facecolor'],
                              edgecolor=self.d['edgecolor'],
                              frameon=self.d['frameon'])

        # set attributes
        self.colorbar = None
        self.mpl = []

    def smooth(self, idata, mode, **kwargs):
        '''Smooth data for contour plots'''
        
        # modify default parameter
        self._set_defaults(kwargs)
        
        if mode == 'convolve':
            # if weight is None, create a weight matrix
            if self.d['convolve_weights'] is None:
                weights = np.ones((self.d['convolve_weightSize'],
                                   self.d['convolve_weightSize']))
                center = (self.d['convolve_weightSize'] - 1) / 2
                for i in range(-(center), center+1, 1):
                    for j in range(-(center), center+1, 1):
                        weights[i][j] /= pow(2.0, max(abs(i),abs(j)))
                self.d['convolve_weights'] = weights
            odata = ndimage.convolve(idata,
                                    weights=self.d['convolve_weights'],
                                    mode=self.d['convolve_mode'],
                                    cval=self.d['convolve_cval'],
                                    origin=self.d['convolve_origin'])
        elif mode == 'fourier':
            odata = ndimage.fourier_gaussian(idata,
                                            sigma=self.d['fourier_sigma'],
                                            n=self.d['fourier_n'],
                                            axis=self.d['fourier_axis'])
        elif mode == 'spline':
            odata = ndimage.spline_filter1d(idata,
                                           order=self.d['spline_order'],
                                           axis=self.d['spline_axis'])
        else:
            if mode != 'gaussian':
                print 'apply Gaussian filter in image_process()'
            odata = ndimage.gaussian_filter(idata,
                                           sigma=self.d['gaussian_sigma'],
                                           order=self.d['gaussian_order'],
                                           mode=self.d['gaussian_mode'],
                                           cval=self.d['gaussian_cval'])
        return odata

    def contour(self, data, smooth=False, mode='gaussian', label=True, inline=True, fontsize=3, **kwargs):
        '''Draw lined contour plots
        If contour_smoothing is True, data is smoothed by convolve,
        fourier_gaussian, spline or gaussian (default).
        If contour_mode is 'convolve' and weight is None,
        the weight matrix is created automatically.

        Parameters
        ----------
        data : numpy 2D array
            Input data
        smooth : Boolean
            Apply smoothing?
        mode : string
            'gaussian', 'spline', 'fourier', 'convolve'
        label : Boolean
            Add labels?
        inline : Boolean
            Lables should be inline?
        fontsize : int
            Size of label font
        Parameters for Nansatmap.smooth()

        Modifies
        ---------
        self.mpl : list
            append QuadContourSet instance

        see also
        ---------
        http://docs.scipy.org/doc/scipy/reference/ndimage.html

        '''
        # smooth data
        if smooth:
            data = self.smooth(data, mode, **kwargs)

        # draw contour lines
        self.mpl.append(Basemap.contour(self, self.x, self.y, data, **kwargs))
        
        # add lables to the contour lines
        if label:
            plt.clabel(self.mpl[-1],
                       inline=inline,
                       fontsize=fontsize)

    def contourf(self, data, smooth=False, mode='gaussian', **kwargs):
        '''Draw filled contour plots

        If contour_smoothing is True, data is smoothed by convolve,
        fourier_gaussian, spline or gaussian (default).
        If contour_mode is 'convolve' and weight is None,
        the weight matrix is created automatically.

        Parameters
        ----------
        data : numpy 2D array
        data : numpy 2D array
            Input data
        smooth : Boolean
            Apply smoothing?
        mode : string
            'gaussian', 'spline', 'fourier', 'convolve'
        Parameters for Nansatmap.smooth()

        Modifies
        ---------
        self.mpl : list
            append QuadContourSet instance

        see also
        ---------
        http://docs.scipy.org/doc/scipy/reference/ndimage.html

        '''
        # smooth data
        if smooth:
            data = self.smooth(data, mode, **kwargs)
            
        # draw filled contour
        self.mpl.append(Basemap.contourf(self, self.x, self.y, data, **kwargs))
        self.colorbar = len(self.mpl)-1
                           
    def quiver(self, dataX, dataY, quivectors=30, **kwargs):
        '''Draw quiver plots

        Parameters
        ----------
        dataX :  numpy array
            Input data with X-component
        dataY :  numpy array
            Input data with Y-component
        quivectors : int
            Number of vectors along both dimentions
        Parameters for Basemap.quiver()

        Modifies
        ---------
        self.mpl : list
            append matplotlib.quiver.Quiver instance

        Raises
        --------
        TypeError : occurs if data is not a list or len(data) is not two

        '''
        # subsample for quiver plot
        step0 = dataX.shape[0]/quivectors
        step1 = dataX.shape[1]/quivectors
        dataX2 = dataX[::step0, ::step1]
        dataY2 = dataY[::step0, ::step1]
        lon2 = self.lon[::step0, ::step1]
        lat2 = self.lat[::step0, ::step1]
        x2, y2 = self(lon2, lat2)
        self.mpl.append(Basemap.quiver(self, x2, y2, dataX2, dataY2, **kwargs))

    def pcolormesh(self, data, **kwargs):
        '''Make a pseudo-color plot over the map

        Parameters
        ----------
        data : numpy 2D array
            Input data
        Parameters for Basemap.pcolormesh

        Modifies
        ---------
        self.mpl : list
            append matplotlib.collections.QuadMesh object

        '''
        # Plot a quadrilateral mesh.
        self.mpl.append(Basemap.pcolormesh(self, self.x, self.y, data, **kwargs))
        self.colorbar = len(self.mpl)-1

    def add_colorbar(self, orientation='horisontal', pad=0.01, fontsize=6):
        
        '''Add color bar and title

        Parameters
        ----------
        orientation : str
            'horisontal', 'vertical'
        pad : float
        fontsize : int

        Modifies
        ---------
        Adds colorbar to self.fig
        
        '''
        # add colorbar and set font size
        if self.colorbar is not None:
            cbar = self.fig.colorbar(self.mpl[self.colorbar],
                                     orientation=orientation,
                                     pad=pad)
            imaxes = plt.gca()
            plt.axes(cbar.ax)
            plt.xticks(fontsize=fontsize)
            plt.axes(imaxes)

    def drawgrid(self, fontsize=10, lat_num=5, lon_num=5, lat_labels=[True,False,False,False], lon_labels=[False,False,True,False]):
        '''Prepare map for saving and save

        #. Draw continents if required
        #. Add geo-coordinates if required
        #. Add colorBar / title if required
        #. Save self.fig to a physical file

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
        self.drawparallels(np.arange(self.lat.min(), self.lat.max(),
                           (self.lat.max() - self.lat.min()) / lat_num),
                           labels=lat_labels,
                           fontsize=fontsize)
        self.drawmeridians(np.arange(self.lon.min(), self.lon.max(),
                           (self.lon.max() - self.lon.min()) / lon_num),
                           labels=lon_labels,
                           fontsize=fontsize)
    
    def save(self, fileName, landmask=True, **kwargs):
        '''Prepare map for saving and save

        #. Draw continents if required
        #. Add geo-coordinates if required
        #. Add colorBar / title if required
        #. Save self.fig to a physical file

        Parameters
        -----------
        fileName : string
            name of outputfile
        landmask : Boolean
            Draw landmask?
        Any of NansatMap.__init__() parameters

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        # draw continets
        if landmask:
            self.fillcontinents(color=self.d['continentColor'],
                                lake_color=self.d['lakeColor'])

        # set default extension
        if not((fileName.split('.')[-1] in self.extensionList)):
            fileName = fileName + self.d['DEFAULT_EXTENSION']

        self.fig.savefig(fileName)

    def _set_defaults(self, kwargs):
        '''Check input params and set defaut values

        Look throught default parameters (self.d) and given parameters (kwargs)
        and paste value from input if the key matches

        Parameters
        ----------
        kwargs : dictionary
            parameter names and values

        Returns
        --------
        kwargs : dictionary

        Modifies
        ---------
        self.d

        '''
        keys = kwargs.keys()
        for iKey in keys:
            if iKey in self.d:
                self.d[iKey] = kwargs.pop(iKey)
        return kwargs
