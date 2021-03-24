# Name:    domain.py
# Purpose: Container of Domain class
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov, Aleksander Vines
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
from __future__ import division, absolute_import

import re
import warnings
from xml.etree.ElementTree import ElementTree

import numpy as np

from nansat.utils import add_logger, initial_bearing, haversine, gdal, osr, ogr
from nansat.nsr import NSR
from nansat.vrt import VRT
from nansat.exceptions import NansatProjectionError
from nansat.warnings import NansatFutureWarning

class Domain(object):
    """Container for geographical reference of a raster

    A Domain object describes all attributes of geographical
    reference of a raster:

      * width and height (number of pixels)
      * pixel size (e.g. in decimal degrees or in meters)
      * relation between pixel/line coordinates and geographical
        coordinates (e.g. a linear relation)
      * type of data projection (e.g. geographical or stereographic)

    Parameters
    ----------
    srs : PROJ4 or EPSG or WKT or NSR or osr.SpatialReference()
        Input parameter for nansat.NSR()
    ext : string
        some gdalwarp options + additional options
        [http://www.gdal.org/gdalwarp.html]
        Specifies extent, resolution / size
        Available options: (('-te' or '-lle') and ('-tr' or '-ts'))
        (e.g. '-lle -10 30 55 60 -ts 1000 1000' or
        '-te 100 2000 300 10000 -tr 300 200')
        -tr resolutionx resolutiony
        -ts sizex sizey
        -te xmin ymin xmax ymax
        -lle lonmin latmin lonmax latmax
    ds : GDAL dataset

    Examples
    --------
        >>> d = Domain(srs, ext) #size, extent and spatial reference is given by strings
        >>> d = Domain(ds=GDALDataset) #size, extent copied from input GDAL dataset
        >>> d = Domain(srs, ds=GDALDataset) # spatial reference is given by srs,
            but size and extent is determined from input GDAL dataset

    Notes
    -----
    The core of Domain is a GDAL Dataset. It has no bands, but only
    georeference information: rasterXsize, rasterYsize, GeoTransform and
    Projection or GCPs, etc. which fully describe dimentions and spatial
    reference of the grid.

    There are three ways to store geo-reference in a GDAL dataset:

      * Using GeoTransfrom to define linear relationship between raster
        pixel/line and geographical X/Y coordinates
      * Using GCPs (set of Ground Control Points) to define non-linear
        relationship between pixel/line and X/Y
      * Using Geolocation Array - full grids of X/Y coordinates for
        each pixel of a raster

    The relation between X/Y coordinates of the raster and latitude/longitude
    coordinates is defined by projection type and projection parameters.
    These pieces of information are therefore stored in Domain:

      * Type and parameters of projection +

        * GeoTransform, or
        * GCPs, or
        * GeolocationArrays

    Domain has methods for basic operations with georeference information:

      * creating georeference from input options;
      * fetching corner, border or full grids of X/Y coordinates;
      * making map of the georeferenced grid in a PNG or KML file;
      * and some more...

    The main attribute of Domain is a VRT object self.vrt.
    Nansat inherits from Domain and adds bands to self.vrt

    Raises
    -------
    NansatProjectionError : occurs when Projection() is empty
        despite it is required for creating extentDic.
    OptionError : occures when the arguments are not proper.

    See Also
    ---------
    Nansat.reproject()

    `<http://www.gdal.org/gdalwarp.html>`_

    `<http://trac.osgeo.org/proj/>`_

    `<http://spatialreference.org/>`_

    `<http://www.gdal.org/osr_tutorial.html>`_

    """

    OUTPUT_SEPARATOR = '-' * 40 + '\n'
    KML_BASE = '''<?xml version="1.0" encoding="UTF-8"?>
    <kml xmlns="http://www.opengis.net/kml/2.2"
    xmlns:gx="http://www.google.com/kml/ext/2.2"
    xmlns:kml="http://www.opengis.net/kml/2.2"
    xmlns:atom="http://www.w3.org/2005/Atom">
    {content}
    </kml>'''

    # instance attributes
    vrt = None
    logger = None
    name = None

    def __init__(self, srs=None, ext=None, ds=None, **kwargs):
        """Create Domain from GDALDataset or string options or lat/lon grids"""
        # If too much information is given raise error
        if ds is not None and srs is not None and ext is not None:
            raise ValueError('Ambiguous specification of both dataset, srs- and ext-strings.')

        # choose between input opitons:
        # ds
        # ds and srs
        # srs and ext

        # if only a dataset is given:
        #     copy geo-reference from the dataset
        if ds is not None and srs is None:
            self.vrt = VRT.from_gdal_dataset(ds)

        # If dataset and srs are given (but not ext):
        #   use AutoCreateWarpedVRT to determine bounds and resolution
        elif ds is not None and srs is not None:
            srs = NSR(srs)
            tmp_vrt = gdal.AutoCreateWarpedVRT(ds, None, srs.wkt)
            if tmp_vrt is None:
                raise NansatProjectionError('Could not warp the given dataset to the given SRS.')
            else:
                self.vrt = VRT.from_gdal_dataset(tmp_vrt)

        # If SpatialRef and extent string are given (but not dataset)
        elif srs is not None and ext is not None:
            srs = NSR(srs)
            # create full dictionary of parameters
            extent_dict = Domain._create_extent_dict(ext)

            # convert -lle to -te
            if 'lle' in extent_dict.keys():
                extent_dict = self._convert_extentDic(srs, extent_dict)

            # get size/extent from the created extent dictionary
            geo_transform, raster_x_size, raster_y_size = self._get_geotransform(extent_dict)
            # create VRT object with given geo-reference parameters
            self.vrt = VRT.from_dataset_params(x_size=raster_x_size, y_size=raster_y_size,
                                               geo_transform=geo_transform,
                                               projection=srs.wkt,
                                               gcps=[], gcp_projection='')
        elif 'lat' in kwargs and 'lon' in kwargs:
            warnings.warn('Domain(lon=lon, lat=lat) will be deprectaed!'
                          'Use Domain.from_lonlat()', NansatFutureWarning)
            # create self.vrt from given lat/lon
            self.vrt = VRT.from_lonlat(kwargs['lon'], kwargs['lat'])
        else:
            raise ValueError('"dataset" or "srsString and extentString" '
                              'or "dataset and srsString" are required')

    @classmethod
    def from_lonlat(cls, lon, lat, add_gcps=True):
        """Create Domain object from input longitudes, latitudes arrays

        Parameters
        ----------
        lon : numpy.ndarray
            longitudes
        lat : numpy.ndarray
            latitudes
        add_gcps : bool
            Add GCPs from lon/lat arrays.

        Returns
        -------
            d : Domain

        Examples
        --------
            >>> lon, lat = np.meshgrid(range(10), range(10))
            >>> d1 = Domain.from_lonlat(lon, lat)
            >>> d2 = Domain.from_lonlat(lon, lat, add_gcps=False) # add only geolocation arrays

        """
        d = cls.__new__(cls)
        d.vrt = VRT.from_lonlat(lon, lat, add_gcps)
        return d

    def __repr__(self):
        """Creates string with basic info about the Domain object

        Returns
        ---------
         out_str : str
            size, projection and corner coordinates

        """
        corners_temp = '\t (%6.2f, %6.2f)  (%6.2f, %6.2f)\n'
        wkt, src = self.vrt.get_projection()
        out_str = 'Domain:[%d x %d]\n' % self.shape()[::-1]
        out_str += self.OUTPUT_SEPARATOR
        corners = self.get_corners()
        out_str += 'Projection(%s):\n' % src
        out_str += (NSR(wkt).ExportToPrettyWkt(1) + '\n')
        out_str += self.OUTPUT_SEPARATOR
        out_str += 'Corners (lon, lat):\n'
        out_str += corners_temp % (corners[0][0], corners[1][0], corners[0][2], corners[1][2])
        out_str += corners_temp % (corners[0][1], corners[1][1], corners[0][3], corners[1][3])
        return out_str

    def write_kml(self, xmlFileName=None, kmlFileName=None):
        """Write KML file with domains

        Convert XML-file with domains into KML-file for GoogleEarth
        or write KML-file with the current Domain

        Parameters
        -----------
        xmlFileName : string, optional
            Name of the XML-file to convert. If only this value is given
            - kmlFileName=xmlFileName+'.kml'

        kmlFileName : string, optional
            Name of the KML-file to generate from the current Domain

        """
        xml_filename = xmlFileName
        kml_filename = kmlFileName

        template = '''<Document>
        \t<name>{filename}</name>
        \t\t<Folder><name>{filename}</name><open>1</open>
        {borders}
        \t\t</Folder></Document>'''

        # test input options
        if xml_filename and not kml_filename:
            # if only input XML-file is given - convert it to KML

            # open XML, get all domains
            with open(xml_filename, 'rb') as xml_file:
                xml_domains = list(ElementTree(file=xml_file).getroot())

            # convert domains in XML into list of domains
            domains = [Domain(srs=xml_filename, ext=domain.attrib['name'])
                       for domain in xml_domains]

        elif not xml_filename and kml_filename:
            # if only output KML-file is given
            # then convert the current domain to KML
            domains = [self]

        else:
            # otherwise it is potentially error
            raise ValueError('Either xmlFileName(%s)\
             or kmlFileName(%s) are wrong' % (xml_filename, kml_filename))

        # get border of each domain and join them to a one string
        borders = ''.join([domain._get_border_kml() for domain in domains])
        # open KML, write the modified template
        with open(kml_filename, 'wt') as kml_file:
            kml_content = template.format(name=self.name, filename=kml_filename, borders=borders)
            kml_file.write(self.KML_BASE.format(content=kml_content))

    def _get_border_kml(self):
        """Generate Placemark entry for KML

        Returns
        --------
        kmlEntry : str
            String with the Placemark entry

        """
        klm_entry = '''\t\t\t<Placemark>
        \t\t\t\t<name>{name}</name>
        \t\t\t\t<Style>
        \t\t\t\t\t<LineStyle><color>ffffffff</color></LineStyle>
        \t\t\t\t\t<PolyStyle><fill>0</fill></PolyStyle>
        \t\t\t\t</Style>
        \t\t\t\t<Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing><coordinates>
        {coordinates}
        </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>'''
        domain_lon, domain_lat = self.get_border()
        # convert Border coordinates into KML-like string
        coordinates = ''.join(['%f,%f,0 ' % (lon, lat) for lon, lat in zip(domain_lon, domain_lat)])
        return klm_entry.format(name=self.name, coordinates=coordinates)

    def write_kml_image(self, kmlFileName, kmlFigureName=None):
        """Create KML file for already projected image

        Write Domain Image into KML-file for GoogleEarth

        Parameters
        -----------
        kmlFileName : str
            Name of the KML-file to generate from the current Domain
        kmlFigureName : str
            Name of the projected image stored in .png format

        Examples
        ---------
            >>> n.undo(100) # cancel previous reprojection
            >>> lons, lats = n.get_corners() # Get corners of the image and the pixel resolution
            >>> srsString = '+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs'
            >>> extentString = '-lle %f %f %f %f -ts 3000 3000'
                % (min(lons), min(lats), max(lons), max(lats))
            >>> d = Domain(srs=srsString, ext=extentString) # Create Domain with
                stereographic projection, corner coordinates and resolution 1000m
            >>> n.reproject(d)
            >>> n.write_figure(filename=figureName, bands=[3], clim=[0,0.15],
                                cmapName='gray', transparency=0)
            >>> n.write_kml_image(kmlFileName=oPath + filename + '.kml',
                                kmlFigureName=figureName) # 6.

        """

        kml_filename = kmlFileName
        kml_figurename = kmlFigureName

        template = '''<GroundOverlay>
        \t<name>{filename}</name>
        \t<Icon>
        \t\t<href>{figurename}</href>
        \t\t<viewBoundScale>0.75</viewBoundScale>
        \t</Icon>
        \t<LatLonBox>
        \t\t<north>{north}</north>
        \t\t<south>{south}</south>
        \t\t<east>{east}</east>
        \t\t<west>{west}</west>
        \t</LatLonBox>
        </GroundOverlay>'''
        # test input options
        if kml_figurename is None:
            raise ValueError('kmlFigureName(%s) is not specified' % kmlFigureName)

        # get corner of the domain and add to KML
        domain_lon, domain_lat = self.get_corners()
        with open(kml_filename, 'wt') as kml_file:
            kml_content = template.format(filename=kml_filename, figurename=kml_figurename,
                                          north=max(domain_lat), south=min(domain_lat),
                                          east=max(domain_lon), west=min(domain_lon))
            kml_file.write(self.KML_BASE.format(content=kml_content))

    def get_geolocation_grids(self, stepSize=1, dst_srs=None):
        """Get longitude and latitude grids representing the full data grid

        If GEOLOCATION is not present in the self.vrt.dataset then grids
        are generated by converting pixel/line of each pixel into lat/lon
        If GEOLOCATION is present in the self.vrt.dataset then grids are read
        from the geolocation bands.

        Parameters
        -----------
        stepSize : int
            Reduction factor if output is desired on a reduced grid size

        Returns
        --------
        longitude : numpy array
            grid with longitudes
        latitude : numpy array
            grid with latitudes
        """
        if dst_srs is None:
            dst_srs = NSR()
        step_size = stepSize
        x_vec = list(range(0, self.vrt.dataset.RasterXSize, step_size))
        y_vec = list(range(0, self.vrt.dataset.RasterYSize, step_size))
        x_grid, y_grid = np.meshgrid(x_vec, y_vec)

        if self.vrt.geolocation is not None and len(self.vrt.geolocation.data) > 0:
            # if the vrt dataset has geolocationArray
            # read lon,lat grids from geolocationArray
            lon_grid, lat_grid = self.vrt.geolocation.get_geolocation_grids()
            lon_arr, lat_arr = lon_grid[y_grid, x_grid], lat_grid[y_grid, x_grid]
        else:
            # generate lon,lat grids using GDAL Transformer
            lon_vec, lat_vec = self.transform_points(x_grid.flatten(), y_grid.flatten(),
                                                     dst_srs=dst_srs)
            lon_arr = lon_vec.reshape(x_grid.shape)
            lat_arr = lat_vec.reshape(x_grid.shape)

        return lon_arr, lat_arr

    def _convert_extentDic(self, dst_srs, extentDic):
        """Convert -lle option (lat/lon) to -te (proper coordinate system)

        Source SRS from LAT/LON projection and target SRS from dstWKT.
        Create osr.CoordinateTransformation based on these SRSs and
        convert given values in degrees to the destination coordinate
        system given by WKT.
        Add key 'te' and the converted values into the extentDic.

        Parameters
        -----------
        dst_srs : NSR
            Destination Spatial Reference
        extentDic : dictionary
            dictionary with 'lle' key

        Returns
        --------
        extentDic : dictionary
            input dictionary + 'te' key and its values

        """
        coorTrans = osr.CoordinateTransformation(NSR(), dst_srs)

        # convert lat/lon given by 'lle' to the target coordinate system and
        # add key 'te' and the converted values to extentDic
        x1, y1, _ = coorTrans.TransformPoint(extentDic['lle'][0], extentDic['lle'][3])
        x2, y2, _ = coorTrans.TransformPoint(extentDic['lle'][2], extentDic['lle'][3])
        x3, y3, _ = coorTrans.TransformPoint(extentDic['lle'][2], extentDic['lle'][1])
        x4, y4, _ = coorTrans.TransformPoint(extentDic['lle'][0], extentDic['lle'][1])

        minX = min([x1, x2, x3, x4])
        maxX = max([x1, x2, x3, x4])
        minY = min([y1, y2, y3, y4])
        maxY = max([y1, y2, y3, y4])

        extentDic['te'] = [minX, minY, maxX, maxY]

        return extentDic

    @staticmethod
    def _add_to_dict(extent, option):
        """Convert options to list of float values and add to <extent> dict"""
        try:
            parameters = [float(el.strip()) for el in option[1:]]
        except ValueError:
            raise ValueError('Input values must be int or float')

        key = option[0].strip().replace('-', '')
        extent[key] = parameters
        return key, extent

    @staticmethod
    def _validate_ts_tr(options):
        example = '<-tr x_resolution y_resolution> or <-ts width height>'
        Domain._check_size(len(options), 2, ('-ts', '-tr'), example)
        if options[0] <= 0 or options[1] <= 0:
            raise ValueError('Resolution or width and height must be bigger than 0: %s' % example)

    @staticmethod
    def _validate_te_lle(options):
        example = '<-te x_min y_min x_max y_max> or <-lle min_lon min_lat max_lon max_lat>'
        Domain._check_size(len(options), 4, ('-te', '-lle'), example)
        if options[0] >= options[2] or options[1] >= options[3]:
            raise ValueError('Min cannot be bigger than max: %s' % example)

    @staticmethod
    def _check_size(params_len, size, names, example):
        if params_len != size:
            raise ValueError('%s and %s requires exactly %s parameters (%s given): %s'
                              % (names[0], names[1], size, params_len, example))

    @staticmethod
    def _gen_regexp(param_1, param_2, size):
        return '(-%s|-%s)%s\s?' % (param_1, param_2, '(\s+[-+]?\d*[.\d*]*)' * size)

    @staticmethod
    def _create_extent_dict(extent_str):
        """Create a dictionary from extentString

        Check if extentString is proper.
            * '-te' and '-lle' take 4 numbers.
            * '-ts' and '-tr' take 2 numbers.
            * the combination should be ('-te' or '-lle') and ('-ts' or '-tr')
        If it is proper, create a dictionary
        Otherwise, raise the error.

        Parameters
        -----------
        extentString : string
            '-te xMin yMin xMax yMax',
            '-tr xResolution yResolution',
            '-ts width height',
            '-lle minlon minlat maxlon maxlat'

        Returns
        --------
        extentDict : dictionary
            has key ('te' or 'lle') and ('tr' or 'ts') and their values.

        Raises
        -------
        ValueError : occurs when the extent_str is improper

        """

        combinations = [('te', 'lle', 4), ('ts', 'tr', 2)]
        extent_dict = {}
        for combination in combinations:
            try:
                option = re.findall(Domain._gen_regexp(*combination), extent_str)[0]
            except IndexError:
                raise ValueError('<extent_dict> must contains exactly 2 parameters '
                                  '("-te" or "-lle") and ("-ts" or "-tr")')
            key, extent_dict = Domain._add_to_dict(extent_dict, option)
            if key is 'te' or key is 'lle':
                Domain._validate_te_lle(extent_dict[key])
            elif key is 'ts' or key is 'tr':
                Domain._validate_ts_tr(extent_dict[key])

        return extent_dict

    def get_border(self, n_points=10, fix_lon=True, **kwargs):
        """Generate two vectors with values of lat/lon for the border of domain

        Parameters
        -----------
        n_points : int, optional
            Number of points on each border
        fix_lon : bool
            Convert longitudes to positive numbers when Domain crosses dateline?

        Returns
        --------
        lonVec, latVec : lists
            vectors with lon/lat values for each point at the border

        """
        x_size, y_size = self.shape()[::-1]
        x_rc_vec = Domain._get_row_col_vector(x_size, n_points)
        y_rc_vec = Domain._get_row_col_vector(y_size, n_points)
        col_vec, row_vec = Domain._compound_row_col_vectors(x_size, y_size, x_rc_vec, y_rc_vec)
        lon, lat = self.transform_points(col_vec, row_vec)
        lon = lon.round(decimals=4)
        lat = lat.round(decimals=4)

        crosses_dateline = np.diff(lon).max() > 100 # arbitrary number, larger than a maximum step
        if crosses_dateline and fix_lon:
            lon[lon < 0] += 360

        return lon, lat

    @staticmethod
    def _compound_row_col_vectors(x_size, y_size, x_vec, y_vec):
        col_vec = (x_vec + [x_size] * len(y_vec) + x_vec[::-1] + [0] * len(y_vec))
        row_vec = ([0] * len(x_vec) + y_vec + [y_size] * len(x_vec) + y_vec[::-1])
        return col_vec, row_vec

    @staticmethod
    def _get_row_col_vector(raster_size, n_points):
        # Int because py3 division returns float
        step = max(1, int(raster_size / n_points))
        rc_vec = list(range(0, raster_size, step))[0:n_points]
        rc_vec.append(raster_size)
        return rc_vec

    def get_border_wkt(self, *args, **kwargs):
        """Creates string with WKT representation of the border polygon

        Returns
        --------
        WKTPolygon : string
            string with WKT representation of the border polygon

        """
        lon_vec, lat_vec = self.get_border(*args, **kwargs)

        ''' The following causes erratic geometry when using
        WKTReader().read(n.get_border_wkt(n_points=1000)) - only commented out
        now since this may cause other problems...
        '''
        warnings.warn("> 180 deg correction to longitudes - disabled..")
        polygon_border = ','.join('%s %s' % (lon, lat) for lon, lat in zip(lon_vec, lat_vec))
        # outer quotes have to be double and inner - single!
        # wktPolygon = "PolygonFromText('POLYGON((%s))')" % polyCont
        wkt = 'POLYGON((%s))' % polygon_border
        return wkt

    def get_border_geometry(self, *args, **kwargs):
        """ Get OGR Geometry of the border Polygon

        Returns
        -------
        OGR Geometry : Polygon

        """

        return ogr.CreateGeometryFromWkt(self.get_border_wkt(*args, **kwargs))

    def get_border_geojson(self, *args, **kwargs):
        """Create border of the Polygon in GeoJson format

        Returns
        -------
        the Polygon border in GeoJson format : str

        """
        return ogr.CreateGeometryFromWkt(self.get_border_wkt(*args, **kwargs)).ExportToJson()

    def overlaps(self, anotherDomain):
        """ Checks if this Domain overlaps another Domain

        Returns
        -------
        overlaps : bool
            True if Domains overlaps, False otherwise

        """

        return self.get_border_geometry().Overlaps(anotherDomain.get_border_geometry())

    def intersects(self, anotherDomain):
        """ Checks if this Domain intersects another Domain

        Returns
        -------
        intersects : bool
            True if Domains intersects, False otherwise

        """

        return self.get_border_geometry().Intersects(anotherDomain.get_border_geometry())

    def contains(self, anotherDomain):
        """ Checks if this Domain fully covers another Domain

        Returns
        -------
        contains : bool
            True if this Domain fully covers another Domain, False otherwise

        """

        return self.get_border_geometry().Contains(anotherDomain.get_border_geometry())

    def get_border_postgis(self, **kwargs):
        """ Get PostGIS formatted string of the border Polygon

        Returns
        -------
        'PolygonFromText(PolygonWKT)' : str

        """

        return "PolygonFromText('%s')" % self.get_border_wkt(**kwargs)

    def get_corners(self):
        """Get coordinates of corners of the Domain

        Returns
        --------
        lonVec, latVec : lists
            vectors with lon/lat values for each corner

        """

        col_vec = [0, 0, self.vrt.dataset.RasterXSize, self.vrt.dataset.RasterXSize]
        row_vec = [0, self.vrt.dataset.RasterYSize, 0, self.vrt.dataset.RasterYSize]
        return self.transform_points(col_vec, row_vec)

    def get_min_max_lon_lat(self):
        """Get minimum and maximum of longitude and latitude geolocation grids

        Returns
        --------
        min_lon, max_lon, min_lat, max_lat, : float
            min/max lon/lat values for the Domain

        """
        lon_grd, lat_grd = self.get_geolocation_grids()
        return lon_grd.min(), lon_grd.max(), lat_grd.min(), lat_grd.max(),

    def get_pixelsize_meters(self):
        """Returns the pixelsize (deltaX, deltaY) of the domain

        For projected domains, the exact result which is constant
        over the domain is returned.
        For geographic (lon-lat) projections, or domains with no geotransform,
        the haversine formula is used to calculate the pixel size
        in the center of the domain.

        Returns
        --------
        delta_x, delta_y : float
            pixel size in X and Y directions given in meters
        """

        srs = osr.SpatialReference(self.vrt.dataset.GetProjection())
        if srs.IsProjected:
            if srs.GetAttrValue('unit') == 'metre':
                geoTransform = self.vrt.dataset.GetGeoTransform()
                delta_x = abs(geoTransform[1])
                delta_y = abs(geoTransform[5])
                return delta_x, delta_y

        # Estimate pixel size in center of domain using haversine formula
        center_col = round(self.vrt.dataset.RasterXSize/2)
        center_row = round(self.vrt.dataset.RasterYSize/2)
        lon00, lat00 = self.transform_points([center_col], [center_row])
        lon01, lat01 = self.transform_points([center_col], [center_row + 1])
        lon10, lat10 = self.transform_points([center_col + 1], [center_row])

        delta_x = haversine(lon00, lat00, lon01, lat01)
        delta_y = haversine(lon00, lat00, lon10, lat10)
        return delta_x[0], delta_y[0]

    @staticmethod
    def _get_geotransform(extent_dict):
        """Get the new coordinates and raster size are calculated based on the given extentDic.

        Parameters
        -----------
        extent_dict : dictionary
            includes 'te' key and 'ts' or 'tr' key

        Raises
        -------
        ValueError : occurs when maxX - minX < 0 or maxY - minY < 0

        Returns
        --------
        coordinate : list with 6 float
            GeoTransform

        raster_x_size and raster_y_size

        """
        width = extent_dict['te'][2] - extent_dict['te'][0]
        height = extent_dict['te'][3] - extent_dict['te'][1]

        if 'tr' in extent_dict.keys():
            resolution_x, resolution_y, raster_x_size, raster_y_size = \
                Domain._transform_tr(width, height, extent_dict['tr'])
        else:
            resolution_x, resolution_y, raster_x_size, raster_y_size = \
                Domain._transform_ts(width, height, extent_dict['ts'])

        # create a list for GeoTransform
        coordinates = [extent_dict['te'][0], resolution_x, 0.0,
                       extent_dict['te'][3], 0.0, resolution_y]

        return coordinates, int(raster_x_size), int(raster_y_size)

    @staticmethod
    def _transform_tr(width, height, tr_arr):
        """Calculate X and Y resolution and raster sizes from the "-tr" parameter

        Parameters
        -----------
        width : float, width of domain calculated from the "-te" extent parameter
        height: float, height of domain calculated from the "-te" extent parameter
        tr_arr: list, [<x_resolution>, <y_resolution>]

        Raises
        -------
        ValueError : occurs when the given resolution is larger than width or height.

        Returns
        --------
        resolution_x, resolution_y, raster_x_size, raster_y_size : float

        """
        resolution_x = tr_arr[0]
        resolution_y = -(tr_arr[1])

        if width < resolution_x or height < resolution_y:
            raise ValueError('"-tr" is too large. width is %s, height is %s ' % (width, height))

        raster_x_size = width / resolution_x
        raster_y_size = abs(height / resolution_y)

        return resolution_x, resolution_y, raster_x_size, raster_y_size

    @staticmethod
    def _transform_ts(width, height, ts_arr):
        raster_x_size, raster_y_size = ts_arr
        resolution_x = width / raster_x_size
        resolution_y = -abs(height / raster_y_size)

        return resolution_x, resolution_y, raster_x_size, raster_y_size

    def transform_points(self, colVector, rowVector, DstToSrc=0, dst_srs=None):
        """Transform given lists of X,Y coordinates into lon/lat or inverse

        Parameters
        -----------
        colVector : lists
            X and Y coordinates in pixel/line or lon/lat  coordinate system
        DstToSrc : 0 or 1

            - 0 - forward transform (pix/line => lon/lat)
            - 1 - inverse transformation

        dst_srs : NSR
            destination spatial reference

        Returns
        --------
        X, Y : lists
            X and Y coordinates in lon/lat or pixel/line coordinate system

        """
        if dst_srs is None:
            dst_srs = NSR()
        return self.vrt.transform_points(colVector, rowVector, dst2src=DstToSrc, dst_srs=dst_srs)

    def azimuth_y(self, reductionFactor=1):
        """Calculate the angle of each pixel position vector with respect to
        the Y-axis (azimuth).

        In general, azimuth is the angle from a reference vector (e.g., the
        direction to North) to the chosen position vector. The azimuth
        increases clockwise from direction to North.
        http://en.wikipedia.org/wiki/Azimuth

        Parameters
        -----------
        reductionFactor : integer
            factor by which the size of the output array is reduced

        Returns
        -------
        azimuth : numpy array
            Values of azimuth in degrees in range 0 - 360

        """

        lon_grd, lat_grd = self.get_geolocation_grids(reductionFactor)
        a = initial_bearing(lon_grd[1:, :], lat_grd[1:, :], lon_grd[:-1:, :], lat_grd[:-1:, :])
        # Repeat last row once to match size of lon-lat grids
        a = np.vstack((a, a[-1, :]))
        return a

    def shape(self):
        """Return Numpy-like shape of Domain object (ySize, xSize)

        Returns
        --------
        shape : tuple of two INT
            Numpy-like shape of Domain object (ySize, xSize)

        """
        return self.vrt.dataset.RasterYSize, self.vrt.dataset.RasterXSize

    def reproject_gcps(self, srs_string=''):
        """Reproject all GCPs to a new spatial reference system

        Necessary before warping an image if the given GCPs
        are in a coordinate system which has a singularity
        in (or near) the destination area (e.g. poles for lonlat GCPs)

        Parameters
        ----------
        srs_string : string
            SRS given as Proj4 string. If empty '+proj=stere' is used

        Notes
        --------
            Reprojects all GCPs to new SRS and updates GCPProjection
        """
        if srs_string == '':
            lon, lat = self.get_border()
            srs_string = '+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=%f +lon_0=%f +no_defs' \
                         % (np.nanmedian(lat), np.nanmedian(lon))
        self.vrt.reproject_gcps(srs_string)
