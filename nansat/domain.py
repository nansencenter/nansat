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
from __future__ import absolute_import
import re
from math import sin, pi, cos, acos, copysign
import string
from xml.etree.ElementTree import ElementTree

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon

from nansat.tools import add_logger, initial_bearing, haversine, gdal, osr, ogr
from nansat.tools import OptionError, ProjectionError
from nansat.nsr import NSR
from nansat.vrt import VRT


class Domain(object):
    '''Container for geographical reference of a raster

    A Domain object describes all attributes of geographical
    reference of a raster:
      * width and height (number of pixels)
      * pixel size (e.g. in decimal degrees or in meters)
      * relation between pixel/line coordinates and geographical
        coordinates (e.g. a linear relation)
      * type of data projection (e.g. geographical or stereographic)

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

    '''
    def __init__(self, srs=None, ext=None, ds=None, lon=None,
                 lat=None, name='', logLevel=None):
        '''Create Domain from GDALDataset or string options or lat/lon grids

        d = Domain(srs, ext)
            Size, extent and spatial reference is given by strings
        d = Domain(ds=GDALDataset):
            Size, extent and spatial reference is copied from input
            GDAL dataset
        d = Domain(srs, ds=GDALDataset):
            Spatial reference is given by srs, but size and extent is
            determined
            from input GDAL dataset
        d = Domain(lon=lonGrid, lat=latGrid)
            Size, extent and spatial reference is given by two grids

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
        lat : Numpy array
            Grid with latitudes
        lon : Numpy array
            Grid with longitudes
        name : string, optional
            Name to be added to the Domain object
        logLevel : int, optional, default=30
            level of logging

        Raises
        -------
        ProjectionError : occurs when Projection() is empty
            despite it is required for creating extentDic.
        OptionError : occures when the arguments are not proper.

        Modifies
        ---------
        self.vrt.datasetset : dataset in memory
            dataset is created based on the input arguments

        See Also
        ---------
        Nansat.reproject()
        [http://www.gdal.org/gdalwarp.html]
        [http://trac.osgeo.org/proj/]
        [http://spatialreference.org/]
        [http://www.gdal.org/ogr/osr_tutorial.html]

        '''
        # set default attributes
        self.logger = add_logger('Nansat', logLevel)
        self.name = name

        self.logger.debug('ds: %s' % str(ds))
        self.logger.debug('srs: %s' % srs)
        self.logger.debug('ext: %s' % ext)

        # If too much information is given raise error
        if ds is not None and srs is not None and ext is not None:
            raise OptionError('Ambiguous specification of both '
                              'dataset, srs- and ext-strings.')

        # choose between input opitons:
        # ds
        # ds and srs
        # srs and ext
        # lon and lat

        # if only a dataset is given:
        #     copy geo-reference from the dataset
        if ds is not None and srs is None:
            self.vrt = VRT(gdalDataset=ds)

        # If dataset and srs are given (but not ext):
        #   use AutoCreateWarpedVRT to determine bounds and resolution
        elif ds is not None and srs is not None:
            srs = NSR(srs)
            tmpVRT = gdal.AutoCreateWarpedVRT(ds, None, srs.wkt)
            if tmpVRT is None:
                raise ProjectionError('Could not warp the given dataset'
                                      'to the given SRS.')
            else:
                self.vrt = VRT(gdalDataset=tmpVRT)

        # If SpatialRef and extent string are given (but not dataset)
        elif srs is not None and ext is not None:
            srs = NSR(srs)
            # create full dictionary of parameters
            extentDic = self._create_extentDic(ext)

            # convert -lle to -te
            if 'lle' in extentDic.keys():
                extentDic = self._convert_extentDic(srs, extentDic)

            # get size/extent from the created extet dictionary
            [geoTransform,
             rasterXSize, rasterYSize] = self._get_geotransform(extentDic)
            # create VRT object with given geo-reference parameters
            self.vrt = VRT(srcGeoTransform=geoTransform,
                           srcProjection=srs.wkt,
                           srcRasterXSize=rasterXSize,
                           srcRasterYSize=rasterYSize)
            self.extentDic = extentDic
        elif lat is not None and lon is not None:
            # create self.vrt from given lat/lon
            self.vrt = VRT(lat=lat, lon=lon)
        else:
            raise OptionError('"dataset" or "srsString and extentString" '
                              'or "dataset and srsString" are required')

        self.logger.debug('vrt.dataset: %s' % str(self.vrt.dataset))

    def __repr__(self):
        '''Creates string with basic info about the Domain object

        Modifies
        ---------
        Print size, projection and corner coordinates

        '''
        outStr = 'Domain:[%d x %d]\n' % (self.vrt.dataset.RasterXSize,
                                         self.vrt.dataset.RasterYSize)
        outStr += '-' * 40 + '\n'
        try:
            corners = self.get_corners()
        except:
            self.logger.error('Cannot read projection from source!')
        else:
            outStr += 'Projection:\n'
            outStr += (NSR(self.vrt.get_projection()).ExportToPrettyWkt(1) +
                       '\n')
            outStr += '-' * 40 + '\n'
            outStr += 'Corners (lon, lat):\n'
            outStr += '\t (%6.2f, %6.2f)  (%6.2f, %6.2f)\n' % (corners[0][0],
                                                               corners[1][0],
                                                               corners[0][2],
                                                               corners[1][2])
            outStr += '\t (%6.2f, %6.2f)  (%6.2f, %6.2f)\n' % (corners[0][1],
                                                               corners[1][1],
                                                               corners[0][3],
                                                               corners[1][3])
        return outStr

    def write_kml(self, xmlFileName=None, kmlFileName=None):
        '''Write KML file with domains

        Convert XML-file with domains into KML-file for GoogleEarth
        or write KML-file with the current Domain

        Parameters
        -----------
        xmlFileName : string, optional
            Name of the XML-file to convert. If only this value is given
            - kmlFileName=xmlFileName+'.kml'

        kmlFileName : string, optional
            Name of the KML-file to generate from the current Domain

        '''
        # test input options
        if xmlFileName is not None and kmlFileName is None:
            # if only input XML-file is given - convert it to KML

            # open XML, get all domains
            xmlFile = file(xmlFileName, 'rb')
            kmlFileName = xmlFileName + '.kml'
            xmlDomains = ElementTree(file=xmlFile).getroot()
            xmlFile.close()

            # convert domains in XML into list of domains
            domains = []
            for xmlDomain in list(xmlDomains):
                # append Domain object to domains list
                domainName = xmlDomain.attrib['name']
                domains.append(Domain(srs=xmlFileName, ext=domainName))

        elif xmlFileName is None and kmlFileName is not None:
            # if only output KML-file is given
            # then convert the current domain to KML
            domains = [self]

        else:
            # otherwise it is potentially error
            raise OptionError('Either xmlFileName(%s)\
             or kmlFileName(%s) are wrong' % (xmlFileName, kmlFileName))

        # open KML, write header
        kmlFile = file(kmlFileName, 'wt')
        kmlFile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        kmlFile.write('<kml xmlns="http://www.opengis.net/kml/2.2" '
                      'xmlns:gx="http://www.google.com/kml/ext/2.2" '
                      'xmlns:kml="http://www.opengis.net/kml/2.2" '
                      'xmlns:atom="http://www.w3.org/2005/Atom">\n')
        kmlFile.write('<Document>\n')
        kmlFile.write('    <name>%s</name>\n' % kmlFileName)
        kmlFile.write('        <Folder><name>%s</name><open>1</open>\n'
                      % kmlFileName)

        # get border of each domain and add to KML
        for domain in list(domains):
            kmlEntry = domain._get_border_kml()
            kmlFile.write(kmlEntry)

        # write footer and close
        kmlFile.write('        </Folder></Document></kml>\n')
        kmlFile.close()

    def write_kml_image(self, kmlFileName=None, kmlFigureName=None):
        '''Create KML file for already projected image

        Write Domain Image into KML-file for GoogleEarth

        Parameters
        -----------
        kmlFileName : string, optional
            Name of the KML-file to generate from the current Domain
        kmlFigureName : string, optional
            Name of the projected image stored in .png format

        Examples
        ---------
        # First of all, reproject an image into Lat/Lon WGS84
          (Simple Cylindrical) projection
        # 1. Cancel previous reprojection
        # 2. Get corners of the image and the pixel resolution
        # 3. Create Domain with stereographic projection,
        #    corner coordinates and resolution 1000m
        # 4. Reproject
        # 5. Write image
        # 6. Write KML for the image
        n.reproject() # 1.
        lons, lats = n.get_corners() # 2.
        srsString = '+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs'
        extentString = '-lle %f %f %f %f -ts 3000 3000'
        % (min(lons), min(lats), max(lons), max(lats))
        d = Domain(srs=srsString, ext=extentString) # 3.
        n.reproject(d) # 4.
        n.write_figure(fileName=figureName, bands=[3], clim=[0,0.15],
                       cmapName='gray', transparency=0) # 5.
        n.write_kml_image(kmlFileName=oPath + fileName + '.kml',
                          kmlFigureName=figureName) # 6.

        '''
        # test input options
        if kmlFileName is None:
            raise OptionError('kmlFileName(%s) is wrong' % (kmlFileName))

        if kmlFigureName is None:
            raise OptionError('kmlFigureName(%s) is not specified'
                              % (kmlFigureName))

        # open KML, write header
        kmlFile = file(kmlFileName, 'wt')
        kmlFile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        kmlFile.write('<kml xmlns="http://www.opengis.net/kml/2.2" '
                      'xmlns:gx="http://www.google.com/kml/ext/2.2" '
                      'xmlns:kml="http://www.opengis.net/kml/2.2" '
                      'xmlns:atom="http://www.w3.org/2005/Atom">\n')
        kmlFile.write('<GroundOverlay>\n')
        kmlFile.write('    <name>%s</name>\n' % kmlFileName)
        kmlFile.write('    <Icon>\n')
        kmlFile.write('        <href>%s</href>\n' % kmlFigureName)
        kmlFile.write('        <viewBoundScale>0.75</viewBoundScale>\n')
        kmlFile.write('    </Icon>\n')

        # get corner of the domain and add to KML
        domainLon, domainLat = self.get_corners()

        kmlFile.write('    <LatLonBox>\n')
        kmlFile.write('        <north>%s</north>\n' % max(domainLat))
        kmlFile.write('        <south>%s</south>\n' % min(domainLat))
        kmlFile.write('        <east>%s</east>\n' % max(domainLon))
        kmlFile.write('        <west>%s</west>\n' % min(domainLon))
        kmlFile.write('    </LatLonBox>\n')

        # write footer and close
        kmlFile.write('</GroundOverlay>\n')
        kmlFile.write('</kml>')
        kmlFile.close()

    def get_geolocation_grids(self, stepSize=1, dstSRS=NSR()):
        '''Get longitude and latitude grids representing the full data grid

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
        '''

        X = range(0, self.vrt.dataset.RasterXSize, stepSize)
        Y = range(0, self.vrt.dataset.RasterYSize, stepSize)
        Xm, Ym = np.meshgrid(X, Y)

        if len(self.vrt.geolocationArray.d) > 0:
            # if the vrt dataset has geolocationArray
            # read lon,lat grids from geolocationArray
            lon, lat = self.vrt.geolocationArray.get_geolocation_grids()
            longitude, latitude = lon[Ym, Xm], lat[Ym, Xm]
        else:
            # generate lon,lat grids using GDAL Transformer
            lonVec, latVec = self.transform_points(Xm.flatten(), Ym.flatten(), dstSRS=dstSRS)
            longitude = lonVec.reshape(Xm.shape)
            latitude = latVec.reshape(Xm.shape)

        return longitude, latitude

    def _convert_extentDic(self, dstSRS, extentDic):
        '''Convert -lle option (lat/lon) to -te (proper coordinate system)

        Source SRS from LAT/LON projection and target SRS from dstWKT.
        Create osr.CoordinateTransformation based on these SRSs and
        convert given values in degrees to the destination coordinate
        system given by WKT.
        Add key 'te' and the converted values into the extentDic.

        Parameters
        -----------
        dstSRS : NSR
            Destination Spatial Reference
        extentDic : dictionary
            dictionary with 'lle' key

        Returns
        --------
        extentDic : dictionary
            input dictionary + 'te' key and its values

        '''
        coorTrans = osr.CoordinateTransformation(NSR(), dstSRS)

        # convert lat/lon given by 'lle' to the target coordinate system and
        # add key 'te' and the converted values to extentDic
        x1, y1, _ = coorTrans.TransformPoint(extentDic['lle'][0],
                                             extentDic['lle'][3])
        x2, y2, _ = coorTrans.TransformPoint(extentDic['lle'][2],
                                             extentDic['lle'][3])
        x3, y3, _ = coorTrans.TransformPoint(extentDic['lle'][2],
                                             extentDic['lle'][1])
        x4, y4, _ = coorTrans.TransformPoint(extentDic['lle'][0],
                                             extentDic['lle'][1])

        minX = min([x1, x2, x3, x4])
        maxX = max([x1, x2, x3, x4])
        minY = min([y1, y2, y3, y4])
        maxY = max([y1, y2, y3, y4])

        extentDic['te'] = [minX, minY, maxX, maxY]

        return extentDic

    def _create_extentDic(self, extentString):
        '''Create a dictionary from extentString

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
        extentDic : dictionary
            has key ('te' or 'lle') and ('tr' or 'ts') and their values.

        Raises
        -------
        OptionError : occurs when the extentString is improper

        '''
        extentDic = {}

        # Find -re text
        str_tr = re.findall('-tr\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?',
                            extentString)
        if str_tr != []:
            # Check the number of -tr elements
            elm_str = str(str_tr[0].rstrip())
            elms_str = elm_str.split(None)
            if len(elms_str) != 3 or elms_str[2] == '-':
                raise OptionError('Domain._create_extentDic():'
                                  '-tr is used as'
                                  '"-tr xResolution yResolution"')
            # Add the key and value to extentDic
            extentString = extentString.replace(str_tr[0], '')
            trElem = str(str_tr).split(None)
            trkey = trElem[0].translate(string.maketrans('', ''), "[]-'")
            if trkey != '':
                elements = []
                for i in range(2):
                    elements.append(float(trElem[i + 1].
                                          translate(string.maketrans('', ''),
                                                    "'[]'")))
                extentDic[trkey] = elements

        # Find -ts text
        str_ts = re.findall('-ts\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?',
                            extentString)
        if str_ts != []:
            # Check the number of -ts elements
            elm_str = str(str_ts[0].rstrip())
            elms_str = elm_str.split(None)
            if len(elms_str) != 3 or elms_str[2] == '-':
                raise OptionError('Domain._create_extentDic(): '
                                  '"-ts" is used as "-ts width height"')
            # Add the key and value to extentDic
            extentString = extentString.replace(str_ts[0], '')
            tsElem = str(str_ts).split(None)
            tskey = tsElem[0].translate(string.maketrans('', ''), "[]-'")
            if tskey != '':
                elements = []
                for i in range(2):
                    elements.append(float(tsElem[i + 1].
                                          translate(string.maketrans('', ''),
                                                    "[]'")))
                extentDic[tskey] = elements

        # Find -te text
        str_te = re.findall('-te\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s'
                            '+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?',
                            extentString)
        if str_te != []:
            # Check the number of -te elements
            elm_str = str(str_te[0].rstrip())
            elms_str = elm_str.split(None)
            if len(elms_str) != 5:
                raise OptionError('Domain._create_extentDic():'
                                  '-te is used as "-te xMin yMin xMax yMax"')
            # Add the key and value to extentDic
            extentString = extentString.replace(str_te[0], '')
            teElem = str(str_te).split(None)
            tekey = teElem[0].translate(string.maketrans('', ''), "[]-'")
            if tekey != '':
                elements = []
                for i in range(4):
                    elements.append(float(teElem[i + 1].
                                          translate(string.maketrans('', ''),
                                                    "[]'")))
                extentDic[tekey] = elements

        # Find -lle text
        str_lle = re.findall('-lle\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s'
                             '+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?',
                             extentString)
        if str_lle != []:
            # Check the number of -lle elements
            elm_str = str(str_lle[0].rstrip())
            elms_str = elm_str.split(None)
            if len(elms_str) != 5:
                raise OptionError('Domain._create_extentDic():'
                                  '-lle is used as '
                                  '"-lle minlon minlat maxlon maxlat"')
            # Add the key and value to extentDic
            extentString = extentString.replace(str_lle[0], '')
            lleElem = str(str_lle).split(None)
            llekey = lleElem[0].translate(string.maketrans('', ''), "[]-'")
            if llekey != '':
                elements = []
                for i in range(4):
                    elements.append(float(lleElem[i + 1].
                                          translate(string.maketrans('', ''),
                                                    "[]'")))
                extentDic[llekey] = elements

        result = re.search('\S', extentString)
        # if there are unnecessary letters, give an error
        if result is not None:
            raise OptionError('Domain._create_extentDic():'
                              'extentString is not redable :',
                              extentString)

        # check if one of '-te' and '-lle' is given
        if ('lle' not in extentDic) and ('te' not in extentDic):
            raise OptionError('Domain._create_extentDic():'
                              '"-lle" or "-te" is required.')
        elif ('lle' in extentDic) and ('te' in extentDic):
            raise OptionError('Domain._create_extentDic():'
                              '"-lle" or "-te" should be chosen.')

        # check if one of '-ts' and '-tr' is given
        if ('ts' not in extentDic) and ('tr' not in extentDic):
            raise OptionError('Domain._create_extentDic():'
                              '"-ts" or "-tr" is required.')
        elif ('ts' in extentDic) and ('tr' in extentDic):
            raise OptionError('Domain._create_extentDic():'
                              '"-ts" or "-tr" should be chosen.')
        return extentDic

    def get_border(self, nPoints=10):
        '''Generate two vectors with values of lat/lon for the border of domain

        Parameters
        -----------
        nPoints : int, optional
            Number of points on each border

        Returns
        --------
        lonVec, latVec : lists
            vectors with lon/lat values for each point at the border

        '''
        # prepare vectors with pixels and lines for upper, left, lower
        # and right borders
        sizes = [self.vrt.dataset.RasterXSize, self.vrt.dataset.RasterYSize]

        rcVector1 = [[], []]
        rcVector2 = [[], []]
        # loop for pixels and lines
        for n in range(0, 2):
            step = max(1, sizes[n] / nPoints)
            rcVector1[n] = range(0, sizes[n], step)[0:nPoints]
            rcVector1[n].append(sizes[n])
            rcVector2[n] = rcVector1[n][:]
            rcVector2[n].reverse()

        # coumpund vectors of pixels (col) and lines (row)
        colVector = (rcVector1[0] + [sizes[0]] * len(rcVector1[1]) +
                     rcVector2[0] + [0] * len(rcVector1[1]))
        rowVector = ([0] * len(rcVector1[0]) + rcVector1[1] +
                     [sizes[1]] * len(rcVector1[0]) + rcVector2[1])

        return self.transform_points(colVector, rowVector)

    def _get_border_kml(self, *args, **kwargs):
        '''Generate Placemark entry for KML

        Returns
        --------
        kmlEntry : String
            String with the Placemark entry

        '''
        domainLon, domainLat = self.get_border(*args, **kwargs)

        # convert Border coordinates into KML-like string
        coordinates = ''
        for lon, lat in zip(domainLon, domainLat):
            coordinates += '%f,%f,0 ' % (lon, lat)

        kmlEntry = ''
        # write placemark: name, style, polygon, coordinates
        kmlEntry += '            <Placemark>\n'
        kmlEntry += '                <name>%s</name>\n' % self.name
        kmlEntry += '                <Style>\n'
        kmlEntry += '                    <LineStyle><color>ffffffff</color>'\
                    '</LineStyle>\n'
        kmlEntry += '                    <PolyStyle><fill>0</fill>'\
                    '</PolyStyle>\n'
        kmlEntry += '                </Style>\n'
        kmlEntry += '                <Polygon><tessellate>1</tessellate>'\
                    '<outerBoundaryIs><LinearRing><coordinates>\n'
        kmlEntry += coordinates + '\n'
        kmlEntry += '            </coordinates></LinearRing>'\
                    '</outerBoundaryIs></Polygon></Placemark>\n'

        return kmlEntry

    def get_border_wkt(self, *args, **kwargs):
        '''Creates string with WKT representation of the border polygon

        Returns
        --------
        WKTPolygon : string
            string with WKT representation of the border polygon

        '''
        lonList, latList = self.get_border(*args, **kwargs)

        # apply > 180 deg correction to longitudes
        for ilon, lon in enumerate(lonList):
            lonList[ilon] = copysign(acos(cos(lon * pi / 180.)) / pi * 180,
                                     sin(lon * pi / 180.))

        polyCont = ','.join(str(lon) + ' ' + str(lat)
                            for lon, lat in zip(lonList, latList))
        # outer quotes have to be double and inner - single!
        # wktPolygon = "PolygonFromText('POLYGON((%s))')" % polyCont
        wkt = 'POLYGON((%s))' % polyCont
        return wkt

    def get_border_geometry(self, *args, **kwargs):
        ''' Get OGR Geometry of the border Polygon

        Returns
        -------
        OGR Geometry, type Polygon

        '''

        return ogr.CreateGeometryFromWkt(self.get_border_wkt(*args, **kwargs))

    def overlaps(self, anotherDomain):
        ''' Checks if this Domain overlaps another Domain

        Returns
        -------
        overlaps : bool
            True if Domains overlaps, False otherwise

        '''

        return self.get_border_geometry().Intersects(
                anotherDomain.get_border_geometry())

    def contains(self, anotherDomain):
        ''' Checks if this Domain fully covers another Domain

        Returns
        -------
        contains : bool
            True if this Domain fully covers another Domain, False otherwise

        '''

        return self.get_border_geometry().Contains(
                anotherDomain.get_border_geometry())

    def get_border_postgis(self):
        ''' Get PostGIS formatted string of the border Polygon

        Returns
        -------
        str : 'PolygonFromText(PolygonWKT)'

        '''

        return "PolygonFromText('%s')" % self.get_border_wkt()

    def get_corners(self):
        '''Get coordinates of corners of the Domain

        Returns
        --------
        lonVec, latVec : lists
            vectors with lon/lat values for each corner

        '''

        colVector = [0, 0, self.vrt.dataset.RasterXSize,
                     self.vrt.dataset.RasterXSize]
        rowVector = [0, self.vrt.dataset.RasterYSize, 0,
                     self.vrt.dataset.RasterYSize]
        return self.transform_points(colVector, rowVector)

    def get_min_max_lat_lon(self):
        '''Get minimum and maximum lat and long values in the geolocation grid

        Returns
        --------
        minLat, maxLat, minLon, maxLon : float
            min/max lon/lat values for the Domain

        '''
        allLongitudes, allLatitudes = self.get_geolocation_grids()
        maxLat = -90
        minLat = 90
        for latitudes in allLatitudes:
            for lat in latitudes:
                if lat > maxLat:
                    maxLat = lat
                if lat < minLat:
                    minLat = lat

        maxLon = -180
        minLon = 180
        for longitudes in allLongitudes:
            for lon in longitudes:
                if lon > maxLon:
                    maxLon = lon
                if lon < minLon:
                    minLon = lon

        return minLat, maxLat, minLon, maxLon

    def get_pixelsize_meters(self):
        '''Returns the pixelsize (deltaX, deltaY) of the domain

        For projected domains, the exact result which is constant
        over the domain is returned.
        For geographic (lon-lat) projections, or domains with no geotransform,
        the haversine formula is used to calculate the pixel size
        in the center of the domain.
        Returns
        --------
        deltaX, deltaY : float
        pixel size in X and Y directions given in meters
        '''

        srs = osr.SpatialReference(self.vrt.dataset.GetProjection())
        if srs.IsProjected:
            if srs.GetAttrValue('unit') == 'metre':
                geoTransform = self.vrt.dataset.GetGeoTransform()
                deltaX = abs(geoTransform[1])
                deltaY = abs(geoTransform[5])
                return deltaX, deltaY

        # Estimate pixel size in center of domain using haversine formula
        centerCol = round(self.vrt.dataset.RasterXSize/2)
        centerRow = round(self.vrt.dataset.RasterYSize/2)
        lon00, lat00 = self.transform_points([centerCol], [centerRow])
        lon01, lat01 = self.transform_points([centerCol], [centerRow + 1])
        lon10, lat10 = self.transform_points([centerCol + 1], [centerRow])

        deltaX = haversine(lon00, lat00, lon01, lat01)
        deltaY = haversine(lon00, lat00, lon10, lat10)
        return deltaX[0], deltaY[0]

    def _get_geotransform(self, extentDic):
        '''
        the new coordinates and raster size are calculated based on
        the given extentDic.

        Parameters
        -----------
        extentDic : dictionary
            includes 'te' key and 'ts' or 'tr' key

        Raises
        -------
        OptionError : occurs when maxX - minX < 0 or maxY - minY < 0
        OptionError : occurs when the given resolution is larger than
                     width or height.

        Returns
        --------
        coordinate : list with 6 float
            GeoTransform

        rasterSize : list with two int
            rasterXSize and rasterYSize

        '''
        # recalculate GeoTransform based on extent option
        minX = extentDic['te'][0]
        minY = extentDic['te'][1]
        maxX = extentDic['te'][2]
        maxY = extentDic['te'][3]
        cornerX = minX
        cornerY = maxY
        width = maxX - minX
        height = maxY - minY
        if width <= 0 or height <= 0:
            raise OptionError('The extent is illegal. '
                              '"-te xMin yMin xMax yMax" ')

        if 'tr' in extentDic.keys():
            resolutionX = extentDic['tr'][0]
            resolutionY = -(extentDic['tr'][1])
            if (width < resolutionX or height < resolutionY):
                raise OptionError('"-tr" is too large. '
                                  'width is %s, height is %s '
                                  % (str(width), str(height)))
            rasterXSize = width / resolutionX
            rasterYSize = abs(height / resolutionY)
        else:
            rasterXSize = extentDic['ts'][0]
            rasterYSize = extentDic['ts'][1]
            resolutionX = width / rasterXSize
            resolutionY = -abs(height / rasterYSize)

        # create a list for GeoTransform
        coordinates = [cornerX, resolutionX, 0.0, cornerY, 0.0, resolutionY]

        return coordinates, int(rasterXSize), int(rasterYSize)

    def transform_points(self, colVector, rowVector, DstToSrc=0,
                         dstSRS=NSR()):

        '''Transform given lists of X,Y coordinates into lon/lat or inverse

        Parameters
        -----------
        colVector : lists
            X and Y coordinates in pixel/line or lon/lat  coordinate system
        DstToSrc : 0 or 1
            0 - forward transform (pix/line => lon/lat)
            1 - inverse transformation
        dstSRS : NSR
            destination spatial reference
            
        Returns
        --------
        X, Y : lists
            X and Y coordinates in lon/lat or pixel/line coordinate system

        '''
        return self.vrt.transform_points(colVector, rowVector,
                                         DstToSrc, dstSRS=dstSRS)

    def azimuth_y(self, reductionFactor=1):
        '''Calculate the angle of each pixel position vector with respect to
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

        '''

        lon, lat = self.get_geolocation_grids(reductionFactor)
        a = initial_bearing(lon[1:, :], lat[1:, :],
                            lon[:-1:, :], lat[:-1:, :])
        # Repeat last row once to match size of lon-lat grids
        a = np.vstack((a, a[-1, :]))
        return a

    def shape(self):
        '''Return Numpy-like shape of Domain object (ySize, xSize)

        Returns
        --------
        shape : tuple of two INT
            Numpy-like shape of Domain object (ySize, xSize)

        '''
        return self.vrt.dataset.RasterYSize, self.vrt.dataset.RasterXSize

    def write_map(self, outputFileName,
                  lonVec=None, latVec=None, lonBorder=10., latBorder=10.,
                  figureSize=(6, 6), dpi=50, projection='cyl', resolution='c',
                  continetsColor='coral', meridians=10, parallels=10,
                  pColor='r', pLine='k', pAlpha=0.5, padding=0.,
                  merLabels=[False, False, False, False],
                  parLabels=[False, False, False, False],
                  pltshow=False,
                  labels=None):
        ''' Create an image with a map of the domain

        Uses Basemap to create a World Map
        Adds a semitransparent patch with outline of the Domain
        Writes to an image file

        Parameters
        -----------
        outputFileName : string
            name of the output file name
        lonVec : [floats] or [[floats]]
            longitudes of patches to display
        latVec : [floats] or [[floats]]
            latitudes of patches to display
        lonBorder : float
            10, horisontal border around patch (degrees of longitude)
        latBorder : float
            10, vertical border around patch (degrees of latitude)
        figureSize : tuple of two integers
            (6, 6), size of the generated figure in inches
        dpi: int
            50, resolution of the output figure (size 6,6 and dpi 50
            produces 300 x 300 figure)
        projection : string, one of Basemap projections
            'cyl', projection of the map
        resolution : string, resolution of the map
            'c', crude
            'l', low
            'i', intermediate
            'h', high
            'f', full
        continetsColor : string or any matplotlib color representation
            'coral', color of continets
        meridians : int
            10, number of meridians to draw
        parallels : int
            10, number of parallels to draw
        pColor : string or any matplotlib color representation
            'r', color of the Domain patch
        pLine : string or any matplotlib color representation
            'k', color of the Domain outline
        pAlpha : float 0 - 1
            0.5, transparency of Domain patch
        padding : float
            0., width of white padding around the map
        merLabels : list of 4 booleans
            where to put meridian labels, see also Basemap.drawmeridians()
        parLables : list of 4 booleans
            where to put parallel labels, see also Basemap.drawparallels()
        labels : list of str
            labels to print on top of patches
        '''
        # if lat/lon vectors are not given as input
        if lonVec is None or latVec is None or len(lonVec) != len(latVec):
            lonVec, latVec = self.get_border()

        # convert vectors to numpy arrays
        lonVec = np.array(lonVec)
        latVec = np.array(latVec)

        # estimate mean/min/max values of lat/lon of the shown area
        # (real lat min max +/- latBorder) and (real lon min max +/- lonBorder)
        minLon = max(-180, lonVec.min() - lonBorder)
        maxLon = min(180, lonVec.max() + lonBorder)
        minLat = max(-90, latVec.min() - latBorder)
        maxLat = min(90, latVec.max() + latBorder)
        meanLon = lonVec.mean()
        meanLat = latVec.mean()

        # generate template map (can be also tmerc)
        plt.figure(num=1, figsize=figureSize, dpi=dpi)
        bmap = Basemap(projection=projection,
                       lat_0=meanLat, lon_0=meanLon,
                       llcrnrlon=minLon, llcrnrlat=minLat,
                       urcrnrlon=maxLon, urcrnrlat=maxLat,
                       resolution=resolution)

        # add content: coastline, continents, meridians, parallels
        bmap.drawcoastlines()
        bmap.fillcontinents(color=continetsColor)
        bmap.drawmeridians(np.linspace(minLon, maxLon, meridians),
                           labels=merLabels, fmt='%2.1f')
        bmap.drawparallels(np.linspace(minLat, maxLat, parallels),
                           labels=parLabels, fmt='%2.1f')

        # convert input lat/lon vectors to arrays of vectors with one row
        # if only one vector was given
        if len(lonVec.shape) == 1:
            lonVec = [lonVec]
            latVec = [latVec]

        for i in range(len(lonVec)):
            # convert lat/lons to map units
            mapX, mapY = bmap(list(lonVec[i].flat), list(latVec[i].flat))

            # from x/y vectors create a Patch to be added to map
            boundary = Polygon(zip(mapX, mapY),
                               alpha=pAlpha, ec=pLine, fc=pColor)

            # add patch to the map
            plt.gca().add_patch(boundary)
            plt.gca().set_aspect('auto')

            if labels is not None and labels[i] is not None:
                plt.text(np.mean(mapX), np.mean(mapY), labels[i],
                         va='center', ha='right', alpha=0.5, fontsize=10)

        # save figure and close
        plt.savefig(outputFileName, bbox_inches='tight',
                    dpi=dpi, pad_inches=padding)
        if pltshow:
            plt.show()
        else:
            plt.close('all')

    def reproject_GCPs(self, srsString=''):
        '''Reproject all GCPs to a new spatial reference system

        Necessary before warping an image if the given GCPs
        are in a coordinate system which has a singularity
        in (or near) the destination area (e.g. poles for lonlat GCPs)

        Parameters
        ----------
        srsString : string
            SRS given as Proj4 string. If empty '+proj=stere' is used

        Modifies
        --------
            Reprojects all GCPs to new SRS and updates GCPProjection
        '''
        if srsString == '':
            lon, lat = self.get_border()
            srsString = '+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=%f +lon_0=%f +no_defs'%(
            np.nanmedian(lat), np.nanmedian(lon)) 
        
        
        self.vrt.reproject_GCPs(srsString)
