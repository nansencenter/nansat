# Name:    domain.py
# Purpose: create domain based on either GDAL dataset or
#          proj4string and extentString
#
# Authors:      Asuka Yamakava, Anton Korosov, Knut-Frode Dagestad
#
# Created:     15.09.2011
# Copyright:   (c) NERSC 2012
# Licence:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details:
# http://www.gnu.org/licenses/

import os.path
import re
import string
import sys

from numpy  import *
from random import random
from xml.etree.ElementTree import *

try:
    from osgeo import gdal, osr
except ImportError:
    import gdal
    import osr

class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass;

class OptionError(Error):
    '''Error for improper options (arguments) '''
    pass;

class ProjectionError(Error):
    '''Cannot get the projection'''
    pass

class Domain():

    def __init__(self, dataset=None, srsString=None, extentString=None,
                 domainName=''):
        '''Construct Domain object

        Arguments should be one of the following:
            * inputDataset
            * srsString, extentString
                (in this case, "-te" or "-lle" and "-tr" or "-ts" are required)
        rasterXsize, rasterYsize, GeoTransform and GeoPorjection
        are fetched / computed based on the arguments.
        Then create empty memDataset and set the these values.

        Parameters
        ----------
        dataset : GDAL dataset, optional
        srsString : string, optional
            proj4string (e.g. "+proj=utm +zone=25 +datum=WGS84 +no_defs")
        extentString: string, optional
            Represent extent, resplution / size
            Available options are "-te","-tr","-ts" and "-lle"
            (e.g. "-lle -10 30 55 60 -ts 1000 1000")
        domainName: string, optional
            Domain name

        Raises
        ------
        ProjectionError: occurs when Projection() is empty
            despite it is required for creating extentDic.
        OptionError: occures when the arguments are not proper.

        Modifies
        --------
        self.memDatasetset: dataset in memory
            dataset created based on the arguments

        See Also
        --------
        http://www.gdal.org/gdalwarp.html

        '''

        # defaults
        gcps = []
        rasterXSize = 0
        rasterYSize = 0
        self.name = domainName

        # test option when only dataset is given
        if (dataset is not None and extentString is None):
            rasterXSize = dataset.RasterXSize
            rasterYSize = dataset.RasterYSize
            transform = dataset.GetGeoTransform()
            projection = dataset.GetProjection()
            if projection == "":
                projection = dataset.GetGCPProjection()
            gcps = dataset.GetGCPs()

        # test option when proj4 and extent string are given
        elif (srsString is not None and extentString is not None):
            # if XML-file and domain name is given - get string from that file
            if os.path.isfile(srsString):
                srsString, extentString, self.name = self._from_xml(
                                                          srsString,
                                                          extentString)
            # import srs from srsString and get the projection
            srs = osr.SpatialReference()
            srs.ImportFromProj4(srsString)
            projection = srs.ExportToWkt()
            # create full dictionary of parameters
            extentDic = self._create_extentDic(extentString)
            if "lle" in extentDic.keys():
                WKT = srs.ExportToWkt()
                if WKT == "":
                    raise ProjectionError("Domain.__init__() : "
                                           "WKT is empty. "
                                           "Check 'srsString'!! ")
                extentDic = self._convertExtentDic(WKT, extentDic)
            transform, rasterSize= self._get_geotransform(extentDic)
            rasterXSize = rasterSize[0]
            rasterYSize = rasterSize[1]

        else:
            raise OptionError("'dataset' or 'srsString and extentString' "
                              "are required")

        if rasterXSize != 0 or rasterYSize != 0:
            # If everything is fine so far
            # Create an empty memDataset and set rasterXSize, rasterYSize,
            # GsoTransform and Projection
            memDriver = gdal.GetDriverByName("MEM")
            self.memDataset = memDriver.Create("warped.mem",
                                               rasterXSize, rasterYSize, 0)
            self.memDataset.SetGeoTransform(transform)
            self.memDataset.Coordinates = None
            if len(gcps) > 0:
                self.memDataset.SetGCPs(gcps, projection)
            else:
                self.memDataset.SetProjection(projection)
        else:
            #if rasterXSize or rasterYSize are bad create empty Domain Object
            self.memDataset = None

    def __repr__(self):
        '''Prints basic info about the Domain object to the terminal

        Modifies
        --------
        Print size, projection and corner coordinates

        '''
        toPrettyWKT = osr.SpatialReference()
        toPrettyWKT.ImportFromWkt(self._get_projection(self.memDataset))
        prettyWKT  = toPrettyWKT.ExportToPrettyWkt(1)
        corners = self._get_corners()
        print '-'*40
        print 'Size: %d x %d' % (self.memDataset.RasterXSize,
                                 self.memDataset.RasterYSize)
        print '-'*40
        print 'Projection:'
        print prettyWKT;
        print '-'*40
        print 'Corners:'
        print 'Lat: %5.2f %5.2f %5.2f %5.3f' % (corners[0][0],
                corners[0][1], corners[0][2], corners[0][3])
        print 'Lon: %5.2f %5.2f %5.2f %5.3f' % (corners[1][0],
                corners[1][1], corners[1][2], corners[1][3])

        return ''

    def write_kml(self, xmlFileName=None, kmlFileName=None):
        '''Write KML file with domains

        Convert XML-file with domains into into KML-file for GoogleEart
        or
        Write KML-file with the current Domain

        Parameters
        ----------
        xmlFileName: string
            Name of the XML-file to convert. If only this value is given
            - kmlFileName=xmlFileName+'.kml'

        kmlFileName: string
            Name of the KML-file to generate from the current Domain

        Modifies
        --------
        Generates KML file

        '''
        # test input options
        if xmlFileName is not None and kmlFileName is None:
            # if only input XML-file is given - convert it to KML

            # open XML, get all domains
            xmlFile = file(xmlFileName, "rb")
            kmlFileName = xmlFileName + '.kml'

            xmlDomains = ElementTree(file=xmlFile).getroot()
            # convert domains in XML into list of domains
            domains = []
            for xmlDomain in list(xmlDomains):
                # append Domain object to domains list
                domainName = xmlDomain.attrib['name']
                domains.append(Domain(None, xmlFileName, domainName))

        elif xmlFileName is None and kmlFileName is not None:
            # if only output KML-file is given convert the current domain to KML
            domains = [self]

        else:
            # otherwise it is potentially error
            raise

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

    def _convertExtentDic(self, dstWKT, extentDic):
        ''' Convert values in degree of lat/lon to proper coordinate system

        Get source SRS from element XML and target SRS from dstWKT.
        Create osr.CoordinateTransformation based on these SRSs and
        convert given values in degree to the destination coordinate system
        given by WKT.
        Add key "te" and the converted values into the extentDic.

        Parameters
        ----------
        dstWKT: WKT
            destination WKT
        extentDic: dictionary
            dictionary with "lle" key

        Returns
        -------
        extentDic: dictionary
            dictionary with "te" key and its values

        Modifies
        --------
        Degree of lat/lon is converted to the given coordinate system.
        The converted values are added to the extentDic with "-te" key.

        '''
        # Set source SRS from XML element
        srcSRS = osr.SpatialReference()
        srcString = "+proj=latlong +datum=WGS84 +no_defs"
        srcSRS.ImportFromProj4(srcString)

        # Set destination SRS from dstWKT
        dstSRS = osr.SpatialReference()
        dstSRS.ImportFromWkt(dstWKT)

        CoorTrans = osr.CoordinateTransformation(srcSRS, dstSRS)

        # convert lat/lon given by "lle" to the target coordinate system and
        # add key "te" and the converted values to extentDic
        x1, y1, z1 = CoorTrans.TransformPoint(extentDic["lle"][0],
                                              extentDic["lle"][2])
        x2, y2, z2 = CoorTrans.TransformPoint(extentDic["lle"][0],
                                              extentDic["lle"][3])
        x3, y3, z3 = CoorTrans.TransformPoint(extentDic["lle"][1],
                                              extentDic["lle"][2])
        x4, y4, z4 = CoorTrans.TransformPoint(extentDic["lle"][1],
                                              extentDic["lle"][3])
        minX = min([x1, x2, x3, x4])
        maxX = max([x1, x2, x3, x4])
        minY = min([y1, y2, y3, y4])
        maxY = max([y1, y2, y3, y4])

        val = [minX, minY, maxX, maxY]
        extentDic["te"] = val

        return extentDic

    def _create_extentDic(self, extentString):
        '''Create a dictionary from extentString

        Check if extentString is proper.
            * "-te" and "-lle" take 4 numbers.
            * "-ts" and "-tr" take 2 numbers.
            * the combination should be ("-te" or "-lle") and ("-ts" or "-tr")
        If it is proper, create a dictionary
        Otherwise, raise the error.

        Parameters
        ----------
        extentString: string
            "-te xMin yMin xMax yMax",
            "-tr xResolution yResolution",
            "-ts width height",
            "-lle lonWest lonEast latNorth latSouth"

        Returns
        -------
        optionDic: dictionary
            has key ("te" or "lle") and ("tr" or "ts") and their values.

        Raises
        ------
        OptionError: occurs when the extentString is improper

        '''
        optionDic = {}

        # Find -re text
        str_tr = re.findall('-tr\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?',
                            extentString)
        if str_tr != []:
            # Check the number of -tr elements
            elm_str = str(str_tr[0].rstrip())
            elms_str = elm_str.split(None)
            if len(elms_str) != 3 or elms_str[2] == "-":
                raise OptionError("Domain._create_extentDic(): "
                                  "-tr is used as "
                                  "'-tr xResolution yResolution'")
            # Add the key and value to optionDic
            extentString = extentString.replace(str_tr[0],"")
            trElem = str(str_tr).split(None)
            trkey = trElem[0].translate(string.maketrans("", ""), "[]-'")
            if trkey != "":
                elements = []
                for i in range(2):
                    elements.append(float(trElem[i+1].\
                                          translate(string.maketrans("", ""),
                                          "[]'")))
                optionDic[trkey] = elements

        # Find -ts text
        str_ts = re.findall('-ts\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?',
                            extentString)
        if str_ts != []:
            # Check the number of -ts elements
            elm_str = str(str_ts[0].rstrip())
            elms_str = elm_str.split(None)
            if len(elms_str) != 3 or elms_str[2] == "-":
                raise OptionError("Domain._create_extentDic(): "
                                  "-ts is used as '-ts width height'")
            # Add the key and value to optionDic
            extentString = extentString.replace(str_ts[0],"")
            tsElem = str(str_ts).split(None)
            tskey = tsElem[0].translate(string.maketrans("", ""), "[]-'")
            if tskey != "":
                elements = []
                for i in range(2):
                    elements.append(float(tsElem[i+1].\
                                          translate(string.maketrans("", ""),
                                          "[]'")))
                optionDic[tskey] = elements

        # Find -te text
        str_te = re.findall('-te\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s'
                            '+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?',
                            extentString)
        if str_te != []:
            # Check the number of -te elements
            elm_str = str(str_te[0].rstrip())
            elms_str = elm_str.split(None)
            if len(elms_str) != 5:
                raise OptionError("Domain._create_extentDic(): "
                                  "-te is used as '-te xMin yMin xMax yMax'")
            # Add the key and value to optionDic
            extentString = extentString.replace(str_te[0],"")
            teElem = str(str_te).split(None)
            tekey = teElem[0].translate(string.maketrans("", ""), "[]-'")
            if tekey!= "":
                elements = []
                for i in range(4):
                    elements.append(float(teElem[i+1].\
                                          translate(string.maketrans("", ""),
                                          "[]'")))
                optionDic[tekey] = elements

        # Find -lle text
        str_lle = re.findall('-lle\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s'
                             '+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?',
                             extentString)
        if str_lle != []:
            # Check the number of -lle elements
            elm_str = str(str_lle[0].rstrip())
            elms_str = elm_str.split(None)
            if len(elms_str) != 5:
                raise OptionError("Domain._create_extentDic(): "
                                  "-lle is used as "
                                  "'-lle lonWest lonEast latNorth latSouth'")
            # Add the key and value to optionDic
            extentString = extentString.replace(str_lle[0], "")
            lleElem = str(str_lle).split(None)
            llekey = lleElem[0].translate(string.maketrans("", ""), "[]-'")
            if llekey != "":
                elements = []
                for i in range(4):
                    elements.append(float(lleElem[i+1].\
                                          translate(string.maketrans("", ""),
                                          "[]'")))
                optionDic[llekey] = elements

        result = re.search("\S", extentString)

        # if there are unnecessary letters, give an error
        if result != None:
            raise OptionError("Domain._create_extentDic(): "
                              "extentString is not redable : ", extentString)

        # check if one of "-te" and "-lle" is given
        if ("lle" not in optionDic) and ("te" not in optionDic):
            raise OptionError("Domain._create_extentDic(): "
                              "'-lle' or '-te' is required.")
        elif ("lle" in optionDic) and ("te" in optionDic):
            raise OptionError("Domain._create_extentDic(): "
                              "'-lle' or '-te' should be chosen.")

        # check if one of "-ts" and "-tr" is given
        if ("ts" not in optionDic) and ("tr" not in optionDic):
            raise OptionError("Domain._create_extentDic(): "
                              "'-ts' or '-tr' is required.")
        elif ("ts" in optionDic) and ("tr" in optionDic):
            raise OptionError("Domain._create_extentDic(): "
                              "'-ts' or '-tr' should be chosen.")
        return optionDic

    def _from_xml(self, srsString, extentString):
        ''' Read strings from the given xml file

        Parameters
        ----------
        srsString: file name
            name of the input XML-file
        extentString: string
            name of the domain

        Returns
        -------
        srsString: string
            proj4 string of the destination
        extentString: string
            extent string of the destination
        name: string
            domain name

        Raises
        ------
        OptionError: occures when the given extentString is not
        in the XML-file

         '''
        # open file
        fd = file(srsString, "rb")
        # get root element
        domains = ElementTree(file=fd).getroot()
        # iterate over domains to find the required one
        # for domain in domains.getchildren():
        for domain in list(domains):
            # if the domain name is the same as the given one
            if domain.attrib['name'] == extentString:
                # get contents of the tags
                name = extentString[:]
                srsString = domain.find('srsString').text
                extentString = domain.find('extentString').text
                break
            if domain == list(domains)[-1]:
                raise OptionError("extentString is improper")

        return srsString, extentString, name

    def _get_border(self, nPoints=10):
        '''Generate two vectors with values of lat/lon for the border of domain

        Parameters
        ----------
        nPoints: int
            Number of points on each border

        Returns
        -------
        lonVec, latVec: lists
            vectors with lon/lat values for each point at the border

        '''
        # prepare vectors with pixels and lines for upper, left, lower
        # and right borders
        sizes = [self.memDataset.RasterXSize, self.memDataset.RasterYSize]

        rcVector1 = [[],[]]
        rcVector2 = [[],[]]
        # loop for pixels and lines
        for n in range(0, 2):
            step = sizes[n] / nPoints
            rcVector1[n] = range(0, sizes[n], step)[0:nPoints]
            rcVector1[n].append(sizes[n])
            rcVector2[n] = rcVector1[n][:]
            rcVector2[n].reverse()

        # coumpund vectors of pixels (col) and lines (row)
        colVector = rcVector1[0] + [sizes[0]] * nPoints + \
                    rcVector2[0] + [0] * nPoints
        rowVector = [0] * nPoints + rcVector1[1] + \
                    [sizes[1]] * nPoints + rcVector2[1]

        return self._transfrom_points(colVector, rowVector)

    def _get_border_kml(self):
        '''Generate Placemark entry for KML

        Returns
        -------
        kmlEntry: String
            String with the Placemark entry

        '''
        domainLon, domainLat = self._get_border()

        # convert Border coordinates into KML-like string
        coordinates = ''
        for lon,lat in zip(domainLon, domainLat):
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
        kmlEntry += coordinates+'\n'
        kmlEntry += '            </coordinates></LinearRing>'\
                    '</outerBoundaryIs></Polygon></Placemark>\n'

        return kmlEntry

    def _get_border_polygon(self):
        '''Creates string with WKT representation of the border polygon
        (this method is not used. Delete??)

        Returns
        -------
        WKTPolygon: string
            string with WKT representation of the border polygon

        '''
        lonList, latList = self._get_border()
        polyCont = ','.join(str(lon) + ' ' +str(lat) \
                   for lon,lat in zip(lonList, latList))
        WKTPolygon = 'PolygonFromText("POLYGON((%s))")' % polyCont
        return WKTPolygon

    def _get_corners(self):
        '''Get coordinates of corners of the domain

        Returns
        -------
        lonVec, latVec: lists
            vectors with lon/lat values for each corner

        '''
        sizes = [self.memDataset.RasterXSize, self.memDataset.RasterYSize]
        colVector = [0, 0, self.memDataset.RasterXSize,
                     self.memDataset.RasterXSize]
        rowVector = [0, self.memDataset.RasterYSize, 0,
                     self.memDataset.RasterXSize]
        return self._transfrom_points(colVector, rowVector)

    def _get_geotransform(self, extentDic):
        '''
        the new coordinates and raster size are calculated based on
        the given extentDic.

        Parameters
        ----------
        extentDic : dictionary
            includes "te" key and "ts" or "tr" key

        Raises
        ------
        OptionError: occurs when maxX - minX < 0 or maxY - minY < 0
        OptionError: occurs when the given resolution is larger than
                     width or height.

        Returns
        -------
        coordinate: list with 6 float
            GeoTransform

        rasterSize : list with two int
            rasterXSize and rasterYSize

        '''
        # recalculate GeoTransform based on extent option
        minX = extentDic["te"][0]
        minY = extentDic["te"][1]
        maxX = extentDic["te"][2]
        maxY = extentDic["te"][3]
        cornerX = minX
        cornerY = maxY
        width = maxX - minX
        height = maxY - minY
        if width <= 0 or height <= 0:
            raise OptionError("The extent is illigal. "
                              "'-te xMin yMin xMax yMax' ")

        if "tr" in extentDic.keys():
            resolutionX = extentDic["tr"][0]
            resolutionY = -(extentDic["tr"][1])
            if (width < resolutionX or height < resolutionY):
                raise OptionError("'-tr' is too large. "
                                  "width is " + width + "  height is "+ height)
            rasterXSize = width / resolutionX
            rasterYSize = abs(height / resolutionY)
        else:
            rasterXSize = extentDic["ts"][0]
            rasterYSize = extentDic["ts"][1]
            resolutionX = width / rasterXSize
            resolutionY = abs(height / rasterYSize)

        # create a list for GeoTransform
        coordinates = [cornerX, resolutionX, 0.0, cornerY, 0.0, resolutionY]
        rasterSize = [int(rasterXSize), int(rasterYSize)]

        return coordinates, rasterSize

    def _get_projection(self, dataset):
        '''get projection form dataset

        Get projection from GetProjection() or GetGCPProjection().
        If both are empty, raise error

        Return
        ------
        projection : projection or GCPprojection

        Raises
        ------
        ProjectionError: occurrs when the projection is empty.

        '''
        projection = dataset.GetProjection()
        if projection != "":
            return projection

        projection = dataset.GetGCPProjection()
        if projection != "":
            return projection
        else:
            raise ProjectionError()

    def _transfrom_points(self, colVector, rowVector):
        '''Transform given lists of X,Y coordinates into lat/lon

        Parameters
        ----------
        colVector: lists
            X and Y coordinates with any coordinate system

        Returns
        -------
        lonVector, latVector: lists
            X and Y coordinates in degree of lat/lon

        '''
        # get source SRS (either Projection or GCPProjection)
        srcWKT = self.memDataset.GetProjection()
        if srcWKT == '':
            srcWKT = self.memDataset.GetGCPProjection()

        # prepare target WKT (pure lat/lon)
        proj4string = "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs"
        dstSRS = osr.SpatialReference()
        dstSRS.ImportFromProj4(proj4string)
        dstWKT = dstSRS.ExportToWkt()

        # create transformer
        transformer = gdal.Transformer(self.memDataset, None,
                                       ['SRC_SRS='+srcWKT, 'DST_SRS='+dstWKT])

        # use the transformer to convert pixel/line into lat/lon
        latVector = []
        lonVector = []
        for pixel, line in zip(colVector, rowVector):
            succ,point = transformer.TransformPoint(0, pixel, line)
            lonVector.append(point[0])
            latVector.append(point[1])

        return lonVector, latVector

