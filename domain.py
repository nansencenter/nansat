#-------------------------------------------------------------------------------
# Name:        domain
# Purpose:
#
# Author:      asumak/antonk
#
# Created:     15.09.2011
# Copyright:   (c) NERSC 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
import sys
import string
import re

try:
    from osgeo import gdal
except ImportError:
    import gdal

try:
    from osgeo import osr
except ImportError:
    import osr

from xml.etree.ElementTree import *

from numpy  import *
from random import random

import os.path

class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass;

class OptionError(Error):
    '''Error for options (arguments)'''
    pass;

class ProjectionError(Error):
    ''' Error for projection '''
    pass;


class Domain():

    def __init__(self, ds = None, srsString = None, extentString = None, name=''):
        '''Construct Domain object

        Arguments should be one of the following:
            * inputDS
            * inputDS, extentString
            * srsString, extentString
                (in this case, "-te" or "-lle" and "-tr" or "-ts" are required)
        rasterXsize, rasterYsize, GeoTransform and GeoPorjection
        are fetched / computed based on the arguments.
        Then create empty memDs and set the these values.

        Args:
            ds : dataset
            srsString : proj4string
                        (e.g. "+proj=utm +zone=25 +datum=WGS84 +no_defs")
            extentString : extent string option
                           Available options are "-te","-tr","-ts" and "-lle"
                           (e.g. "-lle -10 30 55 60 -ts 1000 1000")
            name: the domain name

        Raises:
            ProjectionError: occurs when Projection() is empty despite it is required for creating extentDic.

        Side effects:
            set attributes : memDs

        '''

        #defaults
        gcps = [];
        Xsize = 0;
        Ysize = 0;
        self.name = name;

        #test option when only dataset is given
        if ds != None and extentString == None:
            Xsize = ds.RasterXSize;
            Ysize = ds.RasterYSize;
            transform = ds.GetGeoTransform();
            projection = ds.GetProjection();
            if projection == "":
                projection = ds.GetGCPProjection();
            gcps = ds.GetGCPs();

        #test option when dataset and extent string are given (should be excluded?)
        elif ds != None and extentString != None:
            extentDic = self.create_extentDic(extentString);
            try:
                projection = self.get_projection_ds(ds);
            except ProjectionError:
                    raise ProjectionError(\
                        "Both ds.GetProjection() and ds.GetGCPProjection() are empty")
            except:
                error_type, error_value, traceback = sys.exc_info();
                print "Domain.__init__.convertGeoSystem() : Unknown error"
                print (error_type);
                print (error_value);
                sys.exit();

            extentDic = self.convertGeoSystem(projection, extentDic);
            transform, size= self.get_geotransform(extentDic, dataset = ds);
            Xsize = size[0];
            Ysize = size[1];

        #test option when proj4 and extent string are given
        elif srsString != None and extentString != None:
            #if XML-file and domain name is given - get string from that file
            if os.path.isfile(srsString):
                srsString, extentString, self.name = self.from_xml(srsString, extentString);

            #import srs from srsString and get the projection
            srs = osr.SpatialReference();
            srs.ImportFromProj4(srsString);
            projection = srs.ExportToWkt();

            #create full dictionary of parameters
            extentDic = self.create_extentDic(extentString);

            if "lle" in extentDic.keys():
                try:
                    wkt = self.get_wkt_srs(srs);
                except ProjectionError:
                        raise ProjectionError(\
                            "Domain.__init__() : WKT is empty."\
                            " Check 'srsString'!! ")
                except:
                    error_type, error_value, traceback = sys.exc_info();
                    print "Domain.__init__.convertGeoSystem() : Unknown error"
                    print (error_type);
                    print (error_value);
                    sys.exit();
                extentDic = self.convertGeoSystem(wkt, extentDic);

            transform, size= self.get_geotransform(extentDic);
            Xsize = size[0];
            Ysize = size[1];

        if Xsize != 0 or Ysize != 0:
            # If everything is fine so far
            # Create an empty memDs and set Xsize, Ysize, GsoTransform and Projection
            memDriver = gdal.GetDriverByName("MEM");
            self.memDs = memDriver.Create("warped.mem", Xsize, Ysize, 0);
            self.memDs.SetGeoTransform(transform);
            self.memDs.Coordinates = None;
            if len(gcps) > 0:
                self.memDs.SetGCPs(gcps, projection);
            else:
                self.memDs.SetProjection(projection);
        else:
            #if Xsize or Ysize are bad create empty Domain Object
            self.memDs = None;

    def create_extentDic(self, extentString):
        '''Create a dictionary from extentString

        Check if extentString is proper.
            * "-te" and "-lle" take 4 numbers.
            * "-ts" and "-tr" take 2 numbers.
            * the combination should be ("-te" or "-lle") and ("-ts" or "-tr")
        If not, raise the error message
        If ok, create the dictionary

        Args:
            extentString : string
                            "-te xMin yMin xMax yMax",
                            "-tr xResolution yResolution",
                            "-ts width height",
                            "-lle lonWest lonEast latNorth latSouth"

        Return:
            dictionary

        Raises:
            OptionError: occurs when the extentString is not proper


        '''
        optionDic = {};

        # Find -re text
        str_tr = re.findall('-tr\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?', extentString);
        if str_tr != []:
            # Check the number of -tr elements
            elm_str = str(str_tr[0].rstrip());
            elms_str = elm_str.split(None);
            if len(elms_str) != 3 or elms_str[2] == "-":
                raise OptionError("Domain.create_extentDic(): -tr is used as '-tr xResolution yResolution'");
            # Add the key and value to optionDic
            extentString = extentString.replace(str_tr[0],"");
            trElem = str(str_tr).split(None);
            trkey = trElem[0].translate(string.maketrans("", ""), "[]-'");
            if trkey!= "":
                elements = [];
                for i in range(2):
                    elements.append(float(trElem[i+1].translate\
                                           (string.maketrans("", ""), "[]'")));
                optionDic[trkey] = elements;

        # Find -ts text
        str_ts = re.findall('-ts\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?', extentString);
        if str_ts != []:
            # Check the number of -ts elements
            elm_str = str(str_ts[0].rstrip());
            elms_str = elm_str.split(None);
            if len(elms_str) != 3 or elms_str[2] == "-":
                raise OptionError("Domain.create_extentDic(): -ts is used as '-ts width height'");
            # Add the key and value to optionDic
            extentString = extentString.replace(str_ts[0],"");
            tsElem = str(str_ts).split(None);
            tskey = tsElem[0].translate(string.maketrans("", ""), "[]-'");
            if tskey != "":
                elements = [];
                for i in range(2):
                    elements.append(float(tsElem[i+1].translate\
                                            (string.maketrans("", ""), "[]'")));
                optionDic[tskey] = elements;

        # Find -te text
        str_te = re.findall('-te\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?', extentString);
        if str_te != []:
            # Check the number of -te elements
            elm_str = str(str_te[0].rstrip());
            elms_str = elm_str.split(None);
            if len(elms_str) != 5:
                raise OptionError("Domain.create_extentDic(): -te is used as '-te xMin yMin xMax yMax'");
            # Add the key and value to optionDic
            extentString = extentString.replace(str_te[0],"");
            teElem = str(str_te).split(None);
            tekey = teElem[0].translate(string.maketrans("", ""), "[]-'");
            if tekey!= "":
                elements = [];
                for i in range(4):
                    elements.append(float(teElem[i+1].translate\
                                          (string.maketrans("", ""), "[]'")));
                optionDic[tekey] = elements;

        # Find -lle text
        str_lle = re.findall('-lle\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s+[-+]?\d*[.\d*]*\s?', extentString);
        if str_lle != []:
            # Check the number of -lle elements
            elm_str = str(str_lle[0].rstrip());
            elms_str = elm_str.split(None);
            if len(elms_str) != 5:
                raise OptionError("Domain.create_extentDic(): -lle is used as '-lle lonWest lonEast latNorth latSouth'");
            # Add the key and value to optionDic
            extentString = extentString.replace(str_lle[0],"");
            lleElem = str(str_lle).split(None);
            llekey = lleElem[0].translate(string.maketrans("", ""), "[]-'");
            if llekey != "":
                elements = [];
                for i in range(4):
                    elements.append(float(lleElem[i+1].translate\
                                            (string.maketrans("", ""), "[]'")));
                optionDic[llekey] = elements;

        result = re.search("\S", extentString);

        # if there are unnecessary letters, give an error
        if result != None:
            raise OptionError("Domain.create_extentDic(): extentString is not redable : ", extentString);

        # check if one of "-te" and "-lle" is given
        if ("lle" not in optionDic) and ("te" not in optionDic):
            raise OptionError("Domain.create_extentDic(): '-lle' or '-te' is required.");
        elif ("lle" in optionDic) and ("te" in optionDic):
            raise OptionError("Domain.create_extentDic(): '-lle' or '-te' should be chosen.");

        # check if one of "-ts" and "-tr" is given
        if ("ts" not in optionDic) and ("tr" not in optionDic):
            raise OptionError("Domain.create_extentDic(): '-ts' or '-tr' is required.");
        elif ("ts" in optionDic) and ("tr" in optionDic):
            raise OptionError("Domain.create_extentDic(): '-ts' or '-tr' should be chosen.");

        return optionDic;

    def convertGeoSystem(self, trgWkt, extentDic):
        ''' Convert values in degree of lat/lon to proper coordinate system

        Get source and target SRS from element XML.
        Create osr.CoordinateTransformation based on these SRSs and
        convert given values in degree to the target coordinate system.
        if "lle" is given,
        key "te" and the converted values are added to the extentDic.

        Args:
            extentDic : dictionary

        Returns:
            extentDic : dictionary

        '''
        # Set source SRS from XML element
        srcSRS = osr.SpatialReference();
        srcString = "+proj=latlong +datum=WGS84 +no_defs";
        srcSRS.ImportFromProj4(srcString);

        # Set target SRS from XML element
        trgSRS = osr.SpatialReference();
        trgSRS.ImportFromWkt(trgWkt);

        CoorTrans = osr.CoordinateTransformation(srcSRS, trgSRS);

        # convert lat/lon given by "lle" to the target coordinate system and
        # add key "te" and the converted values to extentDic
        if "lle" in extentDic.keys():
            x1, y1, z1 = CoorTrans.TransformPoint(extentDic["lle"][0], \
                                                  extentDic["lle"][2]);
            x2, y2, z2 = CoorTrans.TransformPoint(extentDic["lle"][0], \
                                                  extentDic["lle"][3]);
            x3, y3, z3 = CoorTrans.TransformPoint(extentDic["lle"][1], \
                                                  extentDic["lle"][2]);
            x4, y4, z4 = CoorTrans.TransformPoint(extentDic["lle"][1], \
                                                  extentDic["lle"][3]);
            minx = min([x1, x2, x3, x4]);
            maxx = max([x1, x2, x3, x4]);
            miny = min([y1, y2, y3, y4]);
            maxy = max([y1, y2, y3, y4]);

            val = [minx, miny, maxx, maxy];
            extentDic["te"] = val;
            print "(domain 312) lle --> te : ", val;

        return extentDic;

    def get_geotransform(self, extent, dataset = None):
        '''
        Get geotransform coordinate

        get geotransform coordinate from ds.GetGeoTransform().
        If extentDic is given,
        the new coordinates are calculated based on the given extentDic.

        Args:
            ds : dataset
            extentDic : dictionary

        Raises:
            OptionError: occurs when the Xsize or Ysise is zero becuause of Xsize < resolutionX or Ysize < resolutionY

        Return:
            coordinate: list

        '''
        # Get GeoTransform from the given dataset
        if dataset != None:
            geoTransCoord = dataset.GetGeoTransform();
            cornerX = geoTransCoord[0];
            cornerY = geoTransCoord[3];
            resolutionX = geoTransCoord[1];
            resolutionY = geoTransCoord[5];
            sizeX = dataset.RasterXSize;
            sizeY = dataset.RasterYSize;
            width = resolutionX * sizeX;
            height = abs(resolutionY * sizeY);

        # recalculate GeoTransform if extent option is given
        if "te" in extent.keys():
            xmin = extent["te"][0];
            ymin = extent["te"][1];
            xmax = extent["te"][2];
            ymax = extent["te"][3];
            cornerX = xmin;
            #cornerY = ymin;
            cornerY = ymax;
            width = xmax - xmin;
            height = ymax - ymin;
            if width <= 0 or height <= 0:
                raise OptionError("The extent is illigal. '-te xMin yMin xMax yMax' ")

        if "tr" in extent.keys():
            resolutionX = extent["tr"][0];
            resolutionY = -(extent["tr"][1]);
            sizeX = width / resolutionX;
            sizeY = abs(height / resolutionY);

        if "ts" in extent.keys():
            sizeX = extent["ts"][0];
            sizeY = extent["ts"][1];
            resolutionX = width / sizeX;
            resolutionY = -(height) / sizeY;

        # create a list for GeoTransform
        coordinates = [cornerX, resolutionX, 0.0, cornerY, 0.0, resolutionY];
        size = [int(sizeX), int(sizeY)];

        return coordinates, size;

    #def get_GCPs(self, mat, lon, lat, proj4string, grid_spacing, tie_point_reduction_factor=10):

    def from_xml(self, srsString, extentString):
        ''' Read strings from the given xml file

        Args:
            srsString: name of the input XML-file
            extentString: name of the domain

        Returns:
            srsString: the proj4 string
            extentString: the extent string
         '''
        #open file
        fd = file(srsString, "rb");
        #get root element
        domains = ElementTree(file=fd).getroot();
        #iterate over domains to find the required one
        #for domain in domains.getchildren():
        for domain in list(domains):
            #if the domain name is the same as the given one
            if domain.attrib['name'] == extentString:
                #get contents of the tags
                name = extentString[:];
                srsString = domain.find('srsString').text;
                extentString = domain.find('extentString').text;
                break;
            if domain == list(domains)[-1]:
                raise TypeError("extentString is unproper");

        return srsString, extentString, name;

    def get_border(self, nPoints=10):
        '''
        Generate two vectors with values of lat/lon for the border
        of domain.

        Args:
            nPoints: Number of points on each border

        Returns:
            lonVec, latVec: vectors with lon/lat values for each
                            point at the border
        '''

        #prepare vectors with pixels and lines for upper, left, lower
        #and right borders
        sizes = [self.memDs.RasterXSize, self.memDs.RasterYSize]

        rcVector1 = [[],[]];
        rcVector2 = [[],[]];
        #loop for pixels and lines
        for n in range(0, 2):
            step = sizes[n] / nPoints;
            rcVector1[n] = range(0, sizes[n], step)[0:nPoints];
            rcVector1[n].append(sizes[n]);
            rcVector2[n]  = rcVector1[n][:];
            rcVector2[n].reverse();

        #coumpund vectors of pixels (col) and lines (row)
        colVector = rcVector1[0] + [sizes[0]] * nPoints + \
                    rcVector2[0] + [0] * nPoints;
        rowVector = [0] * nPoints + rcVector1[1] + \
                    [sizes[1]] * nPoints + rcVector2[1];

        return self.transfrom_points(colVector, rowVector);

    def get_corners(self):
        '''Get coordinates of corners of the domain        Returns:
            lonVec, latVec: vectors with lon/lat values for each corner
        '''

        sizes = [self.memDs.RasterXSize, self.memDs.RasterYSize]
        colVector = [0, 0, self.memDs.RasterXSize, self.memDs.RasterXSize]
        rowVector = [0, self.memDs.RasterYSize, 0, self.memDs.RasterXSize]
        return self.transfrom_points(colVector, rowVector);

    def transfrom_points(self, colVector, rowVector):
        '''Transform given lists of X,Y coordinates into lat/lon'''

        #get source SRS (either Projection or GCPProjection)
        srcWkt = self.memDs.GetProjection()
        if srcWkt == '':
            srcWkt = self.memDs.GetGCPProjection()
        
        #prepare target WKT (pure lat/lon)
        proj4string = "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs"
        trgSRS = osr.SpatialReference();
        trgSRS.ImportFromProj4(proj4string);
        trgWkt = trgSRS.ExportToWkt();

        #create transformer
        transformer = gdal.Transformer(self.memDs, None, ['SRC_SRS='+srcWkt, 'DST_SRS='+trgWkt])

        #use the transformer to convert pixel/line into lat/lon
        latVec = [];
        lonVec = [];
        for pixel, line in zip(colVector, rowVector):
            succ,point = transformer.TransformPoint(0, pixel, line);
            lonVec.append(point[0]);
            latVec.append(point[1]);

        return lonVec, latVec

    def get_border_kml(self):
        '''Generate Placemark entry for KML

           Args:
                domainName: any string with name of the Domain

           Returns:
                kmlEntry: String with the Placemark entry
        '''
        domainLon, domainLat = self.get_border();

        #convert Border coordinates into KML-like string
        coordinates = ''
        for lon,lat in zip(domainLon, domainLat):
            coordinates += '%f,%f,0 ' % (lon, lat)

        kmlEntry = ''
        #write placemark: name, style, polygon, coordinates
        kmlEntry += '            <Placemark>\n';
        kmlEntry += '                <name>%s</name>\n' % self.name;
        kmlEntry += '                <Style>\n';
        #kmlEntry += '                    <LineStyle><color>ff%s</color></LineStyle>\n' % hex(int(random()*16777215))[2:];
        kmlEntry += '                    <LineStyle><color>ffffffff</color></LineStyle>\n';
        kmlEntry += '                    <PolyStyle><fill>0</fill></PolyStyle>\n';
        kmlEntry += '                </Style>\n';
        kmlEntry += '                <Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing><coordinates>\n';
        kmlEntry += coordinates+'\n';
        kmlEntry += '            </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n';

        return kmlEntry
    
    def get_border_polygon(self):
        '''Creates string with WKT representation of the border polygon
        
        Returns:
            wktPolygon: string with WKT representation of the border polygon'''
        
        lonList, latList = self.get_border();
        polyCont = ','.join(str(lon) + ' ' +str(lat) for lon,lat in zip(lonList, latList))
        wktPolygon = 'PolygonFromText("POLYGON((%s))")' % polyCont;
        return wktPolygon

    def write_kml(self, xmlFileName=None, kmlFileName=None):
        '''
        Convert XML-file with domains into into KML-file for GoogleEart
        or
        Write KML-file with the current Domain

        Args:
            xmlFileName: Name of the XML-file to convert. If only this
            value is given - kmlFileName=xmlFileName+'.kml'

            kmlFileName: Name of the KML-file to generate from the
            current Domain

        Returns: Nothing

        Side: Generates KML file
        '''
        #test input options
        if xmlFileName is not None and kmlFileName is None:
            #if only input XML-file is given - convert it to KML

            #open XML, get all domains
            xmlFile = file(xmlFileName, "rb");
            kmlFileName = xmlFileName + '.kml';

            xmlDomains = ElementTree(file=xmlFile).getroot();
            #convert domains in XML into list of domains
            domains = [];
            for xmlDomain in list(xmlDomains):
                #append Domain object to domains list
                domainName = xmlDomain.attrib['name']
                domains.append(Domain(None, xmlFileName, domainName));

        elif xmlFileName is None and kmlFileName is not None:
            #if only output KML-file is given convert the current domain to KML
            domains = [self]
        else:
            #otherwise it is potentially error
            raise

        #open KML, write header
        kmlFile = file(kmlFileName, 'wt');
        kmlFile.write('<?xml version="1.0" encoding="UTF-8"?>\n');
        kmlFile.write('<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n');
        kmlFile.write('<Document>\n');
        kmlFile.write('    <name>%s</name>\n' % kmlFileName);
        kmlFile.write('        <Folder><name>%s</name><open>1</open>\n' % kmlFileName);

        #get border of each domain and add to KML
        for domain in list(domains):
            kmlEntry = domain.get_border_kml();
            kmlFile.write(kmlEntry);

        #write footer and close
        kmlFile.write('        </Folder></Document></kml>\n')
        kmlFile.close();

    def get_projection_ds(self, ds):
        '''get projection form dataset

            Get projection from GetProjection() or GetGCPProjection().
            If both are empty, raise error

        Raises:
            ProjectionError: occurrs when the projection is empty.

        Return:
            projection
        '''
        projection = ds.GetProjection();
        if projection != "":
            return projection;

        projection = ds.GetGCPProjection();
        if projection != "":
            return projection;
        else:
            raise ProjectionError();

    def get_wkt_srs(self, srs):
        '''get WKT based on SRS

            Get WKT with ExportToWkt().
            If it is empty, raise error

        Raises:
            ProjectionError: occurrs when the WKT is empty.

        Return:
            WKT
        '''
        wkt = srs.ExportToWkt();
        if wkt !="":
            return wkt;
        else:
            raise ProjectionError();
