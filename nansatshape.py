# Name:    nansashape.py
# Purpose: Container of Nansatshape class
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      07.05.2013
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

# import standard and additional libraries
from nansat_tools import *

class Nansatshape():
    ''' Nansatshape class reads and writes ESRI-shape files

        The core of Nansatshape is a OGR. the main functions of the class are
        1. Create empty object in memory and add data (fields and geometory).
        2. Open shape file and read the data.

        Nansatshape support points, line and ploygons. (not mupti-polygon)

    '''
    def __init__(self, fileName=None, layer=0, srs=None, wkbStyle=ogr.wkbPoint):
        '''Create Nansatshape object

        if <fileName> is given:
            Open OGR datasource and copy to self.datasource in memory
            Read a layer from self.datasource and add it to self.layer
        if <fileName> is not given:
            Create empty OGR datasource in memory
            Add empty layer to self.layer

        Parameters
        -----------
        fileName : string
            location of a shape file
        layer : int or string
            if int and a shapefile is given, it is a layer number which is read.
            if string, it is layer name to created or open
        srs : SpatialReference object
        wkbStyle : ogr.wkbPoint, ogr.wkbPoint25D, ogr.wkbLineString,
                   ogr.wkbLineString25D, ogr.wkbPolygon or ogr.wkbPolygon25D

        Creates
        --------
        self.fieldNames : list of field names
            list of field names in the layer
        self.fileDataType : dictionary
            keys are field names and values are ogr data types.
        self.datasource : ogr data source in memory
        self.layer : ogr layer

        '''
        # Create a empty datasource and layer in memory
        if fileName is None:
            self.datasource= ogr.GetDriverByName('Memory').CreateDataSource('wrk')
            # create a new later
            if layer == 0:
                layer = 'NansatLayer'
            self.layer = self.datasource.CreateLayer(layer, srs, wkbStyle)
        # if filename is given, open the file and copy the datasource in memory
        else:
            # Open shapefile and copy the datasource into memory
            ogrDs = ogr.Open(fileName)
            self.datasource= ogr.GetDriverByName('Memory').CopyDataSource(ogrDs, 'OGRDataSource')
            ogrDs.Destroy()
            # Set a layer from the datasource
            if type(layer) is int:
                self.layer = self.datasource.GetLayer(layer)
            else:
                self.layer = self.datasource.GetLayerByName(layer)
            # Set self.fieldNames and self.fieldDataTypes
            layerDef = self.layer.GetLayerDefn()
            feature = self.layer.GetFeature(0)

    def add_features(self, values, coordinates):
        ''' Set field values and / or geometry to each feature

        !! NB !! Muptipolygon is not supported

        Parameters
        ----------
        values :  structured np.array
            field names and values
            e.g. : np.zeros(n, dtype={'names':['INT', 'STR','FLOAT1','FLOAT2'],
                                      'formats':['i4','a10','f8','f8']})
                    n is a number of features.

        coordinates : structured np.array (2 or 3 x n)
            If geometries are points,
            the num of column is 2 and X-Y coordinates are given.
            If geometries are lines /  polygon,
            the num of column is 3 and X-Y coordinates and ID of each points
            are given.
            Names of the array are 'pixel', 'line', 'ID' respectively.

        Modifies
        --------
        self.layer : Set the values and geometry.

        '''
        # get geometry type
        geomType = self.layer.GetGeomType()

        # add pixel and line coordinates info to values
        if geomType == ogr.wkbPoint or geomType == ogr.wkbPoint25D and \
            len(coordinates) == len(values) and(
            'pixel' in coordinates.dtype.names and
            'line' in coordinates.dtype.names):

            # append pix/lin fields to values if values is given
            descr = [('X (pixel)', int), ('Y (line)', int)]
            # Create Zero structrued np.array ( n x m )
            # n is num of feature of values
            # m is num of field of values + 2 (for 'X (pixel)' and 'Y (line)' )
            tmp = np.zeros(values.shape[0],
                           dtype=values.dtype.descr + descr)
            # copy values to tmp
            for name in values.dtype.names:
                tmp[name] = values[name]
            values = tmp

            # set coordinates into values
            for iCoord in range(coordinates.shape[0]):
                values['X (pixel)'][iCoord] = coordinates[iCoord]['pixel']
                values['Y (line)'][iCoord] = coordinates[iCoord]['line']

        # set fieldDenf
        for i, iFieldName in enumerate (values.dtype.names):
            # get data type for each field
            if str(values.dtype[i]).startswith('int'):
                dtype = ogr.OFTInteger
            elif str(values.dtype[i]).startswith('float'):
                dtype = ogr.OFTReal
            else:
                dtype = ogr.OFTString
            # create fields
            field_defn = ogr.FieldDefn(iFieldName, ogr.OFTString)
            field_defn.SetWidth(32)
            self.layer.CreateField(field_defn)

        # set values to each feature
        for iFeature in range(len(values)):
            feature = ogr.Feature(self.layer.GetLayerDefn())
            # create a geometry
            geometry = ogr.Geometry(type=geomType)
            # if geomType is Point, set a point to a geometry
            if geomType == ogr.wkbPoint or geomType == ogr.wkbPoint25D:
                geometry.SetPoint_2D(0,
                                     float(coordinates[iFeature]['pixel']),
                                     float(coordinates[iFeature]['line']))
            # else, set points to a geometry
            else:
                geomRing = ogr.Geometry(type=ogr.wkbLinearRing)
                # add points for a line
                if geomType == ogr.wkbLineString or geomType == ogr.wkbLineString25D:
                    for iGeo in range(len(coordinates)):
                        if coordinates[iGeo]['ID'] == iFeature:
                            geometry.AddPoint(float(coordinates[iGeo]['pixel']),
                                              float(coordinates[iGeo]['line']))
                # add points for a polygon
                elif geomType == ogr.wkbPolygon or geomType == ogr.wkbPolygon25D:
                    for iGeo in range(len(coordinates)):
                        if coordinates[iGeo]['ID'] == iFeature:
                            geomRing.AddPoint(float(coordinates[iGeo]['pixel']),
                                              float(coordinates[iGeo]['line']))
                else:
                    print 'Can not create geometry. Muptipolygon is not supported.'
                # set geomRing if it has points (polygon)
                if geomRing.GetPointCount() != 0:
                    geometry.AddGeometryDirectly(geomRing)
            # set srs
            srs =self.layer.GetSpatialRef()
            geometry.AssignSpatialReference(srs)

            # set geometry
            feature.SetGeometryDirectly(geometry)

            # set field values
            for j, jField in enumerate(values.dtype.names):
                if values.dtype[jField].name.startswith("int"):
                    feature.SetField(jField, int(values[iFeature][jField]))
                elif values.dtype[jField].name.startswith("float"):
                    feature.SetField(jField, float(values[iFeature][jField]))
                else:
                    feature.SetField(jField, str(values[iFeature][jField]))

            self.layer.CreateFeature(feature)

            feature.Destroy()

    def get_corner_points(self, latlon=True):
        '''Get corner points (geometries of featuers) in the layer

        !!NB!!
        if shapefile has SRS, assume that geometry is lon/lat
        if not, assume that the geometry in pix/lin is given
        Muptipolygon is not supported

        Parameters
        ----------
        latlon : bool
            if True, coordinates in lon/lat

        Returns
        --------
        points : tuple or None
            elements of tuple are X-Y coordinates
        latlon : bool

        '''
        # get srs from the layer
        srs = self.layer.GetSpatialRef()

        points = []
        for feature in self.layer:
            geom = feature.GetGeometryRef()
            # if srs is None, get srs from the geometry
            if geom is not None and srs is None:
                srs = geom.GetSpatialReference()
            # if srs is given, assume the coordinates are in lat/lon
            if srs is None:
                latlon = False
            else:
                latlon = True
            # get corner points from geometry
            if geom is not None and (geom.GetGeometryType() == ogr.wkbPoint or
                                     geom.GetGeometryType() == ogr.wkbPoint25D):
                p = geom.GetPoints()[0]
                points.append((p[0], p[1]))
            elif geom is not None and (geom.GetGeometryType() == ogr.wkbLineString or
                                       geom.GetGeometryType() == ogr.wkbLineString25D):
                for iPoint in range(geom.GetPointCount()):
                    p = geom.GetPoints()[iPoint]
                    points.append((p[0], p[1]))
            elif geom is not None and (geom.GetGeometryType() == ogr.wkbPolygon or
                                       geom.GetGeometryType() == ogr.wkbPolygon25D):
                ring = geom.GetGeometryRef(0)
                for iPoint in range(ring.GetPointCount()):
                    p = ring.GetPoint(iPoint)
                    points.append((p[0], p[1]))
            # if shapefile has no geometry or multipolygon, give warning
            else:
                print 'No supported geometry in shape file. (multi-ploygon is not supported)'
        if points == []:
            points = None
        else:
            points = tuple(points)
        return points, latlon
