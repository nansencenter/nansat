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
        # Create attributes
        self.fieldNames = []
        self.fieldDataTypes = {}
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
            for i in range(layerDef.GetFieldCount()):
                self.fieldNames.append(layerDef.GetFieldDefn(i).GetName())
                self.fieldDataTypes[self.fieldNames[-1]] = feature.GetFieldType(i)

    def add_features(self, values=None, coordinates=None, fieldFeatureID=[], geoFeatureID=[], AddPixLine=True):
        ''' Set field values and / or geometry to each feature

        !! NB !! Muptipolygon is not supported

        Parameters
        ----------
        values :  structured np.array
            field names and values
            e.g. : np.zeros(n, dtype={'names':['INT', 'STR','FLOAT1','FLOAT2'],
                                      'formats':['i4','a10','f8','f8']})
                    n is a number of features.

        coordinates : np.array or list with scalar elements (2 x n)
            X-Y coordinates. n is a number of points

        fieldFeateurID : None or list with int elements (the length is n)
            the values are from 0 to n-1

        geoFeateurID : None or list with int elements (the length is n)
            if the data type is point, it will be [0, 1, 2, 3, ..., (num of fields-1)]
            if polygon or line, the numbers represent which geometry each point belongs to.

        AddPixLine : bool
            if True, add 'X (pixel)' and 'Y (line)' as the fields.

        Modifies
        --------
        self.layer : if necessary, add field names. Then set the values and geometry.

        '''
        # get geometry type
        geomType = self.layer.GetGeomType()

        # add pixel and line coordinates info to values
        if coordinates is not None and AddPixLine is True and \
           (geomType == ogr.wkbPoint or geomType == ogr.wkbPoint25D):
            # create values if values is None
            if values is None:
                values = np.zeros(len(coordinates[1]),
                                  dtype={'names':['X (pixel)', 'Y (line)'],
                                         'formats':['i8','i8']})
            # append pix/lin fields to values if values is given
            else:
                descr = [('X (pixel)', int), ('Y (line)', int)]
                # Create Zero structrued np.array ( n x m )
                # n is max( num of feature of values, num of coordinates)
                # m is num of field of values + 2 (for 'X (pixel)' and 'Y (line)' )
                tmp = np.zeros(max(values.shape[0], len(coordinates[1])),
                               dtype=values.dtype.descr + descr)
                # extend values
                for name in values.dtype.names:
                    # if num of coordinates is more than num of featuer of values,
                    # append 0 to values[name].
                    if len(coordinates[1]) > values.shape[0]:
                        tmp[name] = np.append(values[name],
                                              np.zeros(len(coordinates[1]) - values.shape[0]))
                    else:
                        tmp[name] = values[name]
                values = tmp
            # set coordinates into values
            for iCoord in range(len(coordinates[1])):
                values['X (pixel)'][iCoord] = coordinates[0][iCoord]
                values['Y (line)'][iCoord] = coordinates[1][iCoord]

        # if fieldFeatureID is not given, assume that fieldValues
        # are arranged from feature[0] to feature[n-1]
        if values is not None and fieldFeatureID == []:
                fieldFeatureID = range(len(values))

        # if geoFeatureID is None, create new one.
        if coordinates is not None and geoFeatureID == []:
            # if Point, featureID is serial number from 0 to len(coordinates)
            if geomType == ogr.wkbPoint or geomType == ogr.wkbPoint25D:
                geoFeatureID = range(len(coordinates[1]))
            # if line or ploygon, assume it is one geometry and belongs to the 1st feature.
            else:
                geoFeatureID = [0] * len(coordinates[1])

        # Remove duplicate featureID from geoFeatureID + fieldFeatureID
        featureID = list(set(geoFeatureID + fieldFeatureID))

        # convert list to np.array
        FieldFeatureID = np.asarray(fieldFeatureID)
        geoFeatureID = np.asarray(geoFeatureID)

        names = None
        # Tarnspose values
        if values is not None:
            # set self.fieldDataTypes of new fields
            self._make_datatype_dict(values)
            # set new fieldnames to the layer
            names = values.dtype.names
            self._add_fieldnames(names)

        for iFeature in featureID:
            # get values corresponding to iFeature
            if iFeature in fieldFeatureID:
                fieldValues = values[int(np.where(FieldFeatureID==iFeature)[0])]
            else:
                fieldValues = None

            # create a geometry
            geometry = None
            if coordinates is not None and iFeature in geoFeatureID:
                # get geoCoordinates corresponding to iFeature
                mask = geoFeatureID==iFeature
                geoCoords = []
                for iCoordinates in coordinates:
                    geoCoords.append(np.extract(mask, iCoordinates).tolist())
                geometry = ogr.Geometry(type=geomType)
                # if geomType is Point, set a point to a geometry
                if geomType == ogr.wkbPoint or geomType == ogr.wkbPoint25D:
                    geometry.SetPoint_2D(0, geoCoords[0][0], geoCoords[1][0])
                # else, set points to a geometry
                else:
                    geomRing = ogr.Geometry(type=ogr.wkbLinearRing)
                    for iPoint in range(len(geoCoords[1])):
                        # add points for a line
                        if geomType == ogr.wkbLineString or geomType == ogr.wkbLineString25D:
                            geometry.AddPoint(float(geoCoords[0][iPoint]), float(geoCoords[1][iPoint]))
                        # add points for a polygon
                        elif geomType == ogr.wkbPolygon or geomType == ogr.wkbPolygon25D:
                            geomRing.AddPoint(float(geoCoords[0][iPoint]), float(geoCoords[1][iPoint]))
                        else:
                            print 'Can not create geometry. Muptipolygon is not supported.'
                    # set geomRing if it has points (polygon)
                    if geomRing.GetPointCount() != 0:
                        geometry.AddGeometryDirectly(geomRing)
                # set srs
                srs =self.layer.GetSpatialRef()
                geometry.AssignSpatialReference(srs)

            # set field values and geometry to the feature
            self._set2feature(iFeature, fieldNames=names,
                              fieldValues=fieldValues, geometry=geometry)

    def _set2feature(self, featureID, fieldNames=[], fieldValues=None, geometry=None):
        ''' set fields or a geometry to a feature

        Parameters
        ----------
        featureID : int
            feature ID number
        fieldNames : list with string elements
            names of fields
        fieldValues : list (the lenght is a number of fields)
            values corresponding to the fieldNames list
        geometry : geometry object

        Modifies
        --------
        self.layer : set fields or a geometry
            if feature with given featureID has existed, set the values / geometry
            if not, create a new feature then set the values / geometry

        '''
        # get featuer with featureID. if it does not exist, create a new feature
        feature = self.layer.GetFeature(featureID)
        setFeature = True
        if feature is None:
            feature = ogr.Feature(self.layer.GetLayerDefn())
            setFeatuer = False

        # if geometry is given, set geometry
        if geometry is not None:
            feature.SetGeometryDirectly(geometry)

        # if fieldValue is given, set field values
        if fieldValues is not None:
            intType = [np.int, np.int8, np.int16, np.int32, np.int64,
                           np.uint16, np.uint32, np.uint64]
            floatType = [np.float, np.float16, np.float32, np.float64]
            for iField, iFieldName in enumerate(fieldNames):
                val = fieldValues[iField]
                if type(val) in intType:
                    val = int(val)
                elif type(val) in floatType:
                    val = float(val)
                feature.SetField(str(iFieldName), val)

        # Set the feature or add a new feature to the layer
        if setFeature:
            self.layer.SetFeature(feature)
        else:
            self.layer.CreateFeature(feature)
        feature.Destroy()

    def _add_fieldnames(self, fieldNames):
        ''' add a new field names to the layer

        Parameters
        ----------
        fieldNames : list with string elements
            field names

        Modifies
        --------
        self.layer : add new fields
        self.fieldNames : add new field names

        '''
        for iFieldName in fieldNames:
            if not (iFieldName in self.fieldNames):
                field = ogr.FieldDefn(iFieldName, self.fieldDataTypes[iFieldName])
                self.layer.CreateField(field)
                self.fieldNames.append(iFieldName)

    def _make_datatype_dict(self, values, names=None):
        '''set self.fieldDataTypes

        Parameters
        ----------
        values : structured numpy array or list with string or scalar elements
            values of fields
        name : list with string elements
            names of fields


        Modifies
        --------
        self.fieldDataTypes
            add new fields name and their OGR datatype

        '''
        intType = [np.int, np.int8, np.int16, np.int32, np.int64,
                       np.uint16, np.uint32, np.uint64, int]
        floatType = [np.float, np.float16, np.float32, np.float64, float]

        # if field names are in structured numpy array, set them in name list
        if names is None:
            names = values.dtype.names

        for i, iName in enumerate (names):
            if type(values) == np.ndarray:
                val = values.dtype.fields[iName][0]
            else:
                val = values[i][0]
            if not (iName in self.fieldDataTypes.keys()):
                # get ogr datatype for each field
                if type(val) in intType:
                    ogrDataType = ogr.OFTInteger
                elif type(val) in floatType:
                    ogrDataType = ogr.OFTReal
                else:
                    ogrDataType = ogr.OFTString
                # set field name and datatype to self.fieldDataTypes
                self.fieldDataTypes[iName] = ogrDataType

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
