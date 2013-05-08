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
    def __init__(self, fileName=None, layerNum=0, layerName='NansatLayer', srs=None, wkt=None, wkbStyle=ogr.wkbPoint):
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
        layerNum : int
            if a shapefile is given, it is a layer number which is read.
        layerName : String
            layer name to created or open
        srs : SpatialReference object
        wkt : String
            Well known text
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
            # create srs if wkt is given
            if srs is None and wkt is not None:
                srs = osr.SpatialReference()
                srs.ImportFromWkt(wkt)
            # create a new later
            self.layer = self.datasource.CreateLayer(layerName, srs, wkbStyle)
        # if filename is given, open the file and copy the datasource in memory
        else:
            # Open shapefile and copy the datasource into memory
            ogrDs = ogr.Open(fileName)
            self.datasource= ogr.GetDriverByName('Memory').CopyDataSource(ogrDs, 'OGRDataSource')
            ogrDs.Destroy()
            # Set a layer from the datasource
            if layerName == 'NansatLayer':
                self.layer = self.datasource.GetLayer(layerNum)
            else:
                self.layer = self.datasource.GetLayerByName(layerName)
            # Set self.fieldNames and self.fieldDataTypes
            layerDef = self.layer.GetLayerDefn()
            feature = self.layer.GetFeature(0)
            for i in range(layerDef.GetFieldCount()):
                self.fieldNames.append(layerDef.GetFieldDefn(i).GetName())
                self.fieldDataTypes[self.fieldNames[-1]] = feature.GetFieldType(i)

    def set_layer(self, lonlatCoord=None,  pixlinCoord=None, fieldNames=[], fieldValues=None):
        ''' set geometry and / or fields to each feature in self.layer

        Set geometry and fields to featues.
        If it is necessary, new features and fields are automatically added.

        Parameters
        ----------
        lonlatCoord : np.array with scalar elements
            coordinates in degree.
            the size of the array is (2, n), here n is number of coordinates.
        pixlinCoord :  np.array with scalar elements
            pixlinCoord in column/row.
            the size of the array is (2, n), here n is number of coordinates.
        fieldNames : list with string elements
            names of fields
        fieldValues :  list or np.array
            the size of the array is (m, n),
            n is a number of features and m is a number of fields.

        Modifies
        ---------
        self.layer : set geometry and/or fields to featues.

        '''
        # get data type of each field
        if fieldNames != []:
            if type(fieldValues[0]) != list:
                fieldValues = [fieldValues]
            self._make_datatype_dict(fieldNames, fieldValues)

        # add pixel and line coordinates info to fieldNames and fieldValues
        if pixlinCoord is not None:
            # add filedName for pixel and line coordinates
            for iFieldName in ['X (pixel)', 'Y (line)']:
                fieldNames.append(iFieldName)
                self.fieldDataTypes[iFieldName] = ogr.OFTReal
            # add filedValue of pixel and line coordinates
            if fieldValues is None:
                fieldValues = [pixlinCoord[0].tolist(), pixlinCoord[1].tolist()]
            else:
                fieldValues.append(pixlinCoord[0].tolist())
                fieldValues.append(pixlinCoord[1].tolist())

        # set new field names and datatypes to the layer
        self._add_fieldnames(fieldNames)

        # set geometry
        if lonlatCoord is not None:
            self.create_geometry(lonlatCoord)
        elif pixlinCoord is not None:
            self.create_geometry(pixlinCoord)

        # set field values
        self.create_fields(fieldNames, fieldValues)

    def create_geometry(self, coordinates, featureID=None):
        ''' create geometry objects

        !! NB !! Muptipolygon is not supported

        Parameters
        ----------
        coordinates : np.array or list with scalar elements (2 x n)
            n is a number of points
        feateurID : None or list with int elements (the length is n)
            if the data type is point, it will be [0, 1, 2, 3, ..., (num of fields-1)]
            if polygon or line, the numbers represent which geometry each point belongs to.
        srs : osr.SpatialReference

        '''
        # if list is given, conver to numpy array
        if type(coordinates) == list:
            coordinates = np.asarray(coordinates)

        # get geometry type
        geomType = self.layer.GetGeomType()
        # if featureID is None, create new one.
        if featureID == None:
            if geomType == ogr.wkbPoint or geomType == ogr.wkbPoint25D:
                featureID = range(coordinates.shape[1])
            else:
                featureID = [0] * coordinates.shape[1]

        # Remove duplicate feature ID from featureID
        seen = set()
        seen_add = seen.add
        featureNum = [x for x in featureID if x not in seen and not  seen_add(x)]

        # create a geometry
        for i, iFeature in enumerate(featureNum):
            geometry = ogr.Geometry(type=geomType)
            # if geomType is Point, set a point to a geometry
            if geomType == ogr.wkbPoint or geomType == ogr.wkbPoint25D:
                geometry.SetPoint_2D(0, coordinates[0][i], coordinates[1][i])
            # else, set points to a geometry
            else:
                geomRing = ogr.Geometry(type=ogr.wkbLinearRing)
                for iCoord in range(coordinates.shape[1]):
                    # add points for a line
                    if geomType == ogr.wkbLineString or geomType == ogr.wkbLineString25D:
                        if iFeature == featureID[iCoord]:
                            geometry.AddPoint(float(coordinates[0][iCoord]), float(coordinates[1][iCoord]))
                    # add points for a polygon
                    elif geomType == ogr.wkbPolygon or geomType == ogr.wkbPolygon25D:
                        if iFeature == featureID[iCoord]:
                            geomRing.AddPoint(float(coordinates[0][iCoord]), float(coordinates[1][iCoord]))
                    else:
                        print 'Can not create geometry. Muptipolygon is not supported.'
                # set geomRing if it has points (polygon)
                if geomRing.GetPointCount() != 0:
                    geometry.AddGeometryDirectly(geomRing)
            # set srs
            srs =self.layer.GetSpatialRef()
            geometry.AssignSpatialReference(srs)

            # set geometry to each feature
            self._set2feature(iFeature, geometry=geometry)

    def create_fields(self, fieldNames=[], fieldValues=[], featureID=None):
        ''' create fields and set the values

        Parameters
        ----------
        fieldNames : list with string elements (the lenght is m)
            names of fields
        fieldValues :  np.array (m x n) or list with one or more list as the element.
            e.g. : [[A1, A2, ..., An],[B1, B2, ..., Bn], ...,[Z1, Z2, .., Zn]]
            n is a number of features and m is a number of fields.
        feateurID : None or list with int elements (the length is n)
            the values are from 0 to n-1

        Modifies
        --------
        self.layer : add field names and set the values

        '''
        # set self.fieldDataTypes of new fields
        self._make_datatype_dict(fieldNames, fieldValues)

        # Tarnspose vieldValues
        fieldValues = map(list, zip(*fieldValues))

        # Assume that fieldValues are arranged from feature[0] to feature[n-1]
        if featureID == None:
            featureID = range(len(fieldValues))

        # set new fieldnames to the layer
        self._add_fieldnames(fieldNames)

        # set field values to each featuer
        for i, iFeature in enumerate (featureID):
            self._set2feature(iFeature, fieldNames=fieldNames, fieldValues=fieldValues[i])

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

        # if geometry is given, set geometry. otherwise, set field values
        if geometry is not None:
            feature.SetGeometryDirectly(geometry)
        else:
            for iField, iFieldName in enumerate(fieldNames):
                feature.SetField(str(iFieldName), fieldValues[iField])

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

    def _make_datatype_dict(self, names, values):
        '''set self.fieldDataTypes

        Parameters
        ----------
        name : list with string elements
            names of fields
        values : list with string or scalar elements
            values of fields

        Modifies
        --------
        self.fieldDataTypes
            add new fields name and their OGR datatype

        '''
        for i, iName in enumerate (names):
            if not (iName in self.fieldDataTypes.keys()):
                ogrDataType = {
                    "<type 'str'>" : ogr.OFTString,
                    "<type 'int'>" : ogr.OFTInteger,
                    "<type 'int32'>" : ogr.OFTInteger,
                    "<type 'float'>": ogr.OFTReal,
                    }.get(str(type(values[i][0])), ogr.OFTString)
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
