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

        Nansatshape support points, but not line, ploygons, mupti-polygons

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
        wkbStyle : ogr.wkbPoint, ogr.wkbPoint25D

        Creates
        --------
        self.datasource : ogr data source in memory
        self.layer : ogr layer

        '''
        # create random name for the OGR dataset in memory
        allChars = ascii_uppercase + digits
        randomName = ''.join(choice(allChars) for x in range(10))

        # Create a empty datasource and layer in memory
        if fileName is None:
            self.datasource = ogr.GetDriverByName('Memory').CreateDataSource(randomName)
            # create a new later
            if layer == 0:
                layer = 'NansatLayer'
            self.layer = self.datasource.CreateLayer(layer, srs, wkbStyle)
        # if filename is given, open the file and copy the datasource in memory
        else:
            # Open shapefile and copy the datasource into memory
            ogrDs = ogr.Open(fileName)
            self.datasource = ogr.GetDriverByName('Memory').CopyDataSource(ogrDs, randomName)
            ogrDs.Destroy()
            # Set a layer from the datasource
            if type(layer) is int:
                self.layer = self.datasource.GetLayer(layer)
            else:
                self.layer = self.datasource.GetLayerByName(layer)

    def add_features(self, values, coordinates):
        ''' Set field values and / or geometry to each feature

        Loop over given arrays of coordinates and values and create
        corresponding points with data. Only ogr.wkbPoint is supported.

        Parameters
        ----------
        values :  2-D structured numpy array (aka record array)
            data to be stored in the vector layer. Names of the fields
            in the array will become fields in the layer.

        coordinates : np.array (n x 2)
            n rows with x and y coordinates of points. Length of
            coordinates should be equal to length of values.

        Modifies
        --------
        self.layer : Set the values and geometry.

        '''
        # get geometry type
        geomType = self.layer.GetGeomType()

        # return if input file is non-point
        if geomType != ogr.wkbPoint and geomType != ogr.wkbPoint25D:
            print 'Cannot add features to non-point layers'
            return

        # create fields from columns of values
        for i, iFieldName in enumerate(values.dtype.names):
            # get data type for each field
            if str(values.dtype[i]).startswith('int'):
                dtype = ogr.OFTInteger
            elif str(values.dtype[i]).startswith('float'):
                dtype = ogr.OFTReal
            else:
                dtype = ogr.OFTString
            # create field
            fieldDefn = ogr.FieldDefn(iFieldName, ogr.OFTString)
            fieldDefn.SetWidth(32)
            self.layer.CreateField(fieldDefn)

        # set values to each feature
        for iFeature in range(len(values)):
            feature = ogr.Feature(self.layer.GetLayerDefn())
            # create a geometry
            geometry = ogr.Geometry(type=geomType)
            # set a coordinates of the geometry
            geometry.SetPoint_2D(0, float(coordinates[0, iFeature]),
                                 float(coordinates[1, iFeature]))
            # set srs
            srs = self.layer.GetSpatialRef()
            geometry.AssignSpatialReference(srs)

            # set geometry
            feature.SetGeometryDirectly(geometry)

            # set field values
            for jField in values.dtype.names:
                if values.dtype[jField].name.startswith('int'):
                    feature.SetField(jField, int(values[iFeature][jField]))
                elif values.dtype[jField].name.startswith('float'):
                    feature.SetField(jField, float(values[iFeature][jField]))
                else:
                    feature.SetField(jField, str(values[iFeature][jField]))

            self.layer.CreateFeature(feature)
            feature.Destroy()

    def get_points(self, latlon=True):
        '''Get points (geometries of featuers) in the layer

        !!NB!!
        if shapefile has SRS, assume that geometry is lon/lat
        if not, assume that the geometry is given in pix/lin
        only ogr.wkbPoint or ogr.wkbPoint25D is supported

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
            if geom is not None:
                p = geom.GetPoints()[0]
                points.append((p[0], p[1]))

        if points == []:
            points = None
        else:
            points = tuple(points)
        return points, latlon

    def export(self, fileName):
        '''Save as ESRI shape-file'''
        shapeDriver = ogr.GetDriverByName("ESRI Shapefile")
        shapeDriver.CopyDataSource(self.datasource, fileName)
