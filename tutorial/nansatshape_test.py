import inspect, os
import numpy as np

from nansat import Nansat, Domain, Nansatshape

try:
    from osgeo import gdal, osr, ogr
except:
    import gdal, osr, ogr

# input and output file names
iPath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
iFileName = os.path.join(iPath, 'points.shp')
print 'Input file: ', iFileName
oFileName = os.path.join(iPath, 'tmpdata', 'outputmap_')
print 'Output file prefix: ', oFileName


# Create nansatShape object with point geometry
nansatOGR = Nansatshape(wkbStyle=ogr.wkbPoint)
# create point geometries and set them to features
nansatOGR.create_geometry([[10.0, 15.0, 20.0, 25.0, 30.0],[5.0, 45.0, 25.0, 35.0, 10.0]])
# add fields to the feature
nansatOGR.create_fields(fieldNames=['int', 'string', 'float'],
                        fieldValues=[[100, 200, 300, 400, 500], ['S0','S1','S2','S3','S4'], [1.1, 1.2, 1.3, 1.4, 1.5]])
# append a new field to specified features
nansatOGR.create_fields(fieldNames=['string2'], fieldValues=[['SS1','SS3']], featureID=[1,3])
# append a new feature and set a new geometry
nansatOGR.create_geometry([[80.0],[20.0]], featureID=[5])
# set field values
nansatOGR.create_fields(fieldNames=['int', 'string2'], fieldValues=[[60],['SS5']], featureID=[5])
# save to a file
ogr.GetDriverByName("ESRI Shapefile").CopyDataSource(nansatOGR.datasource, oFileName+'01_Points.shp')

# Create a nansatShape object with polygon geometry
nansatOGR = Nansatshape(wkbStyle=ogr.wkbPolygon)
# create polygon geometry and set them to featuers
nansatOGR.create_geometry([[100, 150, 200, 100, 500, 600, 800, 500],[20, 30, 60, 20, 30, 60, 90, 30]],
                           featureID=[0, 0, 0, 0, 1, 1, 1, 1])
# append a polygon geometry to a new feature
nansatOGR.create_geometry([[300,700, 750, 300],[80,50, 60, 80]], featureID=[2,2,2,2])
# save to a file
ogr.GetDriverByName("ESRI Shapefile").CopyDataSource(nansatOGR.datasource, oFileName+'02_Polygons.shp')

# Create a nansatShape object with polygon geometry
nansatOGR = Nansatshape(wkbStyle=ogr.wkbLineString)
# create polygon geometry and set them to featuers
nansatOGR.create_geometry([[100, 150, 200, 500, 600, 800],[20, 30, 60, 30, 60, 90]],
                           featureID=[0, 0, 0, 1, 1, 1])
# append a polygon geometry to a new feature
nansatOGR.create_geometry([[300,700, 750],[80,50, 60]], featureID=[2,2,2])
# add fields to the feature
nansatOGR.create_fields(fieldNames=['int', 'string', 'float', 'string2'], fieldValues=[[100, 200, 300], ['S0','S1','S2'],[1.1, 2.2, 3.3], ['A0','A1','A2'] ])
# replace specified field
nansatOGR.create_fields(fieldNames=['int'], fieldValues=[[500]], featureID=[1])
# save to a file
ogr.GetDriverByName("ESRI Shapefile").CopyDataSource(nansatOGR.datasource, oFileName+'03_Lines.shp')

# Create a nansatShape from shape file
nansatOGR = Nansatshape(iFileName)
# Get corner points (geometries of featuers) in the layer
points, latlon = nansatOGR.get_corner_points(latlon=False)
# print corner points
print '-- Corner points --'
for iPoint in points:
    print iPoint








