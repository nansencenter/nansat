#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      asumak
#
# Created:     22.06.2012
# Copyright:   (c) asumak 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import os
from os.path import join, getsize

try:
    from osgeo import gdal
except ImportError:
    import gdal

try:
    from osgeo.gdalconst import *
except ImportError:
    from gdalconst import *

import matplotlib.pyplot as plt
import string
from scipy.misc import toimage
from scipy import ndimage

from nansat import *


ifolder = "/Home/asumak/data/input/" # Asuka
ofolder = "/Home/asumak/data/output/" # Asuka
#import gdalinfo

iFileName = "MOD021KM.A2011121.1105.005.2011121193140.hdf"
#iFileName = "ASA_WSM_1PNPDK20110108_205958_000000923098_00187_46322_6032.N1"
#iFileName = "MER_FRS_1PNPDK20110503_105638_000001833102_00109_47968_7898.N1"
#iFileName = "asi-n6250-20110215-v5.tif"
#iFileName = "TSX1_SAR__SSC______SL_S_SRA_20071215T112105_20071215T112107"
#iFileName = "HIRLAM_GRIB/DNMI-NEurope_proj2.vrt"
#iFileName = "netCDF/TP6AVEDAILYNPP200609.nc"
#iFileName = "netCDF/westernnorway300m_oceancolor_20110301-20110331.nc"
#iFileName = "netCDF/2005092200_sst_21-24.en.nc"

def create_GeolocationBand(vrt, longitude=None, latitude=None, subDatasets=[] ,logLevel=30):
    # === ADD lat/lon bands ====
    # if shape of lat/lon is different from other bands
    dsShape = (vrt.dataset.RasterYSize, vrt.dataset.RasterXSize)

    if longitude is not None and latitude is not None:
        if longitude.shape != dsShape:
            print "ZOOM"
            zoomFactor = (float(dsShape[0]) / float(longitude.shape[0]),
                          float(dsShape[1]) / float(longitude.shape[1]))
            # create RAW/VRT files with full size lat/lon
            longitude = ndimage.zoom(longitude, zoomFactor)
            latitude = ndimage.zoom(latitude, zoomFactor)
        print "Create VRT"
        vrt.logger.debug('Lat/Lon grids [%f, %f] read', latitude.shape[0], latitude.shape[1])
        vrt.lonVrt = VRT(array=longitude, logLevel=logLevel)
        vrt.latVrt = VRT(array=latitude, logLevel=logLevel)
        # add bands pointing to RAW/VRT
        vrt._create_band(vrt.lonVrt.fileName, 1, 'longitude', {'band_name': 'longitude'})
        vrt._create_band(vrt.latVrt.fileName, 1, 'latitude',  {'band_name': 'latitude'})
    else:
        print "Use original files"
        # else just take lat/lons from original file
        for subDataset in subDatasets:
            if 'longitude' in subDataset[1]:
                vrt._create_band(subDataset[0], 1, 'longitude', {'band_name': 'longitude'})
            if 'latitude' in subDataset[1]:
                vrt._create_band(subDataset[0], 1, 'latitude',  {'band_name': 'latitude'})


def setGeolocation(nObj):
    NonValue = 50000
    # set GeoTransform
    trans = nObj.vrt.dataset.GetGeoTransform()
    transString = str(trans[0]) + " " +str(trans[1]) + " "+str(trans[2]) + " "+str(trans[3]) + " "+str(trans[4]) + " "+str(trans[5])
    nObj.set_metadata("GeoTarnsform", transString)

    # modify projection string
    projection = nObj.raw.dataset.GetProjection()
    if projection == '':
        print "from GCP projection"
        projection = nObj.raw.dataset.GetGCPProjection()

    # if GEOLOCATION is given, set metadata and create VRT
    if nObj.vrt.dataset.GetMetadata('GEOLOCATION') != {}:
        print "from GEOLOCATION"
        # set GEOLOCATION metadata
        geolocDict = nObj.vrt.dataset.GetMetadata('GEOLOCATION')
        print projection
        geolocDict['SRS'] = projection
        nObj.set_metadata(geolocDict)

        # create array from GeolocaionArray
        geoFileNames = [geolocDict['X_DATASET'], geolocDict['Y_DATASET']]
        #set longitude and latitude in the arrays
        arrayList = []
        for fileName in geoFileNames:
            geoDs = gdal.Open(fileName)
            arrayList.append(geoDs.GetRasterBand(1).ReadAsArray())
        # add bands
        subDatasets = [(geoFileNames[0],"longitude"),(geoFileNames[1],"latitude")]
        create_GeolocationBand(nObj.vrt, arrayList[0], arrayList[1], subDatasets)

    # if GCPs are given, create VRT
    elif nObj.vrt.dataset.GetGCPs() is not ():
        print "from GCPs"
        wkvList = ["X_GCPs", "Y_GCPs"]
        arrayList = [] # arrayList[0] for longitude, arrayList[1] for latitude
        for iList in range(2):
            arrayList.append(np.ones((nObj.vrt.dataset.RasterYSize, nObj.vrt.dataset.RasterXSize)) * NonValue)
        # set gcps in the arrays
        gcps = nObj.vrt.dataset.GetGCPs()
        for igcp in gcps:
            try:
                arrayList[0][int(igcp.GCPLine)][int(igcp.GCPPixel)] = ((igcp.GCPX + 180) * 100.0)
                arrayList[1][int(igcp.GCPLine)][int(igcp.GCPPixel)] = ((igcp.GCPY + 180) * 100.0)
            except:
                print ("GCPLine or GCPPixel is over the raster size")
                #self.logger.warning("GCPLine or GCPPixel is over the raster size")

        for i, iWKV in enumerate(wkvList):
            nObj.add_band(array=np.array(arrayList[i], dtype='uint16'), p={"band_name":iWKV, "missing_value":str(NonValue)})
        nObj.set_metadata("gcpProjection", projection)

"""=========================================================================="""

#ds = gdal.Open(ifolder + iFileName)
#ds = gdal.GetDriverByName("VRT").CreateCopy(ofolder+"gdal.vrt", ds)

n1 = Nansat(ifolder + iFileName, logLevel=30)
#n1.vrt.export(ofolder+"n0.vrt")
setGeolocation(n1)

# Create a netCDF file
#print "export to VRT..."
#n1.vrt.export(ofolder+"n1.vrt")
print "export to netCDF..."
n1.export(ofolder+"test.nc", rmMetadata=['band_name'])
#band = n1.vrt.dataset.GetRasterBand(1)
#band.SetMetadataItem("band_name", "BandNameTest")
n1.export(ofolder+"test.nc")

# Open the netCDF file
#print "Open n2..."
n2 = Nansat(ofolder+"test.nc")
#n2 = Nansat(ifolder + iFileName)
#n2.write_figure(ofolder + "fig", bands=3)
n2.vrt.export(ofolder + "n2.vrt")
#band = n2.vrt.dataset.GetRasterBand(1)
#print band.GetMetadata()

