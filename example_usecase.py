#-------------------------------------------------------------------------------
# Name:      example_usecase.py
# Purpose:   Example that shows how to use NANSAT.
#            Calls nansat module that is main of nansat.
#
# Author:      asumak
#
# Created:     29.06.2011
# Copyright:   (c) asumak 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os

try:
    from osgeo import gdal
except ImportError:
    import gdal

try:
    from osgeo.gdalconst import *
except ImportError:
    from gdalconst import *

import matplotlib.pyplot as plt
from scipy.misc import toimage

from nansat import *
'''
ifolder = "/Home/asumak/nansat/data/input/" # Asuka
ofolder = "/Home/asumak/nansat/data/output/" # Asuka
import gdalinfo
#- Uncomment/comment lines below to work with a specific file
##file = "MYD02QKM.A2011182.0040.005.2011182193302.hdf"
#file = "MOD021KM.A2011121.1105.005.2011121193140.hdf"
#file = "MOD02QKM.A2011229.1125.005.2011229195243.hdf"
#file = "MYD02HKM.A2009105.0915.005.2009105193508.hdf"
#file = "MOD021KM.A2010105.2120.005.2010106075131.hdf"
#file = "MYD021KM.A2011228.0925.005.2011229003113.hdf"
##file = "MER_FRS_1PNRAL20100426_091349_000000822088_00494_42632_0001.N1"
#file = "RS2_20090227_063055_0080_SCWA_HHHV_SCW_30853_0000_1897838"
#file = "RS2_20111109_060616_0045_SCNA_HHHV_SGF_164373_9871_6913894"
#file = "gfs.t06z.master.grbf03"
#file = "DNMI-NEurope.grb"
file = "ASA_WSM_1PNPDK20110108_205958_000000923098_00187_46322_6032.N1"
#file = "MER_RR__2PNPDK20100519_101302_000005462089_00323_42962_4158.N1"
#file = "MER_RR__2PNPDK20110921_112637_000025653106_00411_49994_3048.N1"
#file = "MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1"
#file = "MER_FRS_1PNPDK20110503_105638_000001833102_00109_47968_7898.N1"
#file = "MOD02HKM.A2011121.1105.005.2011121193140.hdf"
#file = "MOD02QKM.A2011121.1105.005.2011121193140.hdf"
#file = "MOD021KM.A2011121.1105.005.2011121193140.hdf"

'''
ifolder = "c:/Users/asumak/Data/input/" # Asuka
ofolder = "c:/Users/asumak/Data/output/" # Asuka

from osgeo import gdalinfo
#- Uncomment/comment lines below to work with a specific file
#file = "Aqua/ModisL1/250m/MYD02QKM.A2011182.0040.005.2011182193302.hdf"
#file = "Aqua/ModisL1/250m/MOD02QKM.A2011229.1125.005.2011229195243.hdf"
#file = "Aqua/ModisL1/250m/MOD02QKM.A2011121.1105.005.2011121193140.hdf"
#file = "Aqua/ModisL1/500m/MYD02HKM.A2009105.0915.005.2009105193508.hdf"
#file = "Aqua/ModisL1/1km/MOD021KM.A2010105.2120.005.2010106075131.hdf"
file = "Aqua/ModisL1/1km/MYD021KM.A2011228.0925.005.2011229003113.hdf"
#file = "ENVISAT/MERIS/Level2/MER_RR__2PNPDK20100519_101302_000005462089_00323_42962_4158.N1"
#file = "ENVISAT/MERIS/Level1/MER_FRS_1PNRAL20100426_091349_000000822088_00494_42632_0001.N1"
#file = "ENVISAT/ASAR/ASA_WSM_1PNPDK20110108_205958_000000923098_00187_46322_6032.N1"
#file = "RS2_20090227_063055_0080_SCWA_HHHV_SCW_30853_0000_1897838"
#file = "RS2_20111109_060616_0045_SCNA_HHHV_SGF_164373_9871_6913894"
#file = "NCEP_GRIB/gfs.t06z.master.grbf03"
#file = "HIRLAM_GRIB/DNMI-NEurope.grb"
#file = "MOD021KM.A2011121.1105.005.2011121193140.hdf"
#file = "MOD02HKM.A2011121.1105.005.2011121193140.hdf"
#file = "MOD02QKM.A2011121.1105.005.2011121193140.hdf"


# By setting an environmental variable "NANSAT_TEST_FOLDER", this script can be used unmodified also by Anton and Knut-Frode
# This requires in addition that the images are also stored in subfolders as suggested by Asuka (e.g. "ENVISAT/MERIS" and "NCEP_GRIB")

try:
    folder = os.environ["NANSAT_TEST_FOLDER"]
except:
    print "Environmental variable NANSAT_TEST_FOLDER is not set, using " + ifolder

fileName = ifolder + file;

#fileNameProj = "c:/Users/asumak/Data/input/RS2_20090227_063055_0080_SCWA_HHHV_SCW_30853_0000_1897838"
#gdalinfo.main(['foo', fileName])

#------------------------------------------------------------------------------#
# Open the file
#------------------------------------------------------------------------------#
#mapperName = "radarsat1";
#bandList = [1,2];
#n = Nansat(fileName, mapperName, bandList);

n = Nansat(fileName);

#------------------------------------------------------------------------------#
# show serial bandNo. , name and parameters
#------------------------------------------------------------------------------#
n.list_bands();

#------------------------------------------------------------------------------#
# reduce the size of an image
#------------------------------------------------------------------------------#
# -- "method" should be 'subsampling' (nearest neighbor) or
#                     'blockaveraging' (averaged resampling)
#n.downscale(factor = 10, method = "subsampling");

#------------------------------------------------------------------------------#
# get a band as a GDALRasterBand
#------------------------------------------------------------------------------#
# -- based on band No. (default is bandNo. = 1)
band = n.get_GDALRasterBand(1);

# -- based on ShortName + parameters
#    this is prior to array = n.get_band(bandNo)
#bandIdDic = {"name":"radiance", "wavelength":"13935"};
#bandIdDic = {"name":"radiance"};
#band = n.get_GDALRasterBand(bandID = bandIdDic);

#------------------------------------------------------------------------------#
# Save to a raster band to a figure in PNG
#------------------------------------------------------------------------------#
# -- based on band No. (default is bandNo. = 1)
#n.write_figure(ofolder + "writeFigure.png", 1, thresholdRatio = 0.99);

# -- based on ShortName + parameters
#bandIdList = {"wkv":"radiance", "wavelength":"645"};
#n.write_figure("c:/Users/asumak/Data/output/testRasterBand.png", bandName = bandIdList);

#------------------------------------------------------------------------------#
# export in-memory VRT dataset to a physical file
#------------------------------------------------------------------------------#
n.export_VRT(ofolder + "source_vrt.vrt");
#n.export_VRT();

#------------------------------------------------------------------------------#
# Reprojection
#------------------------------------------------------------------------------#

#n.reproject();

# -- without option
#srsString = "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs";
srsString = "+proj=latlong +datum=WGS84 +no_defs"
#srsString = "+proj=utm +zone=33 +datum=WGS84 +ellps=WGS84 +no_defs";
#srsString = "+proj=utm +zone=33 +datum=WGS84 +no_defs";
#srsString = "+proj=utm +datum=WGS84 +ellps=WGS84 +no_defs";
#srsString = "+proj=stere +lon_0=-25.0 +lat_0=0.0 +k=0.9999079 +no_defs";
#srsString = "+proj=stere +lon_0=-25.0 +lat_0=0.0 +k=0.9999079 +no_defs";
n.reproject(srsString);

# -- with option
# Requirement :  ("-lle" or "-te") and ("-ts" or "-tr")
# "-lle lonwest loneast latnorth latsouth"
# "-te xmin ymin xmax ymax"
#extentString = "-te -1188742 7317234 2301290 9829086 -tr 1163 1004"
#extentString = "-ts 500 500"
#extentString = "-lle 15 20 80 75 -ts 500 500"
#extentString = "-te 212366 8323606 500000 8898211 -tr 1000 1000"
#extentString = "-lle 25.0 30.0 80.0 75.0"
#extentString = "-tr 3500 3500"
#extentString = "-te -15000000 0 20000000 10000000 -tr 10000 10000"
#extentString = "-lle -10 30 55 60.0 -tr 1000 1000"
# ----- resamplingAlg
#       0 = "near", 1 = "bilinear", 2 = "cubic", 3 ="cubic spline", 4 ="lanczos"
#n.reproject(srsString, extentString, resamplingAlg = 0);

#n.reproject();
#band = n.get_GDALRasterBand(1);
#band = n.get_GDALRasterBand(bandNo);

#------------------------------------------------------------------------------#
# Save to a raster band to a figure in PNG
#------------------------------------------------------------------------------#
# -- based on band No. (default is bandNo. = 1)
#n.export_VRT(ofolder + "/warped_vrt.vrt");
#n.write_figure(ofolder + "writeFigure.png", 1, thresholdRatio = 0.99);
#n.write_figure(ofolder + "writeFigure.png", 1, pixValMin = None, pixValMax = None, imageDatatype = np.uint8, thresholdRatio = 0.90);
#n.export_VRT("c:/Users/asumak/Data/output/example.vrt");

#n.write_figure("c:/Users/asumak/Data/output/testRasterBand.png", 1, thresholdRatio = 0.99);

# -- based on ShortName + parameters
#bandIdList = {"polarization":"HV"};
#n.write_figure("c:/Users/asumak/Data/output/testRasterBand", bandName = bandIdList);

#------------------------------------------------------------------------------#
# Save to raster bands to a Tiff file
#------------------------------------------------------------------------------#
# elements in the list is given by n.list_bands()
#bands = [2];
#dataType = gdal.GDT_Int16 # default
#dataType = gdal.GDT_Float32
#n.export(ofolder + "testExport", bands, dataType);














