#-------------------------------------------------------------------------------
# Name:        test.pyat
# Purpose:
#
# Author:      asumak
#
# Created:     04.10.2011
# Copyright:   (c) asumak 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

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

#in_folder = "c:/Users/asumak/Data/input" # Asuka
#out_folder = "c:/Users/asumak/Data/output"

in_folder = "/Home/asumak/nansat/data/input" # Asuka
out_folder = "/Home/asumak/nansat/data/output" # Asuka

# By setting an environmental variable "NANSAT_TEST_FOLDER", this script can be used unmodified also by Anton and Knut-Frode
# This requires in addition that the images are also stored in subfolders as suggested by Asuka (e.g. "ENVISAT/MERIS" and "NCEP_GRIB")
try:
    in_folder = os.environ['NANSAT_TEST_FOLDER']
    out_folder = os.environ['NANSAT_TEST_FOLDER'] + '/output/'
except:
    print "Environmental variable NANSAT_TEST_FOLDER is not set, using " + in_folder

satList = [
"DNMI-NEurope.grb",
"gfs.t06z.master.grbf03",
"MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1",
"MER_FRS_1PNPDK20110503_105638_000001833102_00109_47968_7898.N1",
"MOD02HKM.A2011121.1105.005.2011121193140.hdf",
"MOD02QKM.A2011121.1105.005.2011121193140.hdf",
"MOD021KM.A2011121.1105.005.2011121193140.hdf",
"ASA_WSM_1PNPDK20110108_205958_000000923098_00187_46322_6032.N1",
"RS2_20111109_060616_0045_SCNA_HHHV_SGF_164373_9871_6913894"
#"MER_RR__2PNPDK20110921_112637_000025653106_00411_49994_3048.N1"
]
'''
satList = [
"DNMI-NEurope.grb",
"gfs.t06z.master.grbf03"
]
'''
srsString = "+proj=stereo +lon_0=10.5 +lat_0=65.0 +k=0.9999079 +no_defs";
#srsString = "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs";
'''
# MDIS
satList = [
"MOD021KM.A2011121.1105.005.2011121193140.hdf"
]
#srsString = "+proj=stereo +lon_0=5.0 +lat_0=65.0 +k=0.9999079 +no_defs";
srsString = "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs";
'''
'''
# RS2, NCEP_GRIB, HIRLAM_GRIB
satList = [
"ASA_WSM_1PNPDK20110108_205958_000000923098_00187_46322_6032.N1",
"RS2_20090227_063055_0080_SCWA_HHHV_SCW_30853_0000_1897838",
"gfs.t06z.master.grbf03",
"DNMI-NEurope.grb"
]
srsString = "+proj=utm +zone=25 +datum=WGS84 +no_defs";
'''

extentString = "-lle 9.0 12.0 66.0 64.0 -tr 250 250"

for iFile in range(len(satList)):
    fileName = in_folder + "/" + satList[iFile];

#   Open the file
    n = Nansat(fileName);

#   Show serial bandNo. , name and parameters
    n.list_bands();

#   Save a raster band to a figure in PNG
    outPath = satList[iFile].rsplit(".", 1);
    '''
    print "--Raw--"
    n.write_figure(out_folder + "/" + outPath[0] + "_raw.png", 1, imageSize = 0.1);
#   export in-memory VRT dataset to a physical file
    n.export_VRT(out_folder + "/" + outPath[0] + "_raw.vrt");

#   reproject with default option
    print "--Default--"
    n.reproject();
    n.write_figure(out_folder + "/" + outPath[0] + "_proj0.png", 1, imageSize = 1.0);
#   export in-memory VRT dataset to a physical file
    n.export_VRT(out_folder + "/" + outPath[0] + "_proj0.vrt");

#   reproject with option1
    print "--Option1--"
    n.reproject(srsString);
    n.write_figure(out_folder + "/" + outPath[0] + "_proj1.png", 1, imageSize = 0.1);
#   export in-memory VRT dataset to a physical file
    n.export_VRT(out_folder + "/" + outPath[0] + "_proj1.vrt");
    '''
#   reproject with option2
    print "--Option2--"
    # resamplingAlg
    # 0 = "near", 1 = "bilinear", 2 = "cubic", 3 ="cubic spline", 4 ="lanczos"
    n.reproject(srsString, extentString, resamplingAlg = 3);
    n.write_figure(out_folder + "/" + outPath[0] + "_latlong.png", 1, pixValMin = None, pixValMax = None, imageDatatype = np.uint8, thresholdRatio = 0.700);
    #n.write_figure(out_folder + "/" + outPath[0] + "_proj2.png", 1, pixValMin = None, pixValMax = None);
#   export in-memory VRT dataset to a physical file
    #n.export_VRT(out_folder + "/" + outPath[0] + "_proj2.vrt");
