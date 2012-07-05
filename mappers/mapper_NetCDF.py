#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      asumak
#
# Created:     29.06.2012
# Copyright:   (c) asumak 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from vrt import *
from nansat_tools import Node
import numpy as np

class Mapper(VRT):
    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30):

        if not "NC_GLOBAL#Conventions" in gdalMetadata.keys():
            raise AttributeError("NETCDF BAD MAPPER")

        # Get file names from dataset or subdataset
        if gdalDataset.RasterCount==1:
            fileNames = {"SUBDATASET_1_NAME" : fileName}
        else:
            fileNames = gdalDataset.GetMetadata('SUBDATASETS')
            for iKey in fileNames.keys():
                if not '_NAME' in iKey:
                    del fileNames[iKey]

        # get raster size from the first band
        firstSubDataset = gdal.Open(fileNames["SUBDATASET_1_NAME"])
        bandSize = {"XSize" : firstSubDataset.RasterXSize,
                    "YSize" : firstSubDataset.RasterYSize}

        # add bands with metadata and corresponding values to the empty VRT
        metaDictList = []
        geoFileDict = {}
        for iFileName in fileNames:
            subDataset = gdal.Open(fileNames[iFileName])
            # take bands whose sizes are same as the first band.
            if subDataset.RasterXSize==bandSize["XSize"] and subDataset.RasterYSize==bandSize["YSize"]:
                for iBand in range(subDataset.RasterCount):
                    bandDict = subDataset.GetRasterBand(iBand+1).GetMetadata_Dict()
                    try:
                        source = bandDict.pop("source")
                    except:
                        source = fileNames[iFileName]
                    try:
                        sourceBands = int(bandDict.pop("sourceBands"))
                    except:
                        sourceBands = iBand + 1
                    try:
                        varName = bandDict.pop("NETCDF_VARNAME")
                        if varName.find("X_GCPs") != -1:
                            geoFileDict["X_GCPs"] = source.replace(".vrt", ".raw")
                        elif varName.find("Y_GCPs") != -1:
                            geoFileDict["Y_GCPs"] = source.replace(".vrt", ".raw")
                    except:
                        pass
                    try:
                        bandDict.pop("_FillValue")
                    except:
                        pass
                    metaDictList.append(({'source': source, 'sourceBand': sourceBands, 'wkv': "",'parameters':bandDict}))

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDataset, logLevel=logLevel)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDictList, bandSize)

        # Modify band size for dataset
        VRTXML = self.read_xml()
        node0 = Node.create(VRTXML)
        node0.replaceAttribute("rasterXSize", str(bandSize["XSize"]))
        node0.replaceAttribute("rasterYSize", str(bandSize["YSize"]))
        self.write_xml(str(node0.rawxml()))

        # set projection to dataset
        if 'NC_GLOBAL#GDAL_gcpProjection' in gdalMetadata.keys():
            projection = gdalMetadata['NC_GLOBAL#GDAL_gcpProjection']
            projection = projection.replace("|",",").replace("&quat;",'"')
            self.dataset.SetProjection(projection)
            if 'NC_GLOBAL#GDAL_gcpProjection' in gdalMetadata.keys():
                self.dataset.SetMetadataItem('NC_GLOBAL#GDAL_gcpProjection', projection)
            if 'NC_GLOBAL#GDAL_projection' in gdalMetadata.keys():
                self.dataset.SetMetadataItem('NC_GLOBAL#GDAL_projection', projection)
            if 'NC_GLOBAL#GDAL_SRS' in gdalMetadata.keys():
                self.dataset.SetMetadataItem('NC_GLOBAL#GDAL_SRS', projection)

        # for netCDF created by NANSAT
        if "X_GCPs" in geoFileDict and "Y_GCPs" in geoFileDict:
            # set GCPs
            print "GCPs"
            NonValue = 50000
            fileNameList = [geoFileDict["X_GCPs"], geoFileDict["Y_GCPs"]]
            arrayList = []
            for iFile in fileNameList:
                f = gdal.VSIFOpenL(iFile, "rb")
                gdal.VSIFSeekL(f, 0, 2)
                vsiFileSize = gdal.VSIFTellL(f)
                gdal.VSIFSeekL(f, 0, 0)
                strData = gdal.VSIFReadL(1, vsiFileSize, f)
                gdal.VSIFCloseL(f)
                arrayList.append(np.fromstring(strData, dtype=np.uint16, count=-1, sep='').reshape(self.dataset.RasterYSize, self.dataset.RasterXSize))
            # Get indices which have GCP
            index = np.where(arrayList[0] != NonValue)
            # Create GCPs and get projection
            gcpList = []
            for i in range(len(index[0])):
                igcp = gdal.GCP(arrayList[0][index[0][i]][index[1][i]] / 100.0 - 180.0,
                                arrayList[1][index[0][i]][index[1][i]] / 100.0 - 180.0,
                                0.0, index[1][i] + 0.5, index[0][i] + 0.5,
                                "", str(i+1))
                gcpList.append(igcp)
            self.dataset.SetGCPs(tuple(gcpList), projection)
        elif "NC_GLOBAL#GDAL_X_DATASET" in gdalMetadata and "NC_GLOBAL#GDAL_Y_DATASET" in gdalMetadata:
            """Is it better to use 'try', 'except' and check all keys?? """
            # set GEOLOCATION
            print "GEOLOCATION"
            geolocDict = {'LINE_OFFSET':gdalMetadata['NC_GLOBAL#GDAL_LINE_OFFSET'],
                          'LINE_STEP':gdalMetadata['NC_GLOBAL#GDAL_LINE_STEP'],
                          'PIXEL_OFFSET':gdalMetadata['NC_GLOBAL#GDAL_PIXEL_OFFSET'],
                          'PIXEL_STEP':gdalMetadata['NC_GLOBAL#GDAL_PIXEL_STEP'],
                          'SRS': projection,
                          'X_BAND':gdalMetadata['NC_GLOBAL#GDAL_X_BAND'],
                          'X_DATASET': gdalMetadata['NC_GLOBAL#GDAL_X_DATASET'],
                          'Y_BAND':gdalMetadata['NC_GLOBAL#GDAL_Y_BAND'],
                          'Y_DATASET': gdalMetadata['NC_GLOBAL#GDAL_Y_DATASET']}
            self.dataset.SetMetadata(geolocDict, 'GEOLOCATION')

        elif "NC_GLOBAL#GDAL_GeoTarnsform" in gdalMetadata:
            geoTransString = gdalMetadata['NC_GLOBAL#GDAL_GeoTarnsform'].split(" ")
            geoTransform = (float(geoTransString[0]), float(geoTransString[1]),
                            float(geoTransString[2]), float(geoTransString[3]),
                            float(geoTransString[4]), float(geoTransString[5]))
            self.dataset.SetGeoTransform(geoTransform)

        # for other netCDF files (which are not created by NANSAT)
        else:
            gdalSubDataset = gdal.Open(self.dataset.GetRasterBand(1).GetMetadata_Dict()['source'])
            projection = gdalSubDataset.GetProjection()
            if projection == '':
                projection = gdalSubDataset.GetGCPProjection()
            geoTransform = gdalSubDataset.GetGeoTransform()
            self.dataset.SetProjection(projection)
            self.dataset.SetGeoTransform(geoTransform)
        ##self.export("/Home/asumak/data/output/gcp166.vrt")