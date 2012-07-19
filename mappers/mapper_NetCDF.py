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

        rmMetadatas = ['NETCDF_VARNAME', '_FillValue', '_Unsigned']
        
        # Get file names from dataset or subdataset
        if gdalDataset.RasterCount==1:
            fileNames = [fileName]
        else:
            subDatasets = gdalDataset.GetSubDatasets()
            fileNames = [f[0] for f in subDatasets]

        # get raster size from the first band
        firstSubDataset = gdal.Open(fileNames[0])

        # add bands with metadata and corresponding values to the empty VRT
        metaDict = []
        geoFileDict = {}
        xDatasetSource = ''
        yDatasetSource = ''
        for i, fileName in enumerate(fileNames):
            subDataset = gdal.Open(fileName)
            # take bands whose sizes are same as the first band.
            if (subDataset.RasterXSize == firstSubDataset.RasterXSize and
                        subDataset.RasterYSize == firstSubDataset.RasterYSize):
                if 'GEOLOCATION_X_DATASET' in fileName:
                    xDatasetSource = fileName
                elif 'GEOLOCATION_Y_DATASET' in fileName:
                    yDatasetSource = fileName
                else:
                    for iBand in range(subDataset.RasterCount):
                        bandDict = subDataset.GetRasterBand(iBand+1).GetMetadata_Dict()
                        sourceBands = iBand + 1
                        wkv = bandDict.get('standard_name', '')
                        for rmMetadata in rmMetadatas:
                            if rmMetadata in bandDict:
                                bandDict.pop(rmMetadata)
                        
                        print 'bandDict: ', bandDict
                        metaDict.append(({'source': fileName, 'sourceBand': sourceBands, 'wkv': wkv,'parameters':bandDict}))

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, firstSubDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # set projection to dataset
        projection = gdalMetadata.get('NC_GLOBAL#GDAL_NANSAT_Projection', '')
        self.dataset.SetProjection(self.repare_projection(projection))
        
        # ADD GCPs from metadata
        gcpCount = self.add_gcps_from_metadata(gdalMetadata)
        
        # Find proper bands and insert GEOLOCATION into dataset
        if len(xDatasetSource) > 0 and len(yDatasetSource) > 0:
            self.add_geolocation(Geolocation(xDatasetSource, yDatasetSource))
        elif gcpCount == 0:
            # if no GCPs found and not GEOLOCATION set: Set Geotransform
            geoTransformStr = gdalMetadata.get('NC_GLOBAL#GDAL_NANSAT_GeoTransform', '(0|1|0|0|0|0|1)')
            geoTransform = eval(geoTransformStr.replace('|', ','))
            self.dataset.SetGeoTransform(geoTransform)

    def repare_projection(self, projection):
        '''Replace odd symbols in projection string '|' => ','; '&' => '"' '''
        return projection.replace("|",",").replace("&",'"')

    def add_gcps_from_metadata(self, gdalMetadata):
        '''Get GCPs from strings in metadata and insert in dataset'''
        gcpNames = ['GCPPixel', 'GCPLine', 'GCPX', 'GCPY']
        gcpAllValues = []
        # for all gcp coordinates
        for i, gcpName in enumerate(gcpNames):
            # scan throught metadata and find how many lines with each GCP
            gcpLineCount = 0
            for metaDataItem in gdalMetadata:
                if gcpName in metaDataItem:
                    gcpLineCount += 1
            # concat all lines
            gcpString = ''
            for n in range(0, gcpLineCount):
                gcpLineName = 'NC_GLOBAL#GDAL_NANSAT_%s_%03d' % (gcpName, n)
                gcpString += gdalMetadata[gcpLineName]
            # convert strings to floats
            gcpString = gcpString.strip().replace(' ','')
            gcpValues = []
            # append all gcps from string
            for x in gcpString.split('|'):
                if len(x) > 0:
                    gcpValues.append(float(x))
            #gcpValues = [float(x) for x in gcpString.strip().split('|')]
            gcpAllValues.append(gcpValues)
        
        # create list of GDAL GCPs
        gcps = []
        for i in range(0, len(gcpAllValues[0])):
            gcps.append(gdal.GCP(gcpAllValues[2][i], gcpAllValues[3][i], 0,
                                 gcpAllValues[0][i], gcpAllValues[1][i]))
        
        if len(gcps) > 0:
            # get GCP projection and repare
            projection = self.repare_projection(gdalMetadata.get('NC_GLOBAL#GDAL_NANSAT_GCPProjection', ''))
            # add GCPs to dataset
            self.dataset.SetGCPs(gcps, projection)
            self._remove_geotransform()
        
        return len(gcps)
