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

        # add bands with metadata and corresponding values to the empty VRT
        metaDict = []
        geoFileDict = {}
        xDatasetSource = ''
        yDatasetSource = ''
        projection = ''
        firstXSize = 0
        firstYSize = 0
        for i, fileName in enumerate(fileNames):
            subDataset = gdal.Open(fileName)
            # choose the first dataset whith grid
            if (firstXSize == 0 and firstYSize == 0 and
                    subDataset.RasterXSize > 1 and subDataset.RasterYSize > 1):
                firstXSize = subDataset.RasterXSize
                firstYSize = subDataset.RasterYSize
                firstSubDataset = subDataset
                
            # take bands whose sizes are same as the first band.
            if (subDataset.RasterXSize == firstXSize and
                        subDataset.RasterYSize == firstYSize):
                if 'GEOLOCATION_X_DATASET' in fileName:
                    xDatasetSource = fileName
                elif 'GEOLOCATION_Y_DATASET' in fileName:
                    yDatasetSource = fileName
                else:
                    for iBand in range(subDataset.RasterCount):
                        bandDict = subDataset.GetRasterBand(iBand+1).GetMetadata_Dict()
                        sourceBands = iBand + 1
                        # set wkv and bandname
                        wkv = bandDict.get('standard_name', '')
                        band_name = bandDict.get('NETCDF_VARNAME', '')
                        if len(band_name) > 0:
                            bandDict['band_name'] = band_name
                        # remove non-necessary metadata
                        for rmMetadata in rmMetadatas:
                            if rmMetadata in bandDict:
                                bandDict.pop(rmMetadata)
                        
                        metaDict.append(({'source': fileName, 'sourceBand': sourceBands, 'wkv': wkv,'parameters':bandDict}))
                    # get non-empty projection from one of the bands
                    srs = subDataset.GetProjection()
                    if len(srs) != 0:
                        projection = srs

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, firstSubDataset)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # if projetcion was not set automatically:
        #       get projection from GDAL_NANSAT_Projection
        if len(projection) == 0:
            projection = gdalMetadata.get('NC_GLOBAL#GDAL_NANSAT_Projection', '')
        # if no projection was found in netcdf, generate WGS84 by default
        if len(projection) == 0:
            sr = osr.SpatialReference()
            sr.ImportFromProj4("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs")
            projection = sr.ExportToWkt()
        # set projection
        self.dataset.SetProjection(self.repare_projection(projection))
        
        # ADD GCPs from metadata
        gcpCount = self.add_gcps_from_metadata(gdalMetadata)
        
        # Find proper bands and insert GEOLOCATION into dataset
        if len(xDatasetSource) > 0 and len(yDatasetSource) > 0:
            self.add_geolocation(Geolocation(xDatasetSource, yDatasetSource))
        elif gcpCount == 0:
            # if no GCPs found and not GEOLOCATION set: Set Nansat Geotransform if it is not set automatically
            geoTransform = self.dataset.GetGeoTransform()
            if len(geoTransform) == 0:
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
