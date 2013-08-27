# Name:         mapper_generic.py
# Purpose:      Generic Mapper for L3/L4 satellite or modeling data
# Authors:      Asuka Yamakava, Anton Korosov, Morten Wergeland Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

#import os
from vrt import *
from nansat_tools import Node, latlongSRS
import numpy as np

class Mapper(VRT):
    def __init__(self, fileName, gdalDataset, gdalMetadata, logLevel=30,
                 rmMetadatas = ['NETCDF_VARNAME', '_Unsigned',
                                'ScaleRatio', 'ScaleOffset', 'dods_variable'],
                 **kwargs):
        # Remove 'NC_GLOBAL#' and 'GDAL_' and 'NANSAT_' from keys in gdalDataset
        tmpGdalMetadata = {}
        geoMetadata = {}
        origin_is_nansat = False
        for key in gdalMetadata.keys():
            newKey = key.replace('NC_GLOBAL#', '').replace('GDAL_', '')
            if 'NANSAT_' in newKey:
                geoMetadata[newKey.replace('NANSAT_', '')] = gdalMetadata[key]
                origin_is_nansat = True
            else:
                tmpGdalMetadata[newKey] = gdalMetadata[key]
        gdalMetadata = tmpGdalMetadata
        fileExt = os.path.splitext(fileName)[1]
        
        # Get file names from dataset or subdataset
        subDatasets = gdalDataset.GetSubDatasets()
        if len(subDatasets) == 0:
            fileNames = [fileName]
        else:
            fileNames = [f[0] for f in subDatasets]

        # add bands with metadata and corresponding values to the empty VRT
        metaDict = []
        geoFileDict = {}
        xDatasetSource = ''
        yDatasetSource = ''
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
                # get projection from the first subDataset
                projection = firstSubDataset.GetProjection()

            # take bands whose sizes are same as the first band.
            if (subDataset.RasterXSize == firstXSize and
                        subDataset.RasterYSize == firstYSize):
                if projection == '':
                    projection = subDataset.GetProjection()
                if 'GEOLOCATION_X_DATASET' in fileName or 'longitude' in fileName:
                    xDatasetSource = fileName
                elif 'GEOLOCATION_Y_DATASET' in fileName or 'latitude' in fileName:
                    yDatasetSource = fileName
                else:
                    for iBand in range(subDataset.RasterCount):
                        subBand = subDataset.GetRasterBand(iBand+1)
                        bandMetadata = subBand.GetMetadata_Dict()
                        if 'PixelFunctionType' in bandMetadata:
                            bandMetadata.pop('PixelFunctionType')
                        sourceBands = iBand + 1
                        #sourceBands = i*subDataset.RasterCount + iBand + 1

                        # generate src metadata
                        src = {'SourceFilename': fileName, 'SourceBand': sourceBands}
                        # set scale ratio and scale offset
                        scaleRatio = bandMetadata.get('ScaleRatio',
                                     bandMetadata.get('scale',
                                     bandMetadata.get('scale_factor', '')))
                        if len(scaleRatio) > 0:
                            src['ScaleRatio'] = scaleRatio
                        scaleOffset = bandMetadata.get('ScaleOffset',
                                      bandMetadata.get('offset',
                                      bandMetadata.get('add_offset', '')))
                        if len(scaleOffset) > 0:
                            src['ScaleOffset'] = scaleOffset
                        # sate DataType
                        src['DataType'] = subBand.DataType

                        # generate dst metadata
                        # get all metadata from input band
                        dst = bandMetadata
                        # set wkv and bandname
                        dst['wkv'] = bandMetadata.get('standard_name', '')
                        bandName = bandMetadata.get('NETCDF_VARNAME', '') # could we also use bandMetadata.get('name')?
                        if len(bandName) == 0:
                            bandName = bandMetadata.get('dods_variable', '')
                        if len(bandName) > 0:
                            if origin_is_nansat and fileExt == '.nc':
                                # remove digits added by gdal in exporting to
                                # netcdf...
                                if bandName[-1:].isdigit():
                                    bandName=bandName[:-1]
                                if bandName[-1:].isdigit():
                                    bandName=bandName[:-1]
                                dst['name'] = bandName
                            else:
                                dst['name'] = bandName

                        # remove non-necessary metadata from dst
                        for rmMetadata in rmMetadatas:
                            if rmMetadata in dst:
                                dst.pop(rmMetadata)

                        # append band with src and dst dictionaries
                        metaDict.append({'src': src, 'dst': dst})

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, firstSubDataset, srcMetadata=gdalMetadata)

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        # Create complex data bands from 'xxx_real' and 'xxx_imag' bands
        # using pixelfunctions
        rmBands = []
        for iBand in range(self.dataset.RasterCount):
            iBandName = self.dataset.GetRasterBand(iBand+1).GetMetadataItem('name')
            # find real data band
            if iBandName.find("_real") != -1:
                realBand = iBand
                realDtype = self.dataset.GetRasterBand(realBand+1).GetMetadataItem('DataType')
                bandName = iBandName.replace(iBandName.split('_')[-1], '')[0:-1]
                for jBand in range(self.dataset.RasterCount):
                    jBandName = self.dataset.GetRasterBand(jBand+1).GetMetadataItem('name')
                    # find an imaginary data band corresponding to the real data band
                    # and create complex data band from the bands
                    if jBandName.find(bandName+'_imag') != -1:
                        imagBand = jBand
                        imagDtype = self.dataset.GetRasterBand(imagBand+1).GetMetadataItem('DataType')
                        dst = self.dataset.GetRasterBand(imagBand+1).GetMetadata()
                        dst['name'] = bandName
                        dst['PixelFunctionType'] ='ComplexData'
                        dst['dataType'] = 10
                        src = [{'SourceFilename': fileNames[realBand],
                                'SourceBand':  1,
                                'DataType': realDtype},
                               {'SourceFilename': fileNames[imagBand],
                                 'SourceBand': 1,
                                 'DataType': imagDtype}]
                        self._create_band(src, dst)
                        self.dataset.FlushCache()
                        rmBands.append(realBand+1)
                        rmBands.append(imagBand+1)

        # Delete real and imaginary bands
        if len(rmBands) != 0:
            self.delete_bands(rmBands)

        if len(projection) == 0:
            # projection was not set automatically
            # get projection from GCPProjection
            projection = geoMetadata.get('GCPProjection', '')
        if len(projection) == 0:
            # no projection was found in dataset or metadata:
            # generate WGS84 by default
            projection = latlongSRS.ExportToWkt()
        # set projection
        self.dataset.SetProjection(self.repare_projection(projection))

        # check if GCPs were added from input dataset
        gcpCount = firstSubDataset.GetGCPCount()
        if gcpCount == 0:
            # if no GCPs in input dataset: try to add GCPs from metadata
            gcpCount = self.add_gcps_from_metadata(geoMetadata)

        # Find proper bands and insert GEOLOCATION ARRAY into dataset
        if len(xDatasetSource) > 0 and len(yDatasetSource) > 0:
            self.add_geolocationArray(GeolocationArray(xDatasetSource, yDatasetSource))

        elif gcpCount == 0:
            # if no GCPs found and not GEOLOCATION ARRAY set:
            #   Set Nansat Geotransform if it is not set automatically
            geoTransform = self.dataset.GetGeoTransform()
            if len(geoTransform) == 0:
                geoTransformStr = geoMetadata.get('GeoTransform', '(0|1|0|0|0|0|1)')
                geoTransform = eval(geoTransformStr.replace('|', ','))
                self.dataset.SetGeoTransform(geoTransform)

        if 'start_date' in gdalMetadata:
            self._set_time(parse(gdalMetadata['start_date']))

    def repare_projection(self, projection):
        '''Replace odd symbols in projection string '|' => ','; '&' => '"' '''
        return projection.replace("|",",").replace("&",'"')

    def add_gcps_from_metadata(self, geoMetadata):
        '''Get GCPs from strings in metadata and insert in dataset'''
        gcpNames = ['GCPPixel', 'GCPLine', 'GCPX', 'GCPY']
        gcpAllValues = []
        # for all gcp coordinates
        for i, gcpName in enumerate(gcpNames):
            # scan throught metadata and find how many lines with each GCP
            gcpLineCount = 0
            for metaDataItem in geoMetadata:
                if gcpName in metaDataItem:
                    gcpLineCount += 1
            # concat all lines
            gcpString = ''
            for n in range(0, gcpLineCount):
                gcpLineName = '%s_%03d' % (gcpName, n)
                gcpString += geoMetadata[gcpLineName]
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
            projection = self.repare_projection(geoMetadata.get('GCPProjection', ''))
            # add GCPs to dataset
            self.dataset.SetGCPs(gcps, projection)
            self._remove_geotransform()

        return len(gcps)
