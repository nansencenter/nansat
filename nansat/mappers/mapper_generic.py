# Name:         mapper_generic.py
# Purpose:      Generic Mapper for L3/L4 satellite or modeling data
# Authors:      Asuka Yamakava, Anton Korosov, Morten Wergeland Hansen,
#               Aleksander Vines
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
import os
from dateutil.parser import parse
import datetime

import numpy as np
from netCDF4 import Dataset

from nansat.nsr import NSR
from nansat.geolocation import Geolocation
from nansat.vrt import VRT
from nansat.tools import gdal, parse_time
from nansat.exceptions import WrongMapperError

# TODO: remove WrongMapperError
class Mapper(VRT):
    def __init__(self, inputFileName, gdalDataset, gdalMetadata, logLevel=30,
                 rmMetadatas=['NETCDF_VARNAME', '_Unsigned',
                              'ScaleRatio', 'ScaleOffset', 'dods_variable'],
                 **kwargs):
        # Remove 'NC_GLOBAL#' and 'GDAL_' and 'NANSAT_'
        # from keys in gdalDataset
        tmpGdalMetadata = {}
        geoMetadata = {}
        origin_is_nansat = False
        if not gdalMetadata:
            raise WrongMapperError
        for key in gdalMetadata.keys():
            newKey = key.replace('NC_GLOBAL#', '').replace('GDAL_', '')
            if 'NANSAT_' in newKey:
                geoMetadata[newKey.replace('NANSAT_', '')] = gdalMetadata[key]
                origin_is_nansat = True
            else:
                tmpGdalMetadata[newKey] = gdalMetadata[key]
        gdalMetadata = tmpGdalMetadata
        fileExt = os.path.splitext(inputFileName)[1]

        # Get file names from dataset or subdataset
        subDatasets = gdalDataset.GetSubDatasets()
        if len(subDatasets) == 0:
            filenames = [inputFileName]
        else:
            filenames = [f[0] for f in subDatasets]

        # add bands with metadata and corresponding values to the empty VRT
        metaDict = []
        xDatasetSource = ''
        yDatasetSource = ''
        firstXSize = 0
        firstYSize = 0
        for _, filename in enumerate(filenames):
            subDataset = gdal.Open(filename)
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
                if ('GEOLOCATION_X_DATASET' in filename or
                        'longitude' in filename):
                    xDatasetSource = filename
                elif ('GEOLOCATION_Y_DATASET' in filename or
                        'latitude' in filename):
                    yDatasetSource = filename
                else:
                    for iBand in range(subDataset.RasterCount):
                        subBand = subDataset.GetRasterBand(iBand+1)
                        bandMetadata = subBand.GetMetadata_Dict()
                        if 'PixelFunctionType' in bandMetadata:
                            bandMetadata.pop('PixelFunctionType')
                        sourceBands = iBand + 1
                        # sourceBands = i*subDataset.RasterCount + iBand + 1

                        # generate src metadata
                        src = {'SourceFilename': filename,
                               'SourceBand': sourceBands}
                        # set scale ratio and scale offset
                        scaleRatio = bandMetadata.get(
                            'ScaleRatio',
                            bandMetadata.get(
                                'scale',
                                bandMetadata.get('scale_factor', '')))
                        if len(scaleRatio) > 0:
                            src['ScaleRatio'] = scaleRatio
                        scaleOffset = bandMetadata.get(
                            'ScaleOffset',
                            bandMetadata.get(
                                'offset',
                                bandMetadata.get(
                                    'add_offset', '')))
                        if len(scaleOffset) > 0:
                            src['ScaleOffset'] = scaleOffset
                        # sate DataType
                        src['DataType'] = subBand.DataType

                        # generate dst metadata
                        # get all metadata from input band
                        dst = bandMetadata
                        # set wkv and bandname
                        dst['wkv'] = bandMetadata.get('standard_name', '')
                        # first, try the name metadata
                        if 'name' in bandMetadata:
                            bandName = bandMetadata['name']
                        else:
                            # if it doesn't exist get name from NETCDF_VARNAME
                            bandName = bandMetadata.get('NETCDF_VARNAME', '')
                            if len(bandName) == 0:
                                bandName = bandMetadata.get(
                                            'dods_variable', ''
                                            )

                            # remove digits added by gdal in
                            # exporting to netcdf...
                            if (len(bandName) > 0 and origin_is_nansat and
                                    fileExt == '.nc'):
                                if bandName[-1:].isdigit():
                                    bandName = bandName[:-1]
                                if bandName[-1:].isdigit():
                                    bandName = bandName[:-1]

                        # if still no bandname, create one
                        if len(bandName) == 0:
                            bandName = 'band_%03d' % iBand

                        dst['name'] = bandName

                        # remove non-necessary metadata from dst
                        for rmMetadata in rmMetadatas:
                            if rmMetadata in dst:
                                dst.pop(rmMetadata)

                        # append band with src and dst dictionaries
                        metaDict.append({'src': src, 'dst': dst})

        # create empty VRT dataset with geolocation only
        self._init_from_gdal_dataset(firstSubDataset, metadata=gdalMetadata)

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        self._create_complex_bands(filenames)

        if len(projection) == 0:
            # projection was not set automatically
            # get projection from GCPProjection
            projection = geoMetadata.get('GCPProjection', '')
        if len(projection) == 0:
            # no projection was found in dataset or metadata:
            # generate WGS84 by default
            projection = NSR().wkt
        # fix problem with MET.NO files where a, b given in m and XC/YC in km
        if ('UNIT["kilometre"' in projection and
            ',SPHEROID["Spheroid",6378273,7.331926543631893e-12]' in
                projection):
            projection = projection.replace(
                ',SPHEROID["Spheroid",6378273,7.331926543631893e-12]',
                '')
        # set projection
        self.dataset.SetProjection(self.repare_projection(projection))

        # check if GCPs were added from input dataset
        gcps = firstSubDataset.GetGCPs()
        gcpProjection = firstSubDataset.GetGCPProjection()

        # if no GCPs in input dataset: try to add GCPs from metadata
        if not gcps:
            gcps = self.add_gcps_from_metadata(geoMetadata)
        # if yet no GCPs: try to add GCPs from variables
        if not gcps:
            gcps = self.add_gcps_from_variables(inputFileName)

        if gcps:
            if len(gcpProjection) == 0:
                # get GCP projection and repare
                gcpProjection = self.repare_projection(geoMetadata. get('GCPProjection', ''))
            # add GCPs to dataset
            self.dataset.SetGCPs(gcps, gcpProjection)
            self.dataset.SetProjection('')
            self._remove_geotransform()

        # Find proper bands and insert GEOLOCATION ARRAY into dataset
        if len(xDatasetSource) > 0 and len(yDatasetSource) > 0:
            self._add_geolocation(Geolocation.from_filenames(xDatasetSource, yDatasetSource))

        elif not gcps:
            # if no GCPs found and not GEOLOCATION ARRAY set:
            #   Set Nansat Geotransform if it is not set automatically
            geoTransform = self.dataset.GetGeoTransform()
            if len(geoTransform) == 0:
                geoTransformStr = geoMetadata.get('GeoTransform',
                                                  '(0|1|0|0|0|0|1)')
                geoTransform = eval(geoTransformStr.replace('|', ','))
                self.dataset.SetGeoTransform(geoTransform)

        subMetadata = firstSubDataset.GetMetadata()


        ### GET START TIME from METADATA
        time_coverage_start = None
        if 'start_time' in gdalMetadata:
            time_coverage_start = parse_time(gdalMetadata['start_time'])
        elif 'start_date' in gdalMetadata:
            time_coverage_start = parse_time(gdalMetadata['start_date'])
        elif 'time_coverage_start' in gdalMetadata:
            time_coverage_start = parse_time(
                                        gdalMetadata['time_coverage_start'])

        ### GET END TIME from METADATA
        time_coverage_end = None
        if 'stop_time' in gdalMetadata:
            time_coverage_end = parse_time(gdalMetadata['stop_time'])
        elif 'stop_date' in gdalMetadata:
            time_coverage_end = parse_time(gdalMetadata['stop_date'])
        elif 'time_coverage_stop' in gdalMetadata:
            time_coverage_end = parse_time(
                                        gdalMetadata['time_coverage_stop'])
        elif 'end_time' in gdalMetadata:
            time_coverage_end = parse_time(gdalMetadata['end_time'])
        elif 'end_date' in gdalMetadata:
            time_coverage_end = parse_time(gdalMetadata['end_date'])
        elif 'time_coverage_end' in gdalMetadata:
            time_coverage_end = parse_time(
                                        gdalMetadata['time_coverage_end'])

        ### GET start time from time variable
        if (time_coverage_start is None and 'time#standard_name' in subMetadata and
                 subMetadata['time#standard_name'] == 'time' and 'time#units' in subMetadata):
            # get data from netcdf data
            ncFile = Dataset(inputFileName, 'r')
            time_var = ncFile.variables['time']
            t0 = time_var[0]
            if len(time_var) == 1:
                t1 = t0 + 1
            else:
                t1 = time_var[-1]

            time_units_start = parse(time_var.units, fuzzy=True, ignoretz=True)
            time_units_to_seconds = {'second' : 1.0,
                                     'hour' : 60 * 60.0,
                                     'day' : 24 * 60 * 60.0}
            for key in time_units_to_seconds:
                if key in time_var.units:
                    factor = time_units_to_seconds[key]
                    break

            time_coverage_start = time_units_start + datetime.timedelta(seconds=t0 * factor)
            time_coverage_end = time_units_start + datetime.timedelta(seconds=t1 * factor)

        ## finally set values of time_coverage start and end if available
        if time_coverage_start is not None:
            self.dataset.SetMetadataItem('time_coverage_start',
                                    time_coverage_start.isoformat())
        if time_coverage_end is not None:
            self.dataset.SetMetadataItem('time_coverage_end',
                                    time_coverage_end.isoformat())

        if 'sensor' not in gdalMetadata:
            self.dataset.SetMetadataItem('sensor', 'unknown')
        if 'satellite' not in gdalMetadata:
            self.dataset.SetMetadataItem('satellite', 'unknown')
        if 'source_type' not in gdalMetadata:
            self.dataset.SetMetadataItem('source_type', 'unknown')
        if 'platform' not in gdalMetadata:
            self.dataset.SetMetadataItem('platform', 'unknown')
        if 'instrument' not in gdalMetadata:
            self.dataset.SetMetadataItem('instrument', 'unknown')

        self.logger.info('Use generic mapper - OK!')

    def repare_projection(self, projection):
        '''Replace odd symbols in projection string '|' => ','; '&' => '"' '''
        return projection.replace("|", ",").replace("&", '"')

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
            gcpString = gcpString.strip().replace(' ', '')
            gcpValues = []
            # append all gcps from string
            for x in gcpString.split('|'):
                if len(x) > 0:
                    gcpValues.append(float(x))
            # gcpValues = [float(x) for x in gcpString.strip().split('|')]
            gcpAllValues.append(gcpValues)

        # create list of GDAL GCPs
        gcps = []
        for i in range(0, len(gcpAllValues[0])):
            gcps.append(gdal.GCP(gcpAllValues[2][i], gcpAllValues[3][i], 0,
                                 gcpAllValues[0][i], gcpAllValues[1][i]))

        return gcps

    def add_gcps_from_variables(self, filename):
        ''' Get GCPs from GCPPixel, GCPLine, GCPX, GCPY, GCPZ variables '''
        gcpVariables = ['GCPX', 'GCPY', 'GCPZ', 'GCPPixel', 'GCPLine', ]
        # open input netCDF file for reading GCPs
        try:
            ncFile = Dataset(filename, 'r')
        except (TypeError, IOError) as e:
            self.logger.info('%s' % e)
            return None

        # check if all GCP variables exist in the file
        if not all([var in ncFile.variables for var in gcpVariables]):
            return None

        # get data from GCP variables into array
        varData = [ncFile.variables[var][:] for var in gcpVariables]
        varData = np.array(varData)

        # close input file
        ncFile.close()

        # create list of GDAL.GCPs
        gcps = [gdal.GCP(float(x),
                         float(y),
                         float(z),
                         float(pixel),
                         float(line)) for x, y, z, pixel, line in varData.T]

        return gcps
