# Name:    nansat.py
# Name:  nansat.py
# Purpose: Container of Nansat class
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
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
from __future__ import absolute_import
import os
import glob
import sys
import tempfile
import datetime
import pkgutil
import warnings

import numpy as np
if 'nanmedian' in np.__all__:
    from numpy import nanmedian
else:
    from scipy.stats import nanmedian

from numpy.lib.recfunctions import append_fields
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from netCDF4 import Dataset

from nansat.nsr import NSR
from nansat.domain import Domain
from nansat.figure import Figure
from nansat.vrt import VRT
from nansat.tools import add_logger, gdal
from nansat.tools import OptionError, WrongMapperError, NansatReadError, GDALError
from nansat.tools import parse_time, test_openable
from nansat.node import Node
from nansat.pointbrowser import PointBrowser
import collections
if hasattr(collections, 'OrderedDict'):
    from collections import OrderedDict
else:
    from ordereddict import OrderedDict

# container for all mappers
nansatMappers = None


class Nansat(Domain):
    '''Container for geospatial data, performs all high-level operations

    n = Nansat(fileName) opens the file with satellite or model data for
    reading, adds scientific metadata to bands, and prepares the data for
    further handling.

    The instance of Nansat class (the object <n>) contains information
    about geographical reference of the data (e.g raster size, pixel
    resolution, type of projection, etc) and about bands with values of
    geophysical variables (e.g. water leaving radiance, normalized radar
    cross section, chlrophyll concentraion, etc). The object <n> has methods
    for high-level operations with data. E.g.:
    * reading data from file (Nansat.__getitem__);
    * visualization (Nansat.write_figure);
    * changing geographical reference (Nansat.reproject);
    * exporting (Nansat.export)
    * and much more...

    Nansat inherits from Domain (container of geo-reference information)
    Nansat uses instance of VRT (wraper around GDAL VRT-files)
    Nansat uses instance of Figure (collection of methods for visualization)
    '''

    def __init__(self, fileName='', mapperName='', domain=None,
                 array=None, parameters=None, logLevel=30, **kwargs):
        '''Create Nansat object

        if <fileName> is given:
            Open GDAL dataset,
            Read metadata,
            Generate GDAL VRT file with mapping of variables in memory
            Create logger
            Create Nansat object for perfroming high-level operations
        if <domain> and <array> are given:
            Create VRT object from data in <array>
            Add geolocation from <domain>

        Parameters
        -----------
        fileName : string
            location of the file
        mapperName : string, optional
            name of the mapper from nansat/mappers dir. E.g.
            'ASAR', 'hirlam', 'merisL1', 'merisL2', etc.
        domain : Domain object
            Geo-reference of a new raster
        array : numpy array
            Firts band of a new raster
        parameters : dictionary
            Metadata for the 1st band of a new raster,e.g. name, wkv, units,...
        logLevel : int, optional, default: logging.DEBUG (30)
            Level of logging. See: http://docs.python.org/howto/logging.html
        kwargs : additional arguments for mappers

        Creates
        --------
        self.mapper : str
            name of the used mapper
        self.fileName : file name
            set file name given by the argument
        self.vrt : VRT object
            Wrapper around VRT file and GDAL dataset with satellite raster data
        self.logger : logging.Logger
            logger for output debugging info
        self.name : string
            name of object (for writing KML)

        Examples
        --------
        n = Nansat(filename)
        # opens file for reading. Opening is lazy - no data is read at this
        # point, only metadata that describes the dataset and bands

        n = Nansat(domain=d)
        # create an empty Nansat object. <d> is the Domain object which
        # describes the grid (projection, resolution and extent)

        n = Nansat(domain=d, array=a, parameters=p)
        # create a Nansat object in memory with one band from input array <a>.
        # <p> is a dictionary with metadata for the band

        a = n[1]
        # fetch data from Nansat object from the first band

        a = n['band_name']
        # fetch data from the band which has name 'band_name'

        '''
        # check the arguments
        if fileName == '' and domain is None:
            raise OptionError('Either fileName or domain is required.')

        # create logger
        self.logger = add_logger('Nansat', logLevel)

        # empty dict of VRTs with added bands
        self.addedBands = {}

        # set input file name
        self.fileName = fileName
        # name, for compatibility with some Domain methods
        self.name = os.path.basename(fileName)
        self.path = os.path.dirname(fileName)

        # create self.vrt from a file using mapper or...
        if fileName != '':
            # Make original VRT object with mapping of variables
            self.vrt = self._get_mapper(mapperName, **kwargs)
        # ...create using array, domain, and parameters
        else:
            # Set current VRT object
            self.vrt = VRT(gdalDataset=domain.vrt.dataset)
            self.domain = domain
            self.mapper = ''
            if array is not None:
                # add a band from array
                self.add_band(array=array, parameters=parameters)

        self.logger.debug('Object created from %s ' % self.fileName)

    def __getitem__(self, bandID):
        ''' Returns the band as a NumPy array, by overloading []

        Parameters
        -----------
        bandID : int or str
            If int, array from band with number <bandID> is returned
            If string, array from band with metadata 'name' equal to
            <bandID> is returned

        Returns
        --------
        self.get_GDALRasterBand(bandID).ReadAsArray() : NumPy array

        '''
        # get band
        band = self.get_GDALRasterBand(bandID)
        # get expression from metadata
        expression = band.GetMetadata().get('expression', '')
        # get data
        bandData = band.ReadAsArray()
        if bandData is None:
            raise GDALError('Cannot read array from band %s' % str(bandID))

        # execute expression if any
        if expression != '':
            bandData = eval(expression)

        # Set invalid and missing data to np.nan (for floats only)
        if ('_FillValue' in band.GetMetadata() and
             bandData.dtype.char in np.typecodes['AllFloat']):
            fillValue = float(band.GetMetadata()['_FillValue'])
            bandData[bandData == fillValue] = np.nan
            # quick hack to avoid problem with wrong _FillValue - see issue
            # #123
            if fillValue == 9.96921e+36:
                altFillValue = -10000.
                bandData[bandData == altFillValue] = np.nan

        # replace infs with np.NAN
        if np.size(np.where(np.isinf(bandData))) > 0:
            bandData[np.isinf(bandData)] = np.nan

        # erase out-of-swath pixels with np.Nan (if not integer)
        if (self.has_band('swathmask') and bandData.dtype.char in
                                            np.typecodes['AllFloat']):
            swathmask = self.get_GDALRasterBand('swathmask').ReadAsArray()
            bandData[swathmask == 0] = np.nan

        return bandData

    def __repr__(self):
        '''Creates string with basic info about the Nansat object'''

        outString = '-' * 40 + '\n'
        outString += self.fileName + '\n'
        outString += '-' * 40 + '\n'
        outString += 'Mapper: ' + self.mapper + '\n'
        outString += '-' * 40 + '\n'
        outString += self.list_bands(False)
        outString += '-' * 40 + '\n'
        outString += Domain.__repr__(self)
        return outString

    def add_band(self, array, parameters=None, nomem=False):
        '''Add band from the array to self.vrt

        Create VRT object which contains VRT and RAW binary file and append it
        to self.vrt.bandVRTs

        Parameters
        -----------
        array : ndarray
            band data
        parameters : dictionary
            band metadata: wkv, name, etc. (or for several bands)
        nomem : boolean, saves the vrt to a tempfile if nomem is True

        Modifies
        ---------
        Creates VRT object with VRT-file and RAW-file
        Adds band to the self.vrt

        Examples
        --------
        n.add_band(a, p)
        # add new band from numpy array <a> with metadata <p> in memory
        # Shape of a should be equal to the shape of <n>

        n.add_band(a, p, nomem=True)
        # add new band from an array <a> with metadata <p> but keep it
        # temporarli on disk intead of memory
        '''
        self.add_bands([array], [parameters], nomem)

    def add_bands(self, arrays, parameters=None, nomem=False):
        '''Add band from the array to self.vrt

        Create VRT object which contains VRT and RAW binary file and append it
        to self.vrt.bandVRTs

        Parameters
        -----------
        array : ndarray or list
            band data (or data for several bands)
        parameters : dictionary or list
            band metadata: wkv, name, etc. (or for several bands)
        nomem : boolean, saves the vrt to a tempfile if nomem is True

        Modifies
        ---------
        Creates VRT object with VRT-file and RAW-file
        Adds band to the self.vrt

        Examples
        --------
        n.add_bands([a1, a2], [p1, p2])
        # add two new bands from numpy arrays <a1> and <a2> with metadata in
        # <p1> and <p2>

        '''
        # replace empty parameters with list of None
        if parameters is None:
            parameters = [None] * len(arrays)

        # create VRTs from arrays
        bandVRTs = [VRT(array=array, nomem=nomem) for array in arrays]

        self.vrt = self.vrt.get_super_vrt()

        # add the array band into self.vrt and get bandName
        for bi, bandVRT in enumerate(bandVRTs):
            params = parameters[bi]
            if params is None:
                params = {}
            bandName = self.vrt._create_band(
                {'SourceFilename': bandVRT.fileName,
                 'SourceBand': 1},
                params)
            self.vrt.bandVRTs[bandName] = bandVRT

        self.vrt.dataset.FlushCache()  # required after adding bands

    def bands(self):
        ''' Make a dictionary with all metadata from all bands

        Returns
        --------
        b : dictionary
            key = N, value = dict with all band metadata

        '''
        b = {}
        for iBand in range(self.vrt.dataset.RasterCount):
            b[iBand + 1] = self.get_metadata(bandID=iBand + 1)

        return b

    def has_band(self, band):
        '''Check if self has band with name <band>
        Parameters
        ----------
            band : str
                name or standard_name of the band to check

        Returns
        -------
            True/False if band exists or not

        '''
        bandExists = False
        for b in self.bands():
            if (self.bands()[b]['name'] == band or
                 ('standard_name' in self.bands()[b] and
                  self.bands()[b]['standard_name'] == band)):
                bandExists = True

        return bandExists

    def export(self, fileName, bands=None, rmMetadata=[], addGeolocArray=True,
               addGCPs=True, driver='netCDF', bottomup=False, options=None):
        '''Export Nansat object into netCDF or GTiff file

        Parameters
        -----------
        fileName : str
            output file name
        bands: list (default=None)
            Specify band numbers to export.
            If None, all bands are exported.
        rmMetadata : list
            metadata names for removal before export.
            e.g. ['name', 'colormap', 'source', 'sourceBands']
        addGeolocArray : bool
            add geolocation array datasets to exported file?
        addGCPs : bool
            add GCPs to exported file?
        driver : str
            Name of GDAL driver (format)
        bottomup : bool
            False: Default. Write swath-projected data with rows and columns
                   organized as in the original product.
            True:  Use the default behaviour of GDAL, which is to flip the rows
        options : str or list
            GDAL export options in format of: 'OPT=VAL', or
            ['OPT1=VAL1', 'OP2='VAL2']
            See also http://www.gdal.org/frmt_netcdf.html


        Modifies
        ---------
        Create a netCDF file

        !! NB
        ------
        If number of bands is more than one,
        serial numbers are added at the end of each band name.

        It is possible to fix it by changing
        line.4605 in GDAL/frmts/netcdf/netcdfdataset.cpp :
        'if( nBands > 1 ) sprintf(szBandName,"%s%d",tmpMetadata,iBand);'
        --> 'if( nBands > 1 ) sprintf(szBandName,"%s",tmpMetadata);'

        CreateCopy fails in case the band name has special characters,
        e.g. the slash in 'HH/VV'.

        Metadata strings with special characters are escaped with XML/HTML
        encoding.

        Examples
        --------
        n.export(netcdfile)
        # export all the bands into a netDCF 3 file

        n.export(driver='GTiff')
        # export all bands into a GeoTiff

        '''
        # temporary VRT for exporting
        exportVRT = self.vrt.copy()
        exportVRT.real = []
        exportVRT.imag = []

        # delete unnecessary bands
        rmBands = []
        selfBands = self.bands()
        if bands is not None:
            for selfBand in selfBands:
                # if band number or band name is not listed: mark for removal
                if (selfBand not in bands and
                      selfBands[selfBand]['name'] not in bands):
                    rmBands.append(selfBand)
            # delete bands from VRT
            exportVRT.delete_bands(rmBands)

        # Find complex data band
        complexBands = []
        node0 = Node.create(exportVRT.read_xml())
        for iBand in node0.nodeList('VRTRasterBand'):
            dataType = iBand.getAttribute('dataType')
            if dataType[0] == 'C':
                complexBands.append(int(iBand.getAttribute('band')))

        # if data includes complex data,
        # create two bands from real and imaginary data arrays
        if len(complexBands) != 0:
            for i in complexBands:
                bandMetadataR = self.get_metadata(bandID=i)
                bandMetadataR.pop('dataType')
                if 'PixelFunctionType' in bandMetadataR:
                    bandMetadataR.pop('PixelFunctionType')
                # Copy metadata and modify 'name' for real and imag bands
                bandMetadataI = bandMetadataR.copy()
                bandMetadataR['name'] = bandMetadataR.pop('name') + '_real'
                bandMetadataI['name'] = bandMetadataI.pop('name') + '_imag'
                # Create bands from the real and imaginary numbers
                exportVRT.real.append(VRT(array=self[i].real))
                exportVRT.imag.append(VRT(array=self[i].imag))

                metaDict = [{'src': {
                             'SourceFilename': exportVRT.real[-1].fileName,
                             'SourceBand':  1},
                             'dst': bandMetadataR},
                            {'src': {
                             'SourceFilename': exportVRT.imag[-1].fileName,
                             'SourceBand':  1},
                             'dst': bandMetadataI}]
                exportVRT._create_bands(metaDict)
            # delete the complex bands
            exportVRT.delete_bands(complexBands)

        # add bands with geolocation arrays to the VRT
        if addGeolocArray and len(exportVRT.geolocationArray.d) > 0:
            exportVRT._create_band(
                {'SourceFilename': self.vrt.geolocationArray.d['X_DATASET'],
                 'SourceBand': int(self.vrt.geolocationArray.d['X_BAND'])},
                {'wkv': 'longitude',
                 'name': 'GEOLOCATION_X_DATASET'})
            exportVRT._create_band(
                {'SourceFilename': self.vrt.geolocationArray.d['Y_DATASET'],
                 'SourceBand': int(self.vrt.geolocationArray.d['Y_BAND'])},
                {'wkv': 'latitude',
                 'name': 'GEOLOCATION_Y_DATASET'})

        # manage metadata for each band
        for iBand in range(exportVRT.dataset.RasterCount):
            band = exportVRT.dataset.GetRasterBand(iBand + 1)
            bandMetadata = band.GetMetadata()
            # set NETCDF_VARNAME
            try:
                bandMetadata['NETCDF_VARNAME'] = bandMetadata['name']
            except:
                self.logger.warning('Unable to set NETCDF_VARNAME for band %d'
                                    % (iBand + 1))
            # remove unwanted metadata from bands
            for rmMeta in rmMetadata:
                try:
                    bandMetadata.pop(rmMeta)
                except:
                    self.logger.info('Unable to remove metadata'
                                     '%s from band %d' % (rmMeta, iBand + 1))
            band.SetMetadata(bandMetadata)

        # remove unwanted global metadata
        globMetadata = exportVRT.dataset.GetMetadata()
        for rmMeta in rmMetadata:
            try:
                globMetadata.pop(rmMeta)
            except:
                self.logger.info('Global metadata %s not found' % rmMeta)

        # Apply escaping to metadata strings to preserve special characters (in
        # XML/HTML format)
        globMetadata_escaped = {}
        for key, val in globMetadata.iteritems():
            # Keys not escaped - this may be changed if needed...
            globMetadata_escaped[key] = gdal.EscapeString(val, gdal.CPLES_XML)
        exportVRT.dataset.SetMetadata(globMetadata_escaped)

        # if output filename is same as input one...
        if self.fileName == fileName:
            numOfBands = self.vrt.dataset.RasterCount
            # create VRT from each band and add it
            for iBand in range(numOfBands):
                vrt = VRT(array=self[iBand + 1])
                self.add_band(vrt=vrt)
                metadata = self.get_metadata(bandID=iBand + 1)
                self.set_metadata(key=metadata,
                                  bandID=numOfBands + iBand + 1)

            # remove source bands
            self.vrt.delete_bands(range(1, numOfBands))

        # get CreateCopy() options
        if options is None:
            options = []
        if type(options) == str:
            options = [options]

        # set bottomup option
        if bottomup:
            options += ['WRITE_BOTTOMUP=NO']
        else:
            options += ['WRITE_BOTTOMUP=YES']

        # if GCPs should be added
        gcps = exportVRT.dataset.GetGCPs()
        srs = exportVRT.get_projection()
        addGCPs = addGCPs and driver == 'netCDF' and len(gcps) > 0
        if addGCPs:
            #  remove GeoTransform
            exportVRT._remove_geotransform()
            exportVRT.dataset.SetMetadataItem(
                'NANSAT_GCPProjection', srs.replace(',', '|').replace('"', '&'))
        elif driver == 'GTiff':
            #  remove GeoTransform
            exportVRT._remove_geotransform()
        else:
            # add projection metadata
            exportVRT.dataset.SetMetadataItem(
                'NANSAT_Projection', srs.replace(',', '|').replace('"', '&'))

            # add GeoTransform metadata
            geoTransformStr = str(
                    self.vrt.dataset.GetGeoTransform()).replace(',', '|')
            exportVRT.dataset.SetMetadataItem(
                    'NANSAT_GeoTransform', geoTransformStr)

        # Create an output file using GDAL
        self.logger.debug('Exporting to %s using %s and %s...' % (fileName,
                                                                  driver,
                                                                  options))

        dataset = gdal.GetDriverByName(driver).CreateCopy(fileName,
                                                          exportVRT.dataset,
                                                          options=options)
        dataset = None
        # add GCPs into netCDF file as separate float variables
        if addGCPs:
            self._add_gcps(fileName, gcps, bottomup)

        self.logger.debug('Export - OK!')

    def _add_gcps(self, fileName, gcps, bottomup):
        ''' Add 4 variables with gcps to the generated netCDF file '''
        gcpVariables = ['GCPX', 'GCPY', 'GCPZ', 'GCPPixel', 'GCPLine', ]

        # check if file exists
        if not os.path.exists(fileName):
            self.logger.warning('Cannot add GCPs! %s doesn''t exist!' % fileName)
            return 1

        # open output file for adding GCPs
        ncFile = Dataset(fileName, 'a')

        # get GCP values into single array from GCPs
        gcpValues = np.zeros((5, len(gcps)))
        for i, gcp in enumerate(gcps):
            gcpValues[0, i] = gcp.GCPX
            gcpValues[1, i] = gcp.GCPY
            gcpValues[2, i] = gcp.GCPZ
            gcpValues[3, i] = gcp.GCPPixel
            gcpValues[4, i] = gcp.GCPLine

        # make gcps dimentions
        ncFile.createDimension('gcps', len(gcps))
        # make gcps variables and add data
        for i, var in enumerate(gcpVariables):
            var = ncFile.createVariable(var, 'f4', ('gcps',))
            var[:] = gcpValues[i]

        # write data, close file
        ncFile.close()

    def export2thredds(self, fileName, bands, metadata=None,
                       maskName=None, rmMetadata=[],
                       time=None, createdTime=None):
        ''' Export data into a netCDF formatted for THREDDS server

        Parameters
        -----------
        fileName : str
            output file name
        bands : dict
            {'band_name': {'type'     : '>i1',
                           'scale'    : 0.1,
                           'offset'   : 1000,
                           'metaKey1' : 'meta value 1',
                           'metaKey2' : 'meta value 2'}}
            dictionary sets parameters for band creation
            'type' - string representation of data type in the output band
            'scale' - sets scale_factor and applies scaling
            'offset' - sets 'scale_offset and applies offsetting
            other entries (e.g. 'units': 'K') set other metadata
        metadata : dict
            Glbal metadata to add
        maskName: string;
            if data include a mask band: give the mask name.
            Non-masked value is 64.
            if None: no mask is added
        rmMetadata : list
            unwanted metadata names which will be removed
        time : list with datetime objects
            aqcuisition time of original data. That value will be in time dim
        createdTime : datetime
            date of creation. Will be in metadata 'created'

        !! NB
        ------
        Nansat object (self) has to be projected (with valid GeoTransform and
        valid Spatial reference information) but not wth GCPs

        Examples
        --------
        # create THREDDS formatted netcdf file with all bands and time variable
        n.export2thredds(filename)

        # export only the first band and add global metadata
        n.export2thredds(filename, ['L_469'], {'description': 'example'})

        # export several bands and modify type, scale and offset
        bands = {'L_645' : {'type': '>i2', 'scale': 0.1, 'offset': 0},
                 'L_555' : {'type': '>i2', 'scale': 0.1, 'offset': 0}}
        n.export2thredds(filename, bands)

        '''
        # raise error if self is not projected (has GCPs)
        if len(self.vrt.dataset.GetGCPs()) > 0:
            raise OptionError('Cannot export dataset with GCPS for THREDDS!')

        # replace bands as list with bands as dict
        if type(bands) is list:
            bands = dict.fromkeys(bands, {})

        # Create temporary empty Nansat object with self domain
        data = Nansat(domain=self)

        # get mask (if exist)
        if maskName is not None:
            mask = self[maskName]

        # add required bands to data
        dstBands = {}
        srcBands = [self.bands()[b]['name'] for b in self.bands()]
        for iband in bands:
            # skip non exiting bands
            if iband not in srcBands:
                self.logger.error('%s is not found' % str(iband))
                continue

            array = self[iband]

            # catch None band error
            if array is None:
                raise GDALError('%s is None' % str(iband))

            # set type, scale and offset from input data or by default
            dstBands[iband] = {}
            dstBands[iband]['type'] = bands[iband].get('type',
                                             array.dtype.str.replace('u', 'i'))
            dstBands[iband]['scale'] = float(bands[iband].get('scale', 1.0))
            dstBands[iband]['offset'] = float(bands[iband].get('offset', 0.0))
            if '_FillValue' in bands[iband]:
                dstBands[iband]['_FillValue'] = np.array(
                                            [bands[iband]['_FillValue']],
                                            dtype=dstBands[iband]['type'])[0]

            # mask values with np.nan
            if maskName is not None and iband != maskName:
                array[mask != 64] = np.nan

            # add array to a temporary Nansat object
            bandMetadata = self.get_metadata(bandID=iband)
            data.add_band(array=array, parameters=bandMetadata)
        self.logger.debug('Bands for export: %s' % str(dstBands))

        # get corners of reprojected data
        minLat, maxLat, minLon, maxLon = data.get_min_max_lat_lon()

        # common global attributes:
        if createdTime is None:
            createdTime = (datetime.datetime.utcnow().
                           strftime('%Y-%m-%d %H:%M:%S UTC'))

        globMetadata = {'institution': 'NERSC',
                        'source': 'satellite remote sensing',
                        'creation_date': createdTime,
                        'northernmost_latitude': np.float(maxLat),
                        'southernmost_latitude': np.float(minLat),
                        'westernmost_longitude': np.float(minLon),
                        'easternmost_longitude': np.float(maxLon),
                        'history': ' '}

        # join or replace default by custom global metadata
        if metadata is not None:
            for metaKey in metadata:
                globMetadata[metaKey] = metadata[metaKey]

        # export temporary Nansat object to a temporary netCDF
        fid, tmpName = tempfile.mkstemp(suffix='.nc')
        data.export(tmpName)

        # open files for input and output
        ncI = Dataset(tmpName, 'r')
        ncO = Dataset(fileName, 'w')

        # collect info on dimention names
        dimNames = []
        gridMappingName = None
        for ncIVarName in ncI.variables:
            ncIVar = ncI.variables[ncIVarName]
            dimNames += list(ncIVar.dimensions)
            # get grid_mapping_name
            if hasattr(ncIVar, 'grid_mapping_name'):
                gridMappingName = ncIVar.grid_mapping_name
                gridMappingVarName = ncIVarName
        dimNames = list(set(dimNames))

        # collect info on dimention shapes
        dimShapes = {}
        for dimName in dimNames:
            dimVar = ncI.variables[dimName]
            dimShapes[dimName] = dimVar.shape[0]

        # create dimensions
        for dimName in dimNames:
            ncO.createDimension(dimName, dimShapes[dimName])

        # add time dimention
        ncO.createDimension('time', 1)
        ncOVar = ncO.createVariable('time', '>f8',  ('time', ))
        ncOVar.calendar = 'standard'
        ncOVar.long_name = 'time'
        ncOVar.standard_name = 'time'
        ncOVar.units = 'days since 1900-1-1 0:0:0 +0'
        ncOVar.axis = 'T'

        # get time from Nansat object or from input datetime
        if time is None:
            time = self.time_coverage_start

        # create value of time variable
        td = time - datetime.datetime(1900, 1, 1)
        days = td.days + (float(td.seconds) / 60.0 / 60.0 / 24.0)
        # add date
        ncOVar[:] = days

        # recreate file
        for ncIVarName in ncI.variables:
            ncIVar = ncI.variables[ncIVarName]
            if 'name' in ncIVar.ncattrs():
                ncIVar_name = ncIVar.getncattr('name')
            else:
                ncIVar_name = None

            self.logger.debug('Creating variable: %s' % ncIVarName)
            if ncIVarName in ['x', 'y', 'lon', 'lat']:
                # create simple x/y variables
                ncOVar = ncO.createVariable(ncIVarName, '>f4',
                                            ncIVar.dimensions)
            elif ncIVarName == gridMappingVarName:
                # create projection var
                ncOVar = ncO.createVariable(gridMappingName, ncIVar.dtype.str,
                                            ncIVar.dimensions)
            elif ncIVar_name in dstBands:
                # dont add time-axis to lon/lat grids
                if ncIVar.name in ['lon', 'lat']:
                    dimensions = ncIVar.dimensions
                else:
                    dimensions = ('time', ) + ncIVar.dimensions

                ncOVar = ncO.createVariable(ncIVar.name,
                                            dstBands[ncIVar_name]['type'],
                                            dimensions)

            # copy array from input data
            data = ncIVar[:]

            # copy rounded data from x/y
            if ncIVarName in ['x', 'y']:
                ncOVar[:] = np.floor(data).astype('>f4')
                # add axis=X or axis=Y
                ncOVar.axis = {'x': 'X', 'y': 'Y'}[ncIVarName]

            # copy data from lon/lat
            if ncIVarName in ['lon', 'lat']:
                ncOVar[:] = data.astype('>f4')

            # copy data from variables in the list
            if (len(ncIVar.dimensions) > 0 and ncIVar_name):
                # add offset and scale attributes
                scale = dstBands[ncIVar_name]['scale']
                offset = dstBands[ncIVar_name]['offset']
                if not (offset == 0.0 and scale == 1.0):
                    ncOVar.setncattr('add_offset', offset)
                    ncOVar.setncattr('scale_factor', scale)
                    data = (data - offset) / scale
                # replace non-value by '_FillValue'
                if (ncIVar.name in dstBands):
                    if '_FillValue' in dstBands[ncIVar_name].keys():
                        data[np.isnan(data)] = dstBands[
                                                ncIVar_name]['_FillValue']
                        ncOVar._attributes['_FillValue'] = dstBands[
                                                ncIVar_name]['_FillValue']

                ncOVar[:] = data.astype(dstBands[ncIVar_name]['type'])
                # copy (some) attributes
                for inAttrName in ncIVar.ncattrs():
                    if str(inAttrName) not in rmMetadata + ['dataType',
                                    'SourceFilename', 'SourceBand', '_Unsigned',
                                    'FillValue', 'time', '_FillValue']:
                        ncOVar.setncattr(inAttrName, ncIVar.getncattr(inAttrName))

                # add custom attributes from input parameter bands
                if ncIVar.name in bands:
                    for newAttr in bands[ncIVar.name]:
                        if newAttr not in rmMetadata + ['type', 'scale',
                                                        'offset',
                                                        '_FillValue']:
                            ncOVar._attributes[newAttr] = bands[ncIVar.name][newAttr]
                    # add grid_mapping info
                    if gridMappingName is not None:
                        ncOVar.setncattr('grid_mapping', gridMappingName)

        # copy (some) global attributes
        for globAttr in ncI.ncattrs():
            if not(globAttr.strip().startswith('GDAL')):
                ncO.setncattr(globAttr, ncI.getncattr(globAttr))

        # add common and custom global attributes
        ncO.setncatts(globMetadata)

        # write output file
        ncO.close()

        # close original files
        ncI.close()

        # Delete the temprary netCDF file
        os.close(fid)
        os.remove(tmpName)

        return 0

    def resize(self, factor=1, width=None, height=None,
               pixelsize=None, eResampleAlg=-1):
        '''Proportional resize of the dataset.

        The dataset is resized as (xSize*factor, ySize*factor)
        If desired width, height or pixelsize is specified,
        the scaling factor is calculated accordingly.
        If GCPs are given in a dataset, they are also rewritten.

        Parameters
        -----------
        factor : float, optional, default=1
            Scaling factor for width and height
            > 1 means increasing domain size
            < 1 means decreasing domain size
        width : int, optional
            Desired new width in pixels
        height : int, optional
            Desired new height in pixels
        pixelsize : float, optional
            Desired new pixelsize in meters (approximate).
            A factor is calculated from ratio of the
            current pixelsize to the desired pixelsize.
        eResampleAlg : int (GDALResampleAlg), optional
               -1 : Average (default),
                0 : NearestNeighbour
                1 : Bilinear,
                2 : Cubic,
                3 : CubicSpline,
                4 : Lancoz

        Modifies
        ---------
        self.vrt.dataset : VRT dataset of VRT object
            raster size are modified to downscaled size.
            If GCPs are given in the dataset, they are also overwritten.

        '''
        # get current shape
        rasterYSize = float(self.shape()[0])
        rasterXSize = float(self.shape()[1])

        # estimate factor if pixelsize is given
        if pixelsize is not None:
            deltaX, deltaY = self.get_pixelsize_meters()
            factorX = deltaX / float(pixelsize)
            factorY = deltaY / float(pixelsize)
            factor = (factorX + factorY)/2

        # estimate factor if width or height is given
        if width is not None:
            factor = float(width) / rasterXSize
        if height is not None:
            factor = float(height) / rasterYSize

        # calculate new size
        newRasterYSize = np.round(rasterYSize * factor)
        newRasterXSize = np.round(rasterXSize * factor)

        self.logger.info('New size/factor: (%f, %f)/%f' %
                         (newRasterXSize, newRasterYSize, factor))

        if eResampleAlg <= 0:
            self.vrt = self.vrt.get_subsampled_vrt(newRasterXSize,
                                                   newRasterYSize,
                                                   eResampleAlg)
        else:
            # update size and GeoTranform in XML of the warped VRT object
            self.vrt = self.vrt.get_resized_vrt(newRasterXSize,
                                                newRasterYSize,
                                                eResampleAlg=eResampleAlg)

        # resize gcps
        gcps = self.vrt.vrt.dataset.GetGCPs()
        if len(gcps) > 0:
            gcpPro = self.vrt.vrt.dataset.GetGCPProjection()
            for gcp in gcps:
                gcp.GCPPixel *= factor
                gcp.GCPLine *= factor
            self.vrt.dataset.SetGCPs(gcps, gcpPro)
            self.vrt._remove_geotransform()
        else:
            # change resultion in geotransform to keep spatial extent
            geoTransform = list(self.vrt.vrt.dataset.GetGeoTransform())
            geoTransform[1] = float(geoTransform[1])/factor
            geoTransform[5] = float(geoTransform[5])/factor
            geoTransform = map(float, geoTransform)
            self.vrt.dataset.SetGeoTransform(geoTransform)

        # set global metadata
        subMetaData = self.vrt.vrt.dataset.GetMetadata()
        subMetaData.pop('fileName')
        self.set_metadata(subMetaData)

        return factor

    def get_GDALRasterBand(self, bandID=1):
        ''' Get a GDALRasterBand of a given Nansat object

        If str is given find corresponding band number
        If int is given check if band with this number exists.
        Get a GDALRasterBand from vrt.

        Parameters
        -----------
        bandID : serial number or string, optional (default is 1)
            if number - a band number of the band to fetch
            if string bandID = {'name': bandID}

        Returns
        --------
        GDAL RasterBand

        Example
        -------
        b = n.get_GDALRasterBand(1)
        b = n.get_GDALRasterBand('sigma0')

        '''
        # get band number
        bandNumber = self._get_band_number(bandID)
        # the GDAL RasterBand of the corresponding band is returned
        return self.vrt.dataset.GetRasterBand(bandNumber)

    def list_bands(self, doPrint=True):
        ''' Show band information of the given Nansat object

        Show serial number, longName, name and all parameters
        for each band in the metadata of the given Nansat object.

        Parameters
        -----------
        doPrint : boolean, optional, default=True
            do print, otherwise it is returned as string

        Returns
        --------
        outString : String
            formatted string with bands info

        '''
        # get dictionary of bands metadata
        bands = self.bands()
        outString = ''

        for b in bands:
            # print band number, name
            outString += 'Band : %d %s\n' % (b, bands[b].get('name', ''))
            # print band metadata
            for i in bands[b]:
                outString += '  %s: %s\n' % (i, bands[b][i])
        if doPrint:
            # print to screeen
            print outString
        else:
            return outString

    def reproject(self, dstDomain, eResampleAlg=0, blockSize=None,
                  WorkingDataType=None, tps=None, skip_gcps=1, addmask=True,
                  **kwargs):
        ''' Change projection of the object based on the given Domain

        Create superVRT from self.vrt with AutoCreateWarpedVRT() using
        projection from the dstDomain.
        Modify XML content of the warped vrt using the Domain parameters.
        Generate warpedVRT and replace self.vrt with warpedVRT.
        If current object spans from 0 to 360 and dstDomain is west of 0,
        the object is shifted by 180 westwards.

        Parameters
        -----------
        dstDomain : domain
            destination Domain where projection and resolution are set
        eResampleAlg : int (GDALResampleAlg)
            0 : NearestNeighbour
            1 : Bilinear
            2 : Cubic,
            3 : CubicSpline
            4 : Lancoz
        blockSize : int
            size of blocks for resampling. Large value decrease speed
            but increase accuracy at the edge
        WorkingDataType : int (GDT_int, ...)
            type of data in bands. Shuold be integer for int32 bands
        tps : bool
            Apply Thin Spline Transformation if source or destination has GCPs
            Usage of TPS can also be triggered by setting self.vrt.tps=True
            before calling to reproject.
            This options has priority over self.vrt.tps
        skip_gcps : int
            Using TPS can be very slow if the number of GCPs are large.
            If this parameter is given, only every [skip_gcp] GCP is used,
            improving calculation time at the cost of accuracy.
            If not given explicitly, 'skip_gcps' is fetched from the
            metadata of self, or from dstDomain (as set by mapper or user).
            [defaults to 1 if not specified, i.e. using all GCPs]
        addmask : bool
            If True, add band 'swathmask'. 1 - valid data, 0 no-data.
            This band is used to replace no-data values with np.nan

        Modifies
        ---------
        self.vrt : VRT object with dataset replaced to warpedVRT dataset

        !! NB !!
        ---------
        - Integer data is returnd by integer. Round off to decimal place.
          If you do not want to round off, convert the data types
          to GDT_Float32, GDT_Float64, or GDT_CFloat32.

        See Also
        ---------
        http://www.gdal.org/gdalwarp.html
        '''
        # if self spans from 0 to 360 and dstDomain is west of 0:
        #     shift self westwards by 180 degrees
        # check span
        srcCorners = self.get_corners()
        if round(min(srcCorners[0])) == 0 and round(max(srcCorners[0])) == 360:
            # check intersection of src and dst
            dstCorners = dstDomain.get_corners()
            if min(dstCorners[0]) < 0:
                # shift
                self.vrt = self.vrt.get_shifted_vrt(-180)

        # get projection of destination dataset
        dstSRS = dstDomain.vrt.dataset.GetProjection()

        # get destination GCPs
        dstGCPs = dstDomain.vrt.dataset.GetGCPs()
        if len(dstGCPs) > 0:
            # get projection of destination GCPs
            dstSRS = dstDomain.vrt.dataset.GetGCPProjection()

        xSize = dstDomain.vrt.dataset.RasterXSize
        ySize = dstDomain.vrt.dataset.RasterYSize

        # get geoTransform
        if 'use_gcps' in kwargs.keys() and not (kwargs['use_gcps']):
            corners = dstDomain.get_corners()
            ext = '-lle %0.3f %0.3f %0.3f %0.3f -ts %d %d' % (min(corners[0]),
                                                              min(corners[1]),
                                                              max(corners[0]),
                                                              max(corners[1]),
                                                              xSize, ySize)
            d = Domain(srs=dstSRS, ext=ext)
            geoTransform = d.vrt.dataset.GetGeoTransform()
        else:
            geoTransform = dstDomain.vrt.dataset.GetGeoTransform()

        # set trigger for using TPS
        if tps is True:
            self.vrt.tps = True
        elif tps is False:
            self.vrt.tps = False

        # Reduce number of GCPs for faster reprojection
        # when using TPS (if requested)
        src_skip_gcps = self.vrt.dataset.GetMetadataItem('skip_gcps')
        dst_skip_gcps = dstDomain.vrt.dataset.GetMetadataItem('skip_gcps')
        kwargs['skip_gcps'] = skip_gcps  # default (use all GCPs)
        if dst_skip_gcps is not None:  # ...or use setting from dst
            kwargs['skip_gcps'] = int(dst_skip_gcps)
        if src_skip_gcps is not None:  # ...or use setting from src
            kwargs['skip_gcps'] = int(src_skip_gcps)

        # add band that masks valid values with 1 and nodata with 0
        # after reproject
        if addmask:
            self.vrt = self.vrt.get_super_vrt()
            self.vrt._create_band(
                src=[{
                    'SourceFilename': self.vrt.vrt.fileName,
                    'SourceBand':  1,
                    'DataType': gdal.GDT_Byte}],
                dst={
                    'dataType': gdal.GDT_Byte,
                    'wkv': 'swath_binary_mask',
                    'PixelFunctionType': 'OnesPixelFunc',
                })
            self.vrt.dataset.FlushCache()

        # create Warped VRT
        self.vrt = self.vrt.get_warped_vrt(dstSRS=dstSRS,
                                           dstGCPs=dstGCPs,
                                           eResampleAlg=eResampleAlg,
                                           xSize=xSize, ySize=ySize,
                                           blockSize=blockSize,
                                           geoTransform=geoTransform,
                                           WorkingDataType=WorkingDataType,
                                           **kwargs)

        # set global metadata from subVRT
        subMetaData = self.vrt.vrt.dataset.GetMetadata()
        subMetaData.pop('fileName')
        self.set_metadata(subMetaData)

    def undo(self, steps=1):
        '''Undo reproject, resize, add_band or crop of Nansat object

        Restore the self.vrt from self.vrt.vrt

        Parameters
        -----------
        steps : int
            How many steps back to undo

        Modifies
        --------
        self.vrt

        '''

        self.vrt = self.vrt.get_sub_vrt(steps)

    def watermask(self, mod44path=None, dstDomain=None, **kwargs):
        ''' Create numpy array with watermask (water=1, land=0)

        250 meters resolution watermask from MODIS 44W Product:
        http://www.glcf.umd.edu/data/watermask/

        Watermask is stored as tiles in TIF(LZW) format and a VRT file
        All files are stored in one directory.
        A tarball with compressed TIF and VRT files should be additionally
        downloaded from the Nansat wiki:
        https://svn.nersc.no/nansat/wiki/Nansat/Data/Watermask

        The method :
            Gets the directory either from input parameter or from environment
            variable MOD44WPATH
            Open Nansat object from the VRT file
            Reprojects the watermask onto the current object using reproject()
            or reproject_on_jcps()
            Returns the reprojected Nansat object

        Parameters
        -----------
        mod44path : string, optional, default=None
            path with MOD44W Products and a VRT file
        dstDomain : Domain
            destination domain other than self
        tps : Bool
            Use Thin Spline Transformation in reprojection of watermask?
            See also Nansat.reproject()
        skip_gcps : int
            Factor to reduce the number of GCPs by and increase speed
            See also Nansat.reproject()

        Returns
        --------
        watermask : Nansat object with water mask in current projection

        See also
        ---------
        250 meters resolution watermask from MODIS 44W Product:
            http://www.glcf.umd.edu/data/watermask/

        '''
        mod44DataExist = True
        # check if path is given in input param or in environment
        if mod44path is None:
            mod44path = os.getenv('MOD44WPATH')
        if mod44path is None:
            mod44DataExist = False
        # check if VRT file exist
        elif not os.path.exists(mod44path + '/MOD44W.vrt'):
            mod44DataExist = False
        self.logger.debug('MODPATH: %s' % mod44path)

        if not mod44DataExist:
            raise IOError('250 meters resolution watermask from MODIS '
                    '44W Product does not exist - see Nansat '
                    'documentation to get it (the path is % s)' % mod44path)

        # MOD44W data does exist: open the VRT file in Nansat
        watermask = Nansat(mod44path + '/MOD44W.vrt', mapperName='MOD44W',
                           logLevel=self.logger.level)
        # reproject on self or given Domain
        if dstDomain is None:
            dstDomain = self
        lon, lat = dstDomain.get_border()
        watermask.crop_lonlat([lon.min(), lon.max()], [lat.min(), lat.max()])
        watermask.reproject(dstDomain, addmask=False, **kwargs)

        return watermask

    def write_figure(self, fileName=None, bands=1, clim=None, addDate=False,
                     array_modfunc=None, **kwargs):
        ''' Save a raster band to a figure in graphical format.

        Get numpy array from the band(s) and band information specified
        either by given band number or band id.
        -- If three bands are given, merge them and create PIL image.
        -- If one band is given, create indexed image
        Create Figure object and:
        Adjust the array brightness and contrast using the given min/max or
        histogram.
        Apply logarithmic scaling of color tone.
        Generate and append legend.
        Save the PIL output image in PNG or any other graphical format.
        If the filename extension is 'tif', the figure file is converted
        to GeoTiff

        Parameters
        -----------
        fileName : string, optional
            Output file name. if one of extensions 'png', 'PNG', 'tif',
            'TIF', 'bmp', 'BMP', 'jpg', 'JPG', 'jpeg', 'JPEG' is included,
            specified file is crated. otherwise, 'png' file is created.
            if None, the figure object is returned.
            if True, the figure is shown
        bands : integer or string or list (elements are integer or string),
            default = 1
            the size of the list has to be 1 or 3.
            if the size is 3, RGB image is created based on the three bands.
            Then the first element is Red, the second is Green,
            and the third is Blue.
        clim : list with two elements or 'hist' to specify range of colormap
            None (default) : min/max values are fetched from WKV,
            fallback-'hist'
            [min, max] : min and max are numbers, or
            [[min, min, min], [max, max, max]]: three bands used
            'hist' : a histogram is used to calculate min and max values
        addDate : boolean
            False (default) : no date will be aded to the caption
            True : the first time of the object will be added to the caption
        array_modfunc : None
            None (default) : figure created using array in provided band
            function : figure created using array modified by provided function
        **kwargs : parameters for Figure().

        Modifies
        ---------
        if fileName is specified, creates image file

        Returns
        -------
        Figure object

        Example
        --------
        #write only indexed image, color limits from WKV or from histogram
        n.write_figure('test.jpg')
        #write only RGB image, color limits from histogram
        n.write_figure('test_rgb_hist.jpg', clim='hist', bands=[1, 2, 3])
        #write indexed image, apply log scaling and gamma correction,
        #add legend and type in title 'Title', increase font size and put 15
        tics
        n.write_figure('r09_log3_leg.jpg', logarithm=True, legend=True,
                                gamma=3, titleString='Title', fontSize=30,
                                numOfTicks=15)
        # write an image to png with transparent Mask set to color
        transparency=[0,0,0], following PIL alpha mask
        n.write_figure(fileName='transparent.png', bands=[3],
               mask_array=wmArray,
               mask_lut={0: [0,0,0]},
               clim=[0,0.15], cmapName='gray', transparency=[0,0,0])

        See also
        --------
        Figure()
        http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps

        '''
        # convert <bands> from integer, or string, or list of strings
        # into list of integers
        if isinstance(bands, list):
            for i, band in enumerate(bands):
                bands[i] = self._get_band_number(band)
        else:
            bands = [self._get_band_number(bands)]

        # == create 3D ARRAY ==
        array = None
        for band in bands:
            # get array from band and reshape to (1,height,width)
            iArray = self[band]
            if array_modfunc:
                iArray = array_modfunc(iArray)
            iArray = iArray.reshape(1, iArray.shape[0], iArray.shape[1])
            # create new 3D array or append band
            if array is None:
                array = iArray
            else:
                array = np.append(array, iArray, axis=0)

        # == CREATE FIGURE object and parse input parameters ==
        fig = Figure(array, **kwargs)
        array = None

        # == PREPARE cmin/cmax ==
        # check if cmin and cmax are given as the arguments
        if 'cmin' in kwargs.keys() and 'cmax' in kwargs.keys():
            clim = [kwargs['cmin'], kwargs['cmax']]

        # try to get clim from WKV if it is not given as the argument
        # if failed clim will be evaluated from histogram
        if clim is None:
            clim = [[], []]
            for i, iBand in enumerate(bands):
                try:
                    defValue = (self.vrt.dataset.GetRasterBand(iBand).
                                GetMetadataItem('minmax').split(' '))
                except:
                    clim = 'hist'
                    break
                clim[0].append(float(defValue[0]))
                clim[1].append(float(defValue[1]))

        # Estimate color min/max from histogram
        if clim == 'hist':
            clim = fig.clim_from_histogram(**kwargs)

        # modify clim to the proper shape [[min], [max]]
        # or [[min, min, min], [max, max, max]]
        if (len(clim) == 2 and
           ((isinstance(clim[0], float)) or (isinstance(clim[0], int))) and
           ((isinstance(clim[1], float)) or (isinstance(clim[1], int)))):
            clim = [[clim[0]], [clim[1]]]

        # if the len(clim) is not same as len(bands), the 1st element is used.
        for i in range(2):
            if len(clim[i]) != len(bands):
                clim[i] = [clim[i][0]] * len(bands)

        self.logger.info('clim: %s ' % clim)

        # == PREPARE caption ==
        if 'caption' in kwargs:
            caption = kwargs['caption']
        else:
            # get longName and units from vrt
            band = self.get_GDALRasterBand(bands[0])
            longName = band.GetMetadata().get('long_name', '')
            units = band.GetMetadata().get('units', '')

            # make caption from longname, units
            caption = longName + ' [' + units + ']'

        # add DATE to caption
        if addDate:
            caption += self.time_coverage_start.strftime(' %Y-%m-%d')

        self.logger.info('caption: %s ' % caption)

        # == PROCESS figure ==
        fig.process(cmin=clim[0], cmax=clim[1], caption=caption)

        # == finally SAVE to a image file or SHOW ==
        if fileName is not None:
            if type(fileName) == bool and fileName:
                try:
                    if plt.get_backend() == 'agg':
                        plt.switch_backend('QT4Agg')
                except:
                    fig.pilImg.show()
                else:
                    sz = fig.pilImg.size
                    imgArray = np.array(fig.pilImg.im)
                    if fig.pilImg.getbands() == ('P',):
                        imgArray.resize(sz[1], sz[0])
                    elif fig.pilImg.getbands() == ('R', 'G', 'B'):
                        imgArray.resize(sz[1], sz[0], 3)
                    plt.imshow(imgArray)
                    plt.show()

            elif type(fileName) in [str, unicode]:
                fig.save(fileName, **kwargs)
                # If tiff image, convert to GeoTiff
                if fileName[-3:] == 'tif':
                    self.vrt.copyproj(fileName)
            else:
                raise OptionError('%s is of wrong type %s' %
                                  (str(fileName), str(type(fileName))))
        return fig

    def write_geotiffimage(self, fileName, bandID=1):
        ''' Writes an 8-bit GeoTiff image for a given band.

        The output GeoTiff image is convenient e.g. for display in a GIS tool.
        Colormap is fetched from the metadata item 'colormap'.
            Fallback colormap is 'jet'.
        Color limits are fetched from the metadata item 'minmax'.
            If 'minmax' is not specified, min and max of raster is used.

        The method can be replaced by using nansat.write_figure(),
        however, write_figure uses PIL which does not allow
        Tiff compression, giving much larger files

        Parameters
        -----------
        fileName : string
        bandID : integer or string(default = 1)

        '''
        bandNo = self._get_band_number(bandID)
        band = self.get_GDALRasterBand(bandID)
        minmax = band.GetMetadataItem('minmax')
        # Get min and max from band histogram if not given (from wkv)
        if minmax is None:
            (rmin, rmax) = band.ComputeRasterMinMax(1)
            minmax = str(rmin) + ' ' + str(rmax)

        bMin = float(minmax.split(' ')[0])
        bMax = float(minmax.split(' ')[1])
        # Make colormap from WKV information
        try:
            colormap = band.GetMetadataItem('colormap')
        except:
            colormap = 'jet'
        cmap = cm.get_cmap(colormap, 256)
        cmap = cmap(np.arange(256)) * 255
        colorTable = gdal.ColorTable()
        for i in range(cmap.shape[0]):
            colorEntry = (int(cmap[i, 0]), int(cmap[i, 1]),
                          int(cmap[i, 2]), int(cmap[i, 3]))
            colorTable.SetColorEntry(i, colorEntry)
        # Write Tiff image, with data scaled to values between 0 and 255
        outDataset = gdal.GetDriverByName('Gtiff').Create(fileName,
                                                          band.XSize,
                                                          band.YSize, 1,
                                                          gdal.GDT_Byte,
                                                          ['COMPRESS=LZW'])
        data = self.__getitem__(bandNo)
        scaledData = ((data - bMin) / (bMax - bMin)) * 255
        outDataset.GetRasterBand(1).WriteArray(scaledData)
        outDataset.GetRasterBand(1).SetMetadata(band.GetMetadata())
        try:
            outDataset.GetRasterBand(1).SetColorTable(colorTable)
        except:
            # Happens after reprojection, a possible bug?
            print 'Could not set color table'
            print colorTable
        outDataset = None
        self.vrt.copyproj(fileName)

    @property
    def time_coverage_start(self):
        return parse_time(self.get_metadata('time_coverage_start'))

    @property
    def time_coverage_end(self):
        return parse_time(self.get_metadata('time_coverage_end'))

    def get_metadata(self, key=None, bandID=None):
        ''' Get metadata from self.vrt.dataset

        Parameters
        ----------
        key : string, optional
            name of the metadata key. If not givem all metadata is returned
        bandID : int or str, optional
            number or name of band to get metadata from.
            If not given, global metadata is returned

        Returns
        --------
        a string with metadata if key is given and found
        an empty string if key is given and not found
        a dictionary with all metadata if key is not given

        '''
        # get all metadata from dataset or from band
        if bandID is None:
            metadata = self.vrt.dataset.GetMetadata()
        else:
            metadata = self.get_GDALRasterBand(bandID).GetMetadata()

        # get all metadata or from a key
        if key is not None:
            try:
                metadata = metadata[key]
            except KeyError:
                raise OptionError('%s does not have metadata %s' % (
                                   self.fileName, key))

        return metadata

    def set_metadata(self, key='', value='', bandID=None):
        ''' Set metadata to self.vrt.dataset

        Parameters
        -----------
        key : string or dictionary with strings
            name of the metadata, or dictionary with metadata names, values
        value : string
            value of metadata
        bandID : int or str
            number or name of band
            Without : global metadata is set

        Modifies
        ---------
        self.vrt.dataset : sets metadata in GDAL current dataset

        '''
        # set all metadata to the dataset or to the band
        if bandID is None:
            metaReceiverVRT = self.vrt.dataset
        else:
            bandNumber = self._get_band_number(bandID)
            metaReceiverVRT = self.vrt.dataset.GetRasterBand(bandNumber)

        # set metadata from dictionary or from single pair key,value
        if type(key) == dict:
            for k in key:
                metaReceiverVRT.SetMetadataItem(k, key[k])
        else:
            metaReceiverVRT.SetMetadataItem(key, value)

    def _get_mapper(self, mapperName, **kwargs):
        ''' Create VRT file in memory (VSI-file) with variable mapping

        If mapperName is given only this mapper will be used,
        else loop over all availble mappers in mapperList to get the
        matching one.
        In the loop :
            If the specific error appears the mapper is not used
            and the next mapper is tested.
            Otherwise the mapper returns VRT.
        If type of the sensor is identified, add mapping variables.
        If all mappers fail, make simple copy of the input DS into a VSI/VRT

        Parameters
        -----------
        mapperName : string, optional (e.g. 'ASAR' or 'merisL2')

        Returns
        --------
        tmpVRT : VRT object
            tmpVRT.dataset is a GDAL VRT dataset

        Raises
        --------
        IOError : occurs if the input file does not exist
        OptionError : occurs if given mapper cannot open the input file
        NansatReadError : occurs if no mapper fits the input file

        '''
        if os.path.isfile(self.fileName):
            # Make sure file exists and can be opened for reading
            # before proceeding
            test_openable(self.fileName)
        else:
            ff = glob.glob(os.path.join(self.fileName, '*.*'))
            for f in ff:
                test_openable(f)
        # lazy import of nansat mappers
        # if nansat mappers were not imported yet
        global nansatMappers
        if nansatMappers is None:
            nansatMappers = _import_mappers()

        # open GDAL dataset. It will be parsed to all mappers for testing
        gdalDataset = None
        if self.fileName[:4] != 'http':
            try:
                gdalDataset = gdal.Open(self.fileName)
            except RuntimeError:
                self.logger.error('GDAL could not open ' + self.fileName +
                                  ', trying to read with Nansat mappers...')
        if gdalDataset is not None:
            # get metadata from the GDAL dataset
            metadata = gdalDataset.GetMetadata()
        else:
            metadata = None

        tmpVRT = None

        importErrors = []
        if mapperName is not '':
            # If a specific mapper is requested, we test only this one.
            # get the module name
            mapperName = 'mapper_' + mapperName.replace('mapper_',
                                                        '').replace('.py',
                                                                    '').lower()
            # check if the mapper is available
            if mapperName not in nansatMappers:
                raise OptionError('Mapper ' + mapperName + ' not found')

            # check if mapper is importbale or raise an ImportError error
            if isinstance(nansatMappers[mapperName], tuple):
                errType, err, traceback = nansatMappers[mapperName]
                # self.logger.error(err, exc_info=(errType, err, traceback))
                raise errType, err, traceback

            # create VRT using the selected mapper
            tmpVRT = nansatMappers[mapperName](self.fileName,
                                               gdalDataset,
                                               metadata,
                                               **kwargs)
            self.mapper = mapperName.replace('mapper_', '')
        else:
            # We test all mappers, import one by one
            for iMapper in nansatMappers:
                # skip non-importable mappers
                if isinstance(nansatMappers[iMapper], tuple):
                    # keep errors to show before use of generic mapper
                    importErrors.append(nansatMappers[iMapper][1])
                    continue

                self.logger.debug('Trying %s...' % iMapper)

                # show all ImportError warnings before trying generic_mapper
                if iMapper == 'mapper_generic' and len(importErrors) > 0:
                    self.logger.error('\nWarning! '
                                      'The following mappers failed:')
                    for ie in importErrors:
                        self.logger.error(importErrors)

                # create a Mapper object and get VRT dataset from it
                try:
                    tmpVRT = nansatMappers[iMapper](self.fileName,
                                                    gdalDataset,
                                                    metadata,
                                                    **kwargs)
                    self.logger.info('Mapper %s - success!' % iMapper)
                    self.mapper = iMapper.replace('mapper_', '')
                    break
                except WrongMapperError:
                    pass

        # if no mapper fits, make simple copy of the input DS into a VSI/VRT
        if tmpVRT is None and gdalDataset is not None:
            self.logger.warning('No mapper fits, returning GDAL bands!')
            tmpVRT = VRT(gdalDataset=gdalDataset)
            for iBand in range(gdalDataset.RasterCount):
                tmpVRT._create_band({'SourceFilename': self.fileName,
                                     'SourceBand': iBand + 1})
                tmpVRT.dataset.FlushCache()
            self.mapper = 'gdal_bands'

        # if GDAL cannot open the file, and no mappers exist which can make VRT
        if tmpVRT is None and gdalDataset is None:
            # check if given data file exists
            if not os.path.isfile(self.fileName):
                raise IOError('%s: File does not exist' % (self.fileName))
            raise NansatReadError('%s: File cannot be read with NANSAT - '
                    'consider writing a mapper' % self.fileName)

        return tmpVRT

    def _get_pixelValue(self, val, defVal):
        if val == '':
            return defVal
        else:
            return val

    def _get_band_number(self, bandID):
        '''Return absolute band number

        Check if given bandID is valid
        Return absolute number of the band in the VRT

        Parameters
        ----------
        bandID : int or str or dict
            if int : checks if such band exists and returns band_id
            if str : finds band with coresponding name
            if dict : finds first band with given metadata

        Returns
        --------
        int : absolute band number

        '''
        bandNumber = 0
        # if bandID is str: create simple dict with seraching criteria
        if type(bandID) == str:
            bandID = {'name': bandID}

        # if bandID is dict: search self.bands with seraching criteria
        if type(bandID) == dict:
            bandsMeta = self.bands()
            for b in bandsMeta:
                numCorrectKeys = 0
                for key in bandID:
                    if (key in bandsMeta[b] and
                            bandID[key] == bandsMeta[b][key]):
                        numCorrectKeys = numCorrectKeys + 1
                    if numCorrectKeys == len(bandID):
                        bandNumber = b
                        break

        # if bandID is int and with bounds: return this number
        if (type(bandID) == int and bandID >= 1 and
                bandID <= self.vrt.dataset.RasterCount):
            bandNumber = bandID

        # if no bandNumber found - raise error
        if bandNumber == 0:
            raise OptionError('Cannot find band %s! '
                              'bandNumber is from 1 to %s'
                              % (str(bandID), self.vrt.dataset.RasterCount))

        return bandNumber

    def get_transect(self, points, bands,
                        lonlat=True,
                        smoothRadius=0,
                        smooth_function=nanmedian,
                        data=None,
                        cornersonly=False):
        '''Get values from transect from given vector of poins

        Parameters
        ----------
        points : 2xN list or array, N (number of points) >= 1
            coordinates [[x1, x2, y2], [y1, y2, y3]]
        bands : list of int or string
            elements of the list are band number or band Name
        lonlat : bool
            If the points in lat/lon, then True.
            If the points in pixel/line, then False.
        smoothRadius: int
            If smootRadius is greater than 0, smooth every transect
            pixel as the median or mean value in a circule with radius
            equal to the given number.
        smooth_function: func
            function for averaging values collected within smooth radius
        data : ndarray
            alternative array with data to take values from

        Returns
        --------
        transect : numpy record array

        '''
        # check if points is 2D array with shape 2xN (N>=1)
        if (len(np.shape(points)) != 2 or
              np.shape(points)[0] != 2 or
              np.shape(points)[1] < 1):
            # points are not 2xN array
            raise OptionError('Input points must be 2xN array with N>0')

        # get names of bands
        bandNames = []
        for band in bands:
            try:
                bandN = self._get_band_number(band)
            except OptionError:
                self.logger.error('Wrong band name %s' % band)
            else:
                bandNames.append(self.bands()[bandN]['name'])

        if data is not None:
            bandNames.append('input')

        # if points in degree, convert them into pix/lin
        if lonlat:
            pix, lin = self.transform_points(points[0], points[1], DstToSrc=1)
        else:
            pix, lin = points[0], points[1]

        if cornersonly:
            pixVector, linVector = pix, lin
        else:
            # full vectors of pixel coordinates based on coordinates of vertices
            pixVector, linVector = [pix[0]], [lin[0]]
            for pn in range(len(pix[1:])):
                px0, px1 = pix[pn], pix[pn+1]
                py0, py1 = lin[pn], lin[pn+1]
                length = np.round(np.hypot(px1-px0, py0-py1))
                pixVector += list(np.linspace(px0, px1, length+1)[1:])
                linVector += list(np.linspace(py0, py1, length+1)[1:])

            # remove out of region points
            pixVector = np.floor(pixVector)
            linVector = np.floor(linVector)
            gpi = ((pixVector >= (0 + smoothRadius)) *
                   (linVector >= (0 + smoothRadius)) *
                   (pixVector < (self.shape()[1] - smoothRadius)) *
                   (linVector < (self.shape()[0] - smoothRadius)))
            pixVector = pixVector[gpi]
            linVector = linVector[gpi]

        # create output transect
        t = np.recarray((len(pixVector)), dtype=[('pixel', int),
                                                ('line', int),
                                                ('lon', float),
                                                ('lat', float), ])

        # add pixel, line, lon, lat values to output
        t['pixel'] = pixVector
        t['line'] = linVector
        t['lon'], t['lat'] = self.transform_points(t['pixel'], t['line'],
                                                   DstToSrc=0)

        # mask for extraction within circular area
        xgrid, ygrid = np.mgrid[0:smoothRadius * 2 + 1, 0:smoothRadius * 2 + 1]
        distance = ((xgrid - smoothRadius) ** 2 +
                    (ygrid - smoothRadius) ** 2) ** 0.5
        mask = distance <= smoothRadius

        # get values from bands or input data
        if len(bandNames) > 0:
            for bandName in bandNames:
                if bandName == 'input':
                    bandArray = data
                else:
                    bandArray = self[bandName]
                # average values from pixel inside a circle
                bandValues = []
                for r, c in zip(t['line'], t['pixel']):
                    subarray = bandArray[r-smoothRadius:r+smoothRadius+1,
                                         c-smoothRadius:c+smoothRadius+1]
                    bandValues.append(smooth_function(subarray[mask]))
                t = append_fields(t, bandName, bandValues).data

        return t

    def digitize_points(self, band=1, **kwargs):

        '''Get coordinates of interactively digitized points

        Parameters
        ----------
        band : int or str
            ID of Nansat band
        **kwargs : keyword arguments for imshow

        Returns
        --------
        points : list
            list of 2xN arrays of points to be used in Nansat.get_transect()

        '''
        if matplotlib.is_interactive():
            warnings.warn('''
        Python is started with -pylab option, transect will not work.
        Please restart python without -pylab.''')
            return []

        data = self[band]
        browser = PointBrowser(data, **kwargs)
        points = browser.get_points()

        return points

    def crop_interactive(self, band=1,**kwargs):
        ''' Interactively select boundary and crop Nansat object

        Parameters
        ----------
        band : int or str
            id of the band to show for interactive selection of boundaries
        **kwargs : keyword arguments for imshow

        Modifies
        --------
        self.vrt : VRT
            superVRT is created with modified SrcRect and DstRect
        Returns
        -------
        extent : (xOff, yOff, xSize, ySize)
            xOff  - X offset in the original dataset
            yOff  - Y offset in the original dataset
            xSize - width of the new dataset
            ySize - height of the new dataset

        Examples
        --------
        # crop a subimage interactively
        from matplotlib import cm
        extent = n.crop_interactive(band=1,cmap=cm.gray)

        '''
        maxwidth = 1000
        resized = False
        if self.shape()[1] > maxwidth:
            factor = self.resize(width=1000)
            resized = True
        else:
            factor = 1
        # use interactive PointBrowser for selecting extent
        try:
           points = self.digitize_points(band=band,**kwargs)[0]
        except:
           if resized:
              self.undo()
           return

        xOff = round(points.min(axis=1)[0] / factor)
        yOff = round(points.min(axis=1)[1] / factor)
        xSize = round((points.max(axis=1)[0] - xOff) / factor)
        ySize = round((points.max(axis=1)[1] - yOff) / factor)
        if resized:
            self.undo()

        return self.crop(xOff, yOff, xSize, ySize)

    def crop_lonlat(self, lonlim, latlim):
        ''' Crop Nansat object to fit into given longitude/latitude limit
        Modifies
        --------
        self.vrt : VRT
            superVRT is created with modified SrcRect and DstRect

        Returns
        -------
        extent : (xOff, yOff, xSize, ySize)
            xOff  - X offset in the original dataset
            yOff  - Y offset in the original dataset
            xSize - width of the new dataset
            ySize - height of the new dataset

        Examples
        --------
        # crop a subimage for given lon/lat limits
        extent = n.crop(lonlim=[-10,10], latlim=[-20,20])

        '''
        crnPix, crnLin = self.transform_points([lonlim[0], lonlim[0],
                                                lonlim[1], lonlim[1]],
                                               [latlim[0], latlim[1],
                                                latlim[0], latlim[1]], 1)
        xOff = round(min(crnPix))
        yOff = round(min(crnLin))
        xSize = round(max(crnPix) - min(crnPix))
        ySize = round(max(crnLin) - min(crnLin))

        return self.crop(xOff, yOff, xSize, ySize)

    def crop(self, xOff, yOff, xSize, ySize):
        '''Crop Nansat object

        Create superVRT, modify the Source Rectangle (SrcRect) and Destination
        Rectangle (DstRect) tags in the VRT file for each band in order
        to take only part of the original image,
        create new GCPs or new GeoTransform for the cropped object.

        Parameters
        ----------
        xOff : int
            pixel offset of subimage
        yOff : int
            line offset of subimage
        xSize : int
            width in pixels of subimage
        ySize : int
            height in pizels of subimage

        Modifies
        --------
        self.vrt : VRT
            superVRT is created with modified SrcRect and DstRect
        Returns
        -------
        extent : (xOff, yOff, xSize, ySize)
            xOff  - X offset in the original dataset
            yOff  - Y offset in the original dataset
            xSize - width of the new dataset
            ySize - height of the new dataset

        Examples
        --------
            # crop a subimage of size 100x200 pix from X/Y offset 10, 20 pix
            extent = n.crop(10, 20, 100, 200)

        '''
        RasterXSize = self.vrt.dataset.RasterXSize
        RasterYSize = self.vrt.dataset.RasterYSize

        # set xSize/ySize if ommited in the call
        if xSize is None:
            xSize = RasterXSize - xOff
        if ySize is None:
            ySize = RasterYSize - yOff

        # test if crop is totally outside
        if (xOff > RasterXSize or (xOff + xSize) < 0 or
                yOff > RasterYSize or (yOff + ySize) < 0):
            raise OptionError('''Cropping region is outside the image!
                               xOff: %d, yOff: %d, xSize: %d, ySize: %d''' %
                              (xOff,  yOff, xSize, ySize))

        # set default values of invalud xOff/yOff and xSize/ySize
        if xOff < 0:
            xSize += xOff
            xOff = 0

        if yOff < 0:
            ySize += yOff
            yOff = 0

        if (xSize + xOff) > RasterXSize:
            xSize = RasterXSize - xOff
        if (ySize + yOff) > RasterYSize:
            ySize = RasterYSize - yOff

        extent = (int(xOff), int(yOff), int(xSize), int(ySize))
        self.logger.debug('xOff: %d, yOff: %d, xSize: %d, ySize: %d' % extent)

        # test if crop is too large
        if (xOff == 0 and xSize == RasterXSize and
                yOff == 0 and ySize == RasterYSize):
            self.logger.error(('WARNING! Cropping region is'
                               'larger or equal to image!'))
            return extent

        # create super VRT and get its XML
        self.vrt = self.vrt.get_super_vrt()
        xml = self.vrt.read_xml()
        node0 = Node.create(xml)

        # change size
        node0.node('VRTDataset').replaceAttribute('rasterXSize', str(xSize))
        node0.node('VRTDataset').replaceAttribute('rasterYSize', str(ySize))

        # replace x/y-Off and x/y-Size
        #   in <SrcRect> and <DstRect> of each source
        for iNode1 in node0.nodeList('VRTRasterBand'):
            iNode2 = iNode1.node('ComplexSource')

            iNode3 = iNode2.node('SrcRect')
            iNode3.replaceAttribute('xOff', str(xOff))
            iNode3.replaceAttribute('yOff', str(yOff))
            iNode3.replaceAttribute('xSize', str(xSize))
            iNode3.replaceAttribute('ySize', str(ySize))

            iNode3 = iNode2.node('DstRect')
            iNode3.replaceAttribute('xSize', str(xSize))
            iNode3.replaceAttribute('ySize', str(ySize))

        # write modified XML
        xml = node0.rawxml()
        self.vrt.write_xml(xml)

        # modify GCPs or GeoTranfrom to fit the new shape of image
        gcps = self.vrt.dataset.GetGCPs()
        if len(gcps) > 0:
            dstGCPs = []
            i = 0
            # keep current GCPs
            for igcp in gcps:
                if (0 < igcp.GCPPixel - xOff and
                        igcp.GCPPixel - xOff < xSize and
                        0 < igcp.GCPLine - yOff and
                        igcp.GCPLine - yOff < ySize):
                    i += 1
                    dstGCPs.append(gdal.GCP(igcp.GCPX, igcp.GCPY, 0,
                                            igcp.GCPPixel - xOff,
                                            igcp.GCPLine - yOff, '', str(i)))
            numOfGCPs = i

            if numOfGCPs < 100:
                # create new 100 GPCs (10 x 10 regular matrix)
                pixArray = []
                linArray = []
                for newPix in np.r_[0:xSize:10j]:
                    for newLin in np.r_[0:ySize:10j]:
                        pixArray.append(newPix + xOff)
                        linArray.append(newLin + yOff)

                lonArray, latArray = self.vrt.transform_points(pixArray,
                                                               linArray,
                        dstSRS=NSR(self.vrt.dataset.GetGCPProjection()))

                for i in range(len(lonArray)):
                    dstGCPs.append(gdal.GCP(lonArray[i], latArray[i], 0,
                                            pixArray[i] - xOff,
                                            linArray[i] - yOff,
                                            '', str(numOfGCPs+i+1)))

            # set new GCPss
            self.vrt.dataset.SetGCPs(dstGCPs,
                                     self.vrt.dataset.GetGCPProjection())
            # remove geotranform which was automatically added
            self.vrt._remove_geotransform()
        else:
            # shift upper left corner coordinates
            geoTransfrom = self.vrt.dataset.GetGeoTransform()
            geoTransfrom = map(float, geoTransfrom)
            geoTransfrom[0] += geoTransfrom[1] * xOff
            geoTransfrom[3] += geoTransfrom[5] * yOff
            self.vrt.dataset.SetGeoTransform(geoTransfrom)

        # set global metadata
        subMetaData = self.vrt.vrt.dataset.GetMetadata()
        subMetaData.pop('fileName')
        self.set_metadata(subMetaData)

        return extent


def _import_mappers(logLevel=None):
    ''' Import available mappers into a dictionary

    Returns
    --------
    nansatMappers : dict
        key  : mapper name
        value: class Mapper(VRT) from the mapper module

    '''
    logger = add_logger('import_mappers', logLevel=logLevel)
    # import built-in mappers
    import nansat.mappers
    mappersPackages = [nansat.mappers]

    # import user-defined mappers (if any)
    try:
        import nansat_mappers
    except:
        pass
    else:
        logger.info('User defined mappers found in %s'
                    % nansat_mappers.__path__)
        mappersPackages = [nansat_mappers, nansat.mappers]

    # create ordered dict for mappers
    nansatMappers = OrderedDict()

    for mappersPackage in mappersPackages:
        logger.debug('From package: %s' % mappersPackage.__path__)
        # scan through modules and load all modules that contain class Mapper
        for finder, name, ispkg in (pkgutil.
                                    iter_modules(mappersPackage.__path__)):
            logger.debug('Loading mapper %s' % name)
            loader = finder.find_module(name)
            # try to import mapper module
            try:
                module = loader.load_module(name)
            except ImportError:
                # keep ImportError instance instead of the mapper
                exc_info = sys.exc_info()
                logger.error('Mapper %s could not be imported'
                             % name, exc_info=exc_info)
                nansatMappers[name] = exc_info
            else:
                # add the imported mapper to nansatMappers
                if hasattr(module, 'Mapper'):
                    nansatMappers[name] = module.Mapper

        # move netcdfcdf mapper to the end
        if 'mapper_netcdfcf' in nansatMappers:
            nansatMappers['mapper_netcdfcf'] = nansatMappers.pop('mapper_netcdfcf')

        # move generic_mapper to the end
        if 'mapper_generic' in nansatMappers:
            nansatMappers['mapper_generic'] = nansatMappers.pop('mapper_generic')

    return nansatMappers
