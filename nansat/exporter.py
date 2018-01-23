# Name:  exporter.py
# Purpose: Container of Exporter class
# Authors:      Anton Korosov, Artem Moiseev
# Created:      22.01.2018
# Copyright:    (c) NERSC 2011 - 2018
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
import tempfile
import datetime
import warnings

import gdal
import numpy as np

from netCDF4 import Dataset

from nansat.vrt import VRT
from nansat.tools import OptionError, GDALError, NansatFutureWarning
from nansat.node import Node


class Exporter(object):
    """Abstract class for export functions """
    EXPORT_FILENAME_WARNING = ('Nansat.export(fileName=...) will be disabled from Nansat 1.1. '
                               'Use Nansat.export(filename=...).')
    EXPORT_RM_METADATA_WARNING = ('Nansat.export(rmMetadata=...) will be disabled from Nansat 1.1. '
                                  'Use Nansat.export(rm_metadata=...).')
    EXPORT_ADD_GEOLOC_WARNING = ('Nansat.export(addGeoloc=...) will be disabled from Nansat 1.1. '
                                 'Use Nansat.export(add_geolocation=...).')
    EXPORT_ADD_GCPS_WARNING = ('Nansat.export(addGCPs=...) has no effect and '
                               'will be disabled from Nansat 1.1. ')
    EXPORT_BOTTOMUP_WARNING = ('Nansat.export(bottomup=...) will be disabled from Nansat 1.1. '
                               'Use Nansat.export(options=[WRITE_BOTTOMUP=NO]')
    EXPORT_CREATED_WARNING = ('Nansat.export2thredds(createdTime=...) will be disabled from Nansat 1.1. '
                               'Use Nansat.export2thredds(created=...')
    EXPORT_MASKNAME_WARNING = ('Nansat.export2thredds(maskName=...) will be disabled from Nansat 1.1. '
                               'Use Nansat.export2thredds(mask_name=...')

    def export(self, filename='', fileName='', bands=None, rm_metadata=None, rmMetadata=None,
               addGeoloc=None, add_geolocation=True,
               addGCPs=None, driver='netCDF', bottomup=None, options=None, hardcopy=False):
        '''Export Nansat object into netCDF or GTiff file

        Parameters
        -----------
        filename : str
            output file name
        bands: list (default=None)
            Specify band numbers to export.
            If None, all bands are exported.
        rm_metadata : list
            metadata names for removal before export.
            e.g. ['name', 'colormap', 'source', 'sourceBands']
        add_geoloc : bool
            add geolocation array datasets to exported file?
        driver : str
            Name of GDAL driver (format)
        bottomup : bool
            False: Write swath-projected data with rows and columns
                   organized as in the original product.
            True:  Use the default behaviour of GDAL, which is to flip the rows
        options : str or list
            GDAL export options in format of: 'OPT=VAL', or
            ['OPT1=VAL1', 'OP2='VAL2']
            See also http://www.gdal.org/frmt_netcdf.html
        hardcopy : bool
            Evaluate all bands just before export?

        Modifies
        ---------
        Create a netCDF file

        Notes
        ------
        If number of bands is more than one, serial numbers are added at the end of each band name.
        It is possible to fix it by changing line.4605 in GDAL/frmts/netcdf/netcdfdataset.cpp :
        'if( nBands > 1 ) sprintf(szBandName,"%s%d",tmpMetadata,iBand);'
        --> 'if( nBands > 1 ) sprintf(szBandName,"%s",tmpMetadata);'

        CreateCopy fails in case the band name has special characters,
        e.g. the slash in 'HH/VV'.

        Metadata strings with special characters are escaped with XML/HTML
        encoding.

        Examples
        --------
        # export all the bands into a netDCF 3 file
        >>> n.export(netcdfile)

        # export all bands into a GeoTiff
        >>> n.export(driver='GTiff')

        '''
        # trigger Nansat future warnings on deprecated features
        if fileName != '':
            warnings.warn(self.EXPORT_FILENAME_WARNING, NansatFutureWarning)
            filename = fileName
        if rmMetadata is not None:
            warnings.warn(self.EXPORT_RM_METADATA_WARNING, NansatFutureWarning)
            rm_metadata = rmMetadata
        if addGeoloc is not None:
            warnings.warn(self.EXPORT_ADD_GEOLOC_WARNING, NansatFutureWarning)
            add_geolocation = addGeoloc
        if addGCPs is not None:
            warnings.warn(self.EXPORT_ADD_GCPS_WARNING, NansatFutureWarning)
        if isinstance(bottomup, bool):
            warnings.warn(self.EXPORT_BOTTOMUP_WARNING, NansatFutureWarning)

        if options is None:
            options = []
        if type(options) == str:
            options = [options]

        # temporary VRT for exporting
        export_vrt = self.vrt.copy()
        export_vrt.leave_few_bands(bands)
        export_vrt.split_complex_bands()
        if addGeoloc:
            export_vrt.create_geolocation_bands()
        export_vrt.fix_band_metadata(rm_metadata)
        export_vrt.fix_global_metadata(rm_metadata)

        # if output filename is the same as input one
        if self.filename == filename or hardcopy:
            export_vrt.hardcopy_bands()

        if driver == 'GTiff':
            options, add_gcps = export_vrt.prepare_export_gtiff(options)
        else:
            options, add_gcps = export_vrt.prepare_export_netcdf(options, bottomup)

        # Create output file using GDAL
        dataset = gdal.GetDriverByName(driver).CreateCopy(filename, export_vrt.dataset, options=options)
        dataset = None
        # add GCPs into netCDF file as separate float variables
        if add_gcps:
            Exporter._add_gcps(filename, export_vrt.dataset.GetGCPs(), bottomup)

        self.logger.debug('Export - OK!')

    @staticmethod
    def _add_gcps(filename, gcps, bottomup):
        ''' Add 4 variables with gcps to the generated netCDF file '''
        gcpVariables = ['GCPX', 'GCPY', 'GCPZ', 'GCPPixel', 'GCPLine', ]

        # check if file exists
        if not os.path.exists(filename):
            warnings.warn('Cannot add GCPs! File %s doesn''t exist!' % filename)

        # open output file for adding GCPs
        ncFile = Dataset(filename, 'a')

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

    def export2thredds(self, filename, bands, metadata=None,
                        mask_name=None, maskName=None, rm_metadata=None, rmMetadata=None,
                        time=None, created=None, createdTime=None):
        ''' Export data into a netCDF formatted for THREDDS server

        Parameters
        ----------
        filename : str
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
        mask_name: str
            if data include a mask band: give the mask name.
            Non-masked value is 64.
            if None: no mask is added
        rm_metadata : list
            unwanted metadata names which will be removed
        time : list with datetime objects
            aqcuisition time of original data. That value will be in time dim
        created : datetime
            date of creation. Will be in metadata 'created'

        Note
        ----
        Nansat object (self) has to be projected (with valid GeoTransform and
        valid Spatial reference information) but not wth GCPs

        Examples
        --------
        # create THREDDS formatted netcdf file with all bands and time variable
        >>> n.export2thredds(filename)

        # export only the first band and add global metadata
        >>> n.export2thredds(filename, ['L_469'], {'description': 'example'})

        # export several bands and modify type, scale and offset
        >>> bands = {'L_645' : {'type': '>i2', 'scale': 0.1, 'offset': 0},
                     'L_555' : {'type': '>i2', 'scale': 0.1, 'offset': 0}}
        >>> n.export2thredds(filename, bands)

        '''
        if rmMetadata is not None:
            warnings.warn(self.EXPORT_RM_METADATA_WARNING, NansatFutureWarning)
            rm_metadata = rmMetadata
        if createdTime is not None:
            warnings.warn(self.EXPORT_CREATED_WARNING, NansatFutureWarning)
            created = createdTime
        if maskName is not None:
            warnings.warn(self.EXPORT_MASKNAME_WARNING, NansatFutureWarning)
            mask_name = maskName

        # TODO:
        # New logics
        # 1. Get mask
        # 1. Prepare global metadata
        # 1. self.export to netCDF with metadata, hardcopy and mask (add mask to vrt.harcopy)
        # 1. self._collect_band_metadata_from_exported
        # 1. self._post_proc_for_thredds

        # raise error if self is not projected (has GCPs)
        if len(self.vrt.dataset.GetGCPs()) > 0:
            raise OptionError('Cannot export dataset with GCPS for THREDDS!')

        # Create temporary empty Nansat object with self domain
        data = self.__class__.__new__(self.__class__)
        data._init_from_domain(self)

        # get mask (if exist)
        if maskName is not None:
            mask = self[maskName]

# TODO: move to Exporter._hardcopy_bands
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
        minLon, maxLon, minLat, maxLat = data.get_min_max_lon_lat()

# TODO: move to Exporter._set_global_metadata
        # common global attributes:
        if createdTime is None:
            createdTime = (datetime.datetime.utcnow().
                           strftime('%Y-%m-%d %H:%M:%S UTC'))
# TODO: move 'NERSC', etc... to constants
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
        fid, tmp_filename = tempfile.mkstemp(suffix='.nc')
        data.export(tmp_filename)

        self._post_proc_thredds(tmp_filename, filename, dstBands)

    def _post_proc_thredds(self, tmp_filename, out_filename, band_metadata):
        """Post processing of file for THREDDS (add time variable and metadata)"""
        # open files for input and output
        ncI = Dataset(tmp_filename, 'r')
        ncO = Dataset(out_filename, 'w')

# TODO: move to Exporter._ctreate_time_dim
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

# TODO: move to constantd in Exporter
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

# TODO: move to Exporter._add_thredds_bands
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
                if ncIVar_name in ['lon', 'lat']:
                    dimensions = ncIVar.dimensions
                else:
                    dimensions = ('time', ) + ncIVar.dimensions

                fill_value = None
                if '_FillValue' in ncIVar.ncattrs():
                    fill_value = ncIVar._FillValue
                if '_FillValue' in dstBands[ncIVar_name]:
                    fill_value = dstBands['_FillValue']
                ncOVar = ncO.createVariable(ncIVar_name,
                                            dstBands[ncIVar_name]['type'],
                                            dimensions, fill_value=fill_value)

            # copy array from input data
            data = ncIVar[:]

            for ncattr in ncIVar.ncattrs():
                if ncattr == '_FillValue':
                    continue
                ncOVar.setncattr(ncattr, ncIVar.getncattr(ncattr))

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

                ncOVar[:] = data.astype(dstBands[ncIVar_name]['type'])
                # copy (some) attributes
                for inAttrName in ncIVar.ncattrs():
                    if str(inAttrName) not in rmMetadata + ['dataType',
                                    'SourceFilename', 'SourceBand', '_Unsigned',
                                    'FillValue', 'time', '_FillValue']:
                        ncOVar.setncattr(inAttrName, ncIVar.getncattr(inAttrName))

                # add custom attributes from input parameter bands
                if ncIVar_name in bands:
                    for newAttr in bands[ncIVar_name]:
                        if newAttr not in rmMetadata + ['type', 'scale',
                                                        'offset',
                                                        '_FillValue']:
                            ncOVar.setncattr(newAttr, bands[ncIVar_name][newAttr])
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
        fid = None
        os.remove(tmpName)

        return 0

