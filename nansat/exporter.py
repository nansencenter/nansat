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
from __future__ import print_function, absolute_import, division

import os
import tempfile
import datetime
import warnings

import gdal
import numpy as np

from netCDF4 import Dataset

from nansat.vrt import VRT
from nansat.node import Node

from nansat.warnings import NansatFutureWarning
from nansat.exceptions import NansatGDALError


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
    DEFAULT_INSTITUTE = 'NERSC'
    DEFAULT_SOURCE = 'satellite remote sensing'

    UNWANTED_METADATA = ['dataType', 'SourceFilename', 'SourceBand', '_Unsigned', 'FillValue',
                                'time', '_FillValue', 'type', 'scale', 'offset']

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
        add_geolocation : bool
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
        if add_geolocation:
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

        if not isinstance(bands, dict):
            raise ValueError('<bands> must be dict!')
        if metadata is None:
            metadata = {}
        # raise error if self is not projected (has GCPs)
        if len(self.vrt.dataset.GetGCPs()) > 0:
            raise ValueError('Cannot export dataset with GCPS for THREDDS!')

        # Create temporary empty Nansat object with self domain
        data = self.__class__.__new__(self.__class__)
        data._init_from_domain(self)

        # get mask (if exist)
        if mask_name is not None:
            mask = self[mask_name]

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
                raise NansatGDALError('%s is None' % str(iband))

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
            bandMetadata = self.get_metadata(band_id=iband)
            data.add_band(array=array, parameters=bandMetadata)
        self.logger.debug('Bands for export: %s' % str(dstBands))

        global_metadata = Exporter._set_global_metadata(created, data, metadata)

        # export temporary Nansat object to a temporary netCDF
        fid, tmp_filename = tempfile.mkstemp(suffix='.nc')
        data.export(tmp_filename, rm_metadata=rm_metadata)
        del data

        self._post_proc_thredds(tmp_filename, filename, bands, dstBands, time, global_metadata)

    def _create_dimensions(self, nc_inp, nc_out, time):
        """Create space and time dimenstions in the destination file"""
        # get time from Nansat object or from input datetime
        if time is None:
            time = self.time_coverage_start

        # collect info on dimension names
        dim_names = []
        grid_mapping_name = None
        for inp_var_name in nc_inp.variables:
            inp_var = nc_inp.variables[inp_var_name]
            dim_names += list(inp_var.dimensions)
            # get grid_mapping_name
            if hasattr(inp_var, 'grid_mapping_name'):
                grid_mapping_name = inp_var.grid_mapping_name
                grid_mapping_var_name = inp_var_name
        dim_names = list(set(dim_names))

        # collect info on dimension shapes
        dim_shapes = {}
        for dim_name in dim_names:
            dim_var = nc_inp.variables[dim_name]
            dim_shapes[dim_name] = dim_var.shape[0]

        # create dimensions
        for dim_name in dim_names:
            nc_out.createDimension(dim_name, dim_shapes[dim_name])

        # add time dimension
        nc_out.createDimension('time', 1)
        out_var = nc_out.createVariable('time', '>f8',  ('time', ))
        out_var.calendar = 'standard'
        out_var.long_name = 'time'
        out_var.standard_name = 'time'
        out_var.units = 'days since 1900-1-1 0:0:0 +0'
        out_var.axis = 'T'
        # create value for time variable
        td = time - datetime.datetime(1900, 1, 1)
        days = td.days + (float(td.seconds) / 60.0 / 60.0 / 24.0)
        # add date
        out_var[:] = days

        return grid_mapping_name, grid_mapping_var_name

    def _post_proc_thredds(self, tmp_filename, out_filename, bands, band_metadata, time, global_metadata):
        """Post processing of file for THREDDS (add time variable and metadata)"""
        # open files for input and output
        nc_inp = Dataset(tmp_filename, 'r')
        nc_out = Dataset(out_filename, 'w')

        grid_mapping_name, grid_mapping_var_name = self._create_dimensions(nc_inp, nc_out, time)

        # recreate file
        for inp_var_key in nc_inp.variables.keys():
            inp_var = nc_inp.variables[inp_var_key]
            if 'name' in inp_var.ncattrs():
                inp_var_name = inp_var.getncattr('name')
            else:
                inp_var_name = inp_var.name
            # create projection var
            if inp_var_name == grid_mapping_var_name:
                out_var = Exporter._copy_nc_var(inp_var, nc_out, grid_mapping_name,
                                                inp_var.dtype.str, inp_var.dimensions)
                continue

            # create simple x/y variables
            if inp_var_name in ['x', 'y', 'lon', 'lat']:
                out_var = Exporter._copy_nc_var(inp_var, nc_out, inp_var_name,
                                                '>f4', inp_var.dimensions)
            # create data var
            elif inp_var_name in band_metadata:
                fill_value = None
                if '_FillValue' in inp_var.ncattrs():
                    fill_value = inp_var._FillValue
                if '_FillValue' in band_metadata[inp_var_name]:
                    fill_value = band_metadata['_FillValue']
                dimensions = ('time', ) + inp_var.dimensions
                out_var = Exporter._copy_nc_var(inp_var, nc_out, inp_var_name,
                                                band_metadata[inp_var_name]['type'],
                                                dimensions, fill_value=fill_value)

            # copy array from input data
            data = inp_var[:]

            # copy rounded data from x/y
            if inp_var_name in ['x', 'y']:
                out_var[:] = np.floor(data).astype('>f4')
                # add axis=X or axis=Y
                out_var.axis = {'x': 'X', 'y': 'Y'}[inp_var_name]

            # copy data from lon/lat
            if inp_var_name in ['lon', 'lat']:
                out_var[:] = data.astype('>f4')

            # copy data from variables in the list
            if inp_var_name in band_metadata:
                # add offset and scale attributes
                scale = band_metadata[inp_var_name]['scale']
                offset = band_metadata[inp_var_name]['offset']
                if not (offset == 0.0 and scale == 1.0):
                    out_var.setncattr('add_offset', offset)
                    out_var.setncattr('scale_factor', scale)
                    data = (data - offset) / scale

                out_var[:] = data.astype(band_metadata[inp_var_name]['type'])

                # add custom attributes from input parameter bands
                if inp_var_name in bands:
                    for newAttr in bands[inp_var_name]:
                        if newAttr not in Exporter.UNWANTED_METADATA:
                            out_var.setncattr(newAttr, bands[inp_var_name][newAttr])
                    # add grid_mapping info
                    if grid_mapping_name is not None:
                        out_var.setncattr('grid_mapping', grid_mapping_name)

        # copy (some) global attributes
        for globAttr in nc_inp.ncattrs():
            if not(globAttr.strip().startswith('GDAL')):
                nc_out.setncattr(globAttr, nc_inp.getncattr(globAttr))

        # add common and custom global attributes
        nc_out.setncatts(global_metadata)

        # write output file
        nc_out.close()

        # close original files
        nc_inp.close()

        # Delete the temprary netCDF file
        fid = None
        os.remove(tmp_filename)

        return 0

    @staticmethod
    def _set_global_metadata(created, data, metadata):
        if created is None:
            created = (datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC'))
        # get corners of reprojected data
        min_lon, max_lon, min_lat, max_lat = data.get_min_max_lon_lat()

        global_metadata = {
                    'institution': Exporter.DEFAULT_INSTITUTE,
                    'source': Exporter.DEFAULT_SOURCE,
                    'creation_date': created,
                    'northernmost_latitude': np.float(max_lat),
                    'southernmost_latitude': np.float(min_lat),
                    'westernmost_longitude': np.float(min_lon),
                    'easternmost_longitude': np.float(max_lon),
                    'history': ' '}
        global_metadata.update(metadata)

        return global_metadata

    @staticmethod
    def _add_gcps(filename, gcps, bottomup):
        """Add 4 variables with gcps to the generated netCDF file"""
        gcp_variables = ['GCPX', 'GCPY', 'GCPZ', 'GCPPixel', 'GCPLine', ]

        # check if file exists
        if not os.path.exists(filename):
            warnings.warn('Cannot add GCPs! File %s doesn''t exist!' % filename)

        # open output file for adding GCPs
        ncFile = Dataset(filename, 'a')

        # get GCP values into single array from GCPs
        gcp_values = np.zeros((5, len(gcps)))
        for i, gcp in enumerate(gcps):
            gcp_values[0, i] = gcp.GCPX
            gcp_values[1, i] = gcp.GCPY
            gcp_values[2, i] = gcp.GCPZ
            gcp_values[3, i] = gcp.GCPPixel
            gcp_values[4, i] = gcp.GCPLine

        # make gcps dimensions
        ncFile.createDimension('gcps', len(gcps))
        # make gcps variables and add data
        for i, var in enumerate(gcp_variables):
            var = ncFile.createVariable(var, 'f4', ('gcps',))
            var[:] = gcp_values[i]

        # write data, close file
        ncFile.close()

    @staticmethod
    def _copy_nc_var(inp_var, nc_out, var_name, var_type, dimensions, fill_value=None):
        """Create new NC variable, set name, type, dimensions and copy attributes"""
        out_var = nc_out.createVariable(var_name, var_type, dimensions, fill_value=fill_value)
        for ncattr in inp_var.ncattrs():
            if str(ncattr) not in Exporter.UNWANTED_METADATA:
                out_var.setncattr(str(ncattr), inp_var.getncattr(ncattr))

        return out_var
