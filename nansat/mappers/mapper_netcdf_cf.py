''' Nansat NetCDF-CF mapper

    Check CF-compliance of your files here:
    http://cfconventions.org/compliance-checker.html
'''

import collections
import datetime
import os
import warnings

import numpy as np

from dateutil.parser import parse
from netCDF4 import Dataset
from osgeo import gdal
from pyproj import CRS

from nansat.vrt import VRT
from nansat.nsr import NSR
from nansat.utils import parse_time

from nansat.exceptions import WrongMapperError, NansatMissingProjectionError

class ContinueI(Exception):
    pass

# List of allowed spatial dimensions
ALLOWED_SPATIAL_DIMENSIONS_X = ['x', 'lon', 'numcells']
ALLOWED_SPATIAL_DIMENSIONS_Y = ['y', 'lat', 'numrows']

class Mapper(VRT):
    """
    """

    def __init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs):

        if not filename.endswith('nc'):
            raise WrongMapperError

        self.input_filename = filename

        if not gdal_metadata:
            raise WrongMapperError

        if 'NC_GLOBAL#GDAL_NANSAT_GCPY_000' in list(gdal_metadata.keys()) or \
                'NC_GLOBAL#GDAL_NANSAT_GCPProjection' in list(gdal_metadata.keys()) or \
                'NC_GLOBAL#NANSAT_GCPY_000' in list(gdal_metadata.keys()) or \
                'NC_GLOBAL#NANSAT_GCPProjection' in list(gdal_metadata.keys()):
            # Nansat generated netcdf of swath data is not standard CF
            raise WrongMapperError

        metadata = VRT._remove_strings_in_metadata_keys(gdal_metadata,
                                                        ['NC_GLOBAL#', 'NANSAT_', 'GDAL_'])

        # Set origin metadata (TODO: agree on keyword...)
        origin = ''
        nans = 'NANSAT'
        if 'origin' in list(metadata.keys()):
            origin = metadata['origin'] + ' '
        for key in list(metadata.keys()):
            if nans in key:
                metadata['origin'] = origin + nans
            # else: Nothing needs to be done, origin stays the same...

        # Check conventions metadata
        if 'Conventions' not in list(metadata.keys()) or 'CF' not in metadata['Conventions']:
            raise WrongMapperError

        # OBS: at this point, generic mapper fails...
        #if metadata.has_key('GCPProjection'):
        #    # Probably Nansat generated netcdf of swath data - see issue #192
        #    raise WrongMapperError

        # Create empty VRT dataset with geo-reference
        self._create_empty(gdal_dataset, metadata)

        # Add bands with metadata and corresponding values to the empty VRT
        self.create_bands(self._band_list(gdal_dataset, metadata, *args, **kwargs))

        # Check size?
        #xsize, ysize = self.ds_size(sub0)

        # Create complex bands from *_real and *_imag bands (the function is in
        # vrt.py)
        self._create_complex_bands(self._get_sub_filenames(gdal_dataset))

        # Set GCMD/DIF compatible metadata if available
        self._set_time_coverage_metadata(metadata)

        # Add global metadata
        metadata.pop("title_no", "")
        metadata.pop("summary_no", "")
        self.dataset.SetMetadata(metadata)

        # Then add remaining GCMD/DIF compatible metadata in inheriting mappers

    def times(self):
        ''' Get times from time variable

        NOTE: This cannot be done with gdal because the time variable is a
        vector

        '''
        ds = Dataset(self.input_filename)

        # Get datetime object of epoch and time_units string
        time_units = self._time_reference(ds=ds)

        # Get all times - consider caching to save time (see
        # nansat/mappers/opendap.py)
        times = ds.variables[self._timevarname(ds=ds)]

        # Create numpy array of np.datetime64 times (provide epoch to save time)
        tt = np.array([self._time_count_to_np_datetime64(tn,
            time_reference=time_units) for tn in times])

        return tt

    def _time_reference(self, ds=None):
        """ Get the time reference of the dataset

        Returns:
        --------
            (epoch, units) : tuple
                Datetime.datetime object of the time epoch
                String representation of time units and epoch
        """
        if not ds:
            ds = Dataset(self.input_filename)
        times = ds.variables[self._timevarname(ds=ds)]
        rt = parse(times.units, fuzzy=True) # This sets timezone to local
        # Remove timezone information from epoch, which defaults to
        # utc (otherwise the timezone should be given in the dataset)
        units = times.units
        epoch = datetime.datetime(rt.year, rt.month, rt.day, rt.hour,
                rt.minute, rt.second)
        return epoch, units

    def _timevarname(self, ds=None):
        if not ds:
            ds = Dataset(self.input_filename)
        timevarname = ''
        std_name = 'time'
        if not std_name in ds.variables.keys():
            for key, var in ds.variables.items():
                if 'standard_name' in var.ncattrs() and var.standard_name == std_name:
                    timevarname = key
        else:
            timevarname = std_name
        return timevarname

    def _time_count_to_np_datetime64(self, time_count, time_reference=None):
        if not time_reference:
            time_reference = self._time_reference()
        time_count = float(time_count)
        time_decimal = time_count - np.floor(time_count)
        if 'second' in time_reference[1]:
            tt = np.datetime64(time_reference[0] +
                    datetime.timedelta(seconds=time_count))
        elif 'minute' in time_reference[1]:
            minutes = np.floor(time_count)
            seconds = time_decimal*60
            tt = np.datetime64(time_reference[0] + datetime.timedelta(minutes=minutes,
                seconds=seconds))
        elif 'hour' in time_reference[1]:
            hours = np.floor(time_count)
            minutes = np.floor(time_decimal*60)
            seconds = (time_decimal*60 - minutes)*60
            tt = np.datetime64(time_reference[0] + datetime.timedelta(hours=hours, minutes=minutes,
                seconds=seconds))
        elif 'day' in time_reference[1]:
            days = np.floor(time_count)
            hours = np.floor(time_decimal*24)
            minutes = np.floor((time_decimal*24 - hours)*60)
            seconds = ((time_decimal*24 - hours)*60 - minutes)*60
            tt = np.datetime64(time_reference[0] + datetime.timedelta(days=days, hours=hours,
                minutes=minutes, seconds=seconds))
        else:
            raise Exception('Check time units..')
        return tt

    def _band_list(self, gdal_dataset, gdal_metadata, netcdf_dim=None, bands=None, *args,
                   **kwargs):
        ''' Create list of dictionaries mapping source and destination
        metadata of bands that should be added to the Nansat object.

        Parameters
        ----------
        gdal_dataset : gdal.Dataset
            The gdal dataset opened in nansat.py
        gdal_metadata : dict
            Dictionary of global metadata
        netcdf_dim : dict
            Dictionary of desired slice of a multi-dimensional array.
            Since gdal only returns 2D bands, a multi-dimensional
            array (x,y,z) is split into z bands accompanied with
            metadata information about the position of the slice along
            the z-axis. The (key, value) pairs represent dimension
            name and the desired value in that dimension (not the
            index), respectively. E.g., for a height dimension of
            [10, 20, 30] m, netcdf_dim = {'height': 20} if you want to
            extract the data from 20 m height.
        bands : list
            List of desired bands following NetCDF-CF standard names.
            NOTE: some datasets have other bands as well, i.e., of
            data not yet implemented in CF. We may at some point
            generalize this to provide a dict with key name and value,
            where the key is, e.g., "standard_name" or "metno_name",
            etc.
        '''
        if netcdf_dim is None:
            netcdf_dim = {}
        if bands is None:
            ds = Dataset(self.input_filename)
            bands = ds.variables.keys()
        metadictlist = []
        for fn in self._get_sub_filenames(gdal_dataset):
            band = fn.split(':')[-1]
            if band not in bands:
                continue
            # Don't process geolocation variables/bands
            if ('GEOLOCATION_X_DATASET' in fn or 'longitude' in fn or
                    'GEOLOCATION_Y_DATASET' in fn or 'latitude' in fn):
                continue
            try:
                metadictlist.append(
                    self._get_band_metadata(fn, netcdf_dim=netcdf_dim))
            except ContinueI:
                continue
        return metadictlist

    @staticmethod
    def _pop_spatial_dimensions(dimension_names):
        """ Pop spatial dimensions (longitude and latitude, or x and y)
        """
        for allowed in ALLOWED_SPATIAL_DIMENSIONS_X:
            try:
                ind_dim_x = [i for i, s in enumerate(dimension_names) if allowed in s.lower()][0]
                dimension_names.pop(ind_dim_x)
            except IndexError:
                continue
        for allowed in ALLOWED_SPATIAL_DIMENSIONS_Y:
            try:
                ind_dim_y = [i for i, s in enumerate(dimension_names) if allowed in s.lower()][0]
                dimension_names.pop(ind_dim_y)
            except IndexError:
                continue

    def _get_index_of_dimensions(self, dimension_names, netcdf_dim, dim_sizes):
        """ Nansat only works with 2D data. For for datasets with
        higher dimensions, we need to pick the correct slice. The
        function returns a dictionary with the index of the wanted
        data slice, as specified in the netcdf_dim dictionary.
        """
        ds = Dataset(self.input_filename)
        index4key = collections.OrderedDict()
        for key in dimension_names:
            if key in netcdf_dim.keys():
                val = netcdf_dim[key]
                if key == 'time':
                    if type(val) != datetime.datetime and type(val) != np.datetime64:
                        raise ValueError
                    if type(val) == datetime.datetime:
                        val = np.datetime64(val).astype('M8[s]')
                    # Get band number from given timestamp
                    index = int(np.argmin(np.abs(self.times() - val)))
                else:
                    index = int(np.argmin(np.abs(ds.variables[key][:] - val)))
                index4key[key] = {
                        'index': index,
                        'size': dim_sizes[key],
                    }
            else:
                index4key[key] = {
                        'index': 0,
                        'size': dim_sizes[key],
                    }
        return index4key

    @staticmethod
    def _calculate_band_number(index4dimension, dimension_names):
        """ Uses index and size values in provided dict to calculate
        gdal band number. The index4dimension dict contains the index
        in the dimensions given by dimension_names, and the size of
        the dimensions.
        """
        class Context:
            band_number = 1
            multiplier = 1

        def get_band_number():
            try:
                name_dim0 = dimension_names.pop(0)
            except:
                return
            Context.band_number += index4dimension[name_dim0]['index']*Context.multiplier
            Context.multiplier *= index4dimension[name_dim0]['size']
            get_band_number()

        get_band_number()

        return Context.band_number

    def _get_dimension_info(self, band):
        """ Get band list from a netCDF4.Dataset. See docs of
        _band_list.
        """
        ds = Dataset(self.input_filename)
        var = ds.variables[band]
        dimension_names = [b.name for b in var.get_dims()]
        dimension_names.reverse()
        dim_sizes = collections.OrderedDict()
        for dim in var.get_dims():
            dim_sizes[dim.name] = dim.size

        return dimension_names, dim_sizes

    def _get_band_number(self, band, netcdf_dim={}):
        """ Get gdal band number from multidimensional dataset.
        """
        if netcdf_dim is None:
            netcdf_dim = {}
        
        dimension_names, dim_sizes = self._get_dimension_info(band)
        self._pop_spatial_dimensions(dimension_names)
        index = self._get_index_of_dimensions(dimension_names, netcdf_dim, dim_sizes)
        band_number = self._calculate_band_number(index, dimension_names)

        return band_number

    def _get_band_metadata(self, fn, netcdf_dim=None):
        """ Gets band metadata using gdal
        """
        if netcdf_dim is None:
            netcdf_dim = {}
        band = fn.split(':')[-1]
        band_number = self._get_band_number(band, netcdf_dim)

        subds = gdal.Open(fn)
        band = subds.GetRasterBand(band_number)
        band_metadata = self._clean_band_metadata(band)
        band_metadata['_FillValue'] = band.GetNoDataValue()

        return self._band_dict(fn, band_number, subds, band=band,
                        band_metadata=band_metadata)

    def _clean_band_metadata(self, band, remove = ['_Unsigned', 'ScaleRatio',
        'ScaleOffset', 'PixelFunctionType']):

        band_metadata = band.GetMetadata_Dict()
        for key in remove:
            if key in band_metadata:
                band_metadata.pop(key)

        return band_metadata

    def _band_dict(self, subfilename, band_num, subds, band=None, band_metadata=None):
        '''
        subfilename : string
            Name of subdataset file
        '''

        if not band:
            try:
                band = subds.GetRasterBand(band_num)
            except RuntimeError as e:
                if 'illegal band' in str(e).lower():
                    warnings.warn('Skipping band due to GDAL error: %s' %str(e))
                    return {}
                else:
                    raise
        if not band_metadata:
            band_metadata = self._clean_band_metadata(band)

        if 'time_iso_8601' not in list(band_metadata.keys()):
            if self._timevarname() in list(band_metadata.keys()):
                timecountname = self._timevarname()
            else:
                timecountname = 'NETCDF_DIM_'+self._timevarname()
            try:
                band_metadata['time_iso_8601'] = self._time_count_to_np_datetime64(
                    band_metadata[timecountname])
            except (ValueError, KeyError) as e:
                # No timing information available for this band - it is
                # probably a constant, such as land area fraction or similar.
                # Then we don't need time for this band...
                # The following warning is not really understandable..
                tttt = 0 # do nothing...
                #warnings.warn(
                #        '%s: %s - %s Continuing without time metadata for band %s'
                #        %(e.__repr__().split('(')[0], str(e), e.__doc__,
                #            band_metadata['NETCDF_VARNAME']))

        # Generate source metadata
        src = {'SourceFilename': subfilename, 'SourceBand': band_num}

        # Set scale ratio
        scaleRatio = band_metadata.get('ScaleRatio',
                band_metadata.get('scale',
                    band_metadata.get('scale_factor', '')))
        if len(scaleRatio) > 0:
            src['ScaleRatio'] = scaleRatio

        # Set scale offset
        scaleOffset = band_metadata.get('ScaleOffset',
                band_metadata.get('offset',
                    band_metadata.get('add_offset', '')))
        if len(scaleOffset) > 0:
            src['ScaleOffset'] = scaleOffset

        # Set data type
        src['DataType'] = band.DataType

        # Generate destination metadata
        # Copy all metadata from input band
        dst = band_metadata
        # Set wkv
        dst['wkv'] = band_metadata.get('standard_name', '')

        # Set band name
        if 'name' in band_metadata:
            bandName = band_metadata['name']
        else:
            # if it doesn't exist get name from NETCDF_VARNAME
            bandName = band_metadata.get('NETCDF_VARNAME', '')
            if len(bandName) == 0:
                bandName = band_metadata.get('dods_variable', '')

            # remove digits added by gdal when exporting to netcdf
            if (len(bandName) > 0 and 'origin' in band_metadata.keys()
                    and 'nansat' in band_metadata['origin'].lower()):
                if bandName[-1:].isdigit():
                    bandName = bandName[:-1]
                if bandName[-1:].isdigit():
                    bandName = bandName[:-1]

        # if still no bandname, create one
        if len(bandName) == 0:
            bandName = 'band_%03d' % band_num

        dst['name'] = bandName

        return {'src': src, 'dst': dst}

    def _create_empty(self, gdal_dataset, gdal_metadata):
        try:
            self._create_empty_from_projection_variable(gdal_dataset, gdal_metadata)
        except (KeyError, AttributeError):
            try:
                self._create_empty_from_subdatasets(gdal_dataset, gdal_metadata)
            except NansatMissingProjectionError:
                # this happens rarely - a special workaround is done in mapper_quikscat
                warnings.warn('GDAL cannot determine the dataset projection - using Nansat ' \
                        'spatial reference WKT, assuming a regular longitude/latitude grid')
                self._create_empty_with_nansat_spatial_reference_WKT(gdal_dataset, gdal_metadata)

    def _create_empty_from_projection_variable(self, gdal_dataset, gdal_metadata):

        # TODO: only open dataset once or close them..
        ds = Dataset(self.input_filename)

        shape=[]
        for key in ds.dimensions.keys():
            shape.append(ds.dimensions[key].size)

        for var in ds.variables.keys():
            if ds[var].shape==tuple(shape):
                break

        grid_mapping = ds.variables[ds.variables[var].grid_mapping]

        grid_mapping_dict = {}
        for index in grid_mapping.ncattrs():
            grid_mapping_dict[index] = grid_mapping.getncattr(index)
        crs = CRS.from_cf(grid_mapping_dict)

        try:
            subdataset = gdal.Open(self._get_sub_filenames(gdal_dataset)[0])
        except IndexError:
            subdataset = gdal_dataset
        self._init_from_dataset_params(
                    x_size = subdataset.RasterXSize,
                    y_size = subdataset.RasterYSize,
                    geo_transform = subdataset.GetGeoTransform(),
                    projection = NSR(crs.to_proj4()).wkt,
                    metadata = gdal_metadata)

    def _create_empty_from_subdatasets(self, gdal_dataset, metadata):
        """ Create empty vrt dataset with projection but no bands
        """
        no_projection = True
        # Check if the main gdal dataset has projection
        if gdal_dataset.GetProjection():
            sub = gdal_dataset
            no_projection = False
        else:
            # Loop subdatasets and check for projection
            for fn in self._get_sub_filenames(gdal_dataset):
                sub = gdal.Open(fn)
                if sub.GetProjection():
                    no_projection = False
                    break
        if no_projection:
            #raise NansatMissingProjectionError
            sub_datasets = gdal_dataset.GetSubDatasets()
            filenames = [f[0] for f in sub_datasets]
            for _, filename in enumerate(filenames):
                if 'longitude' in filename:
                    xDatasetSource = filename
                if 'latitude' in filename:
                    yDatasetSource = filename
            if not "xDatasetSource" in locals():
                raise NansatMissingProjectionError
            if not "yDatasetSource" in locals():
                raise NansatMissingProjectionError
            lon_dataset = gdal.Open(xDatasetSource)
            lon = lon_dataset.GetRasterBand(1).ReadAsArray()
            lat_dataset = gdal.Open(yDatasetSource)
            lat = lat_dataset.GetRasterBand(1).ReadAsArray()
            self._init_from_lonlat(lon, lat, add_gcps=False)
        else:
            # Initialise VRT with subdataset containing projection
            self._init_from_gdal_dataset(sub, metadata=metadata)

    def _create_empty_with_nansat_spatial_reference_WKT(self, gdal_dataset, gdal_metadata):
        """ In this case, gdal cannot find the projection of any
            subdatasets. We therefore assume a regular longitude/latitude grid,
            and set the projection to the Nansat Spatial Reference WKT
            [NSR().wkt], using the first subdataset as source
        """
        # NB: if gdal_dataset has subdatasets, gdal_dataset.RasterXSize equals 512,
        # gdal_dataset.RasterYSize equals 512, gdal_dataset.GetGeoTransform returns
        # (1,0,0,0,1,0). This is wrong.
        fn = self._get_sub_filenames(gdal_dataset)
        if  len(fn) == 0:
            subdataset = gdal_dataset
        else:
            subdataset = gdal.Open(fn[0])

        self._init_from_dataset_params(
                    x_size = subdataset.RasterXSize,
                    y_size = subdataset.RasterYSize,
                    geo_transform = subdataset.GetGeoTransform(),
                    projection = NSR().wkt,
                    metadata = gdal_metadata)

    def _set_time_coverage_metadata(self, gdal_metadata):
        ### GET START TIME from METADATA
        time_coverage_start = None
        if 'time_coverage_start' in gdal_metadata:
            time_coverage_start = parse_time(
                                    gdal_metadata['time_coverage_start'])

        ### GET END TIME from METADATA
        time_coverage_end = None
        if 'time_coverage_end' in gdal_metadata:
            time_coverage_end = parse_time(
                                    gdal_metadata['time_coverage_end'])

        # set time_coverage_start if available
        if time_coverage_start is not None:
            self.dataset.SetMetadataItem('time_coverage_start',
                                    time_coverage_start.isoformat())
        # set time_coverage_end if available
        if time_coverage_end is not None:
            self.dataset.SetMetadataItem('time_coverage_end',
                                    time_coverage_end.isoformat())

