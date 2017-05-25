''' Nansat NetCDF-CF mapper

    Check CF-compliance of your files here:
    http://cfconventions.org/compliance-checker.html
'''
import warnings, os, datetime
import numpy as np
import gdal

from cfunits import Units
from scipy.io.netcdf import netcdf_file
from dateutil.parser import parse
from netCDF4 import Dataset

from nansat.vrt import VRT, GeolocationArray
from nansat.nsr import NSR
from nansat.tools import WrongMapperError, parse_time

class Mapper(VRT):

    def __init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs):

        # test_nansat is failing - this mapper needs more work..
        #raise WrongMapperError

        if not filename.endswith('nc'):
            raise WrongMapperError

        if not gdal_metadata:
            raise WrongMapperError

        if gdal_metadata.has_key('NC_GLOBAL#GDAL_NANSAT_GCPY_000') or \
                gdal_metadata.has_key('NC_GLOBAL#GDAL_NANSAT_GCPProjection'):
            # Probably Nansat generated netcdf of swath data - see issue #192
            raise WrongMapperError
        
        metadata = self._remove_strings_in_metadata_keys(gdal_metadata)

        # Set origin metadata (TODO: agree on keyword...)
        origin = ''
        nans = 'NANSAT'
        if metadata.has_key('origin'):
            origin = metadata['origin'] + ' '
        for key in metadata.keys():
            if nans in key:
                metadata['origin'] =  origin + nans
            # else: Nothing needs to be done, origin stays the same...

        # Check conventions metadata
        if not metadata.has_key('Conventions') or not 'CF' in metadata['Conventions']:
            raise WrongMapperError

        # OBS: at this point, generic mapper fails...
        #if metadata.has_key('GCPProjection'):
        #    # Probably Nansat generated netcdf of swath data - see issue #192
        #    raise WrongMapperError

        # Create empty VRT dataset with geo-reference
        self._create_empty(gdal_dataset, metadata)

        # Add bands with metadata and corresponding values to the empty VRT
        self._create_bands(self._band_list(filename, metadata, *args, **kwargs))

        # Check size?
        #xsize, ysize = self.ds_size(sub0)

        # Create complex bands from *_real and *_imag bands (the function is in
        # vrt.py)
        self._create_complex_bands(self.sub_filenames(gdal_dataset))

        # Set GCMD/DIF compatible metadata if available
        self._set_time_coverage_metadata(metadata)

        # Then add remaining GCMD/DIF compatible metadata in inheriting mappers

    def times(self, filename):
        ''' Get times from time variable 

        NOTE: This cannot be done with gdal because the time variable is a
        vector

        Parameters
        ----------
        filename : string
            Name of the original CF-compliant NetCDF file
        '''
        ds = Dataset(filename)
        # Consider caching to save time - see nansat/mappers/opendap.py

        times = ds.variables['time']
        rt = parse(times.units, fuzzy=True) # Sets timezone to local
        # Remove timezone information from reference time, which defaults to
        # utc (otherwise the timezone should be given in the dataset)
        t0 = datetime.datetime(rt.year, rt.month, rt.day, rt.hour,
                rt.minute, rt.second)
        # Create numpy array of np.datetime64 times
        if 'second' in times.units:
            tt = np.array([np.datetime64(t0 + datetime.timedelta(seconds=tn)) for
                tn in times])
        elif 'hour' in times.units:
            tt = np.array([np.datetime64(t0 + datetime.timedelta(hours=int(tn))) for
                tn in times])
        else:
            raise Exception('Check time units..')

        return tt

    def _band_list(self, filename, gdal_metadata, netcdf_dim={}, bands=[]):
        ''' Create list of dictionaries mapping source and destination metadata
        of bands that should be added to the Nansat object.

        Parameters
        ----------
        filename : string
            Name of the original CF-compliant NetCDF file
        gdal_metadata : dict
            Dictionary of global metadata
        netcdf_dim : dict
            Dictionary of desired slice of a multi-dimensional array. Since
            gdal only returns 2D bands, a multi-dimensional array (x,y,z) is
            split into z bands accompanied with metadata information about the
            position of the slice along the z-axis.
        bands : list
            List of desired bands following NetCDF-CF standard names. NOTE:
            some datasets have other bands as well, i.e., of data not yet
            implemented in CF. We may at some point generalize this to provide
            a dict with key name and value, where the key is, e.g.,
            "standard_name" or "metno_name", etc.
        '''
        class ContinueI(Exception):
            pass
        class BreakI(Exception):
            pass
        metadictlist = []
        gdal_dataset = gdal.Open(filename)
        for fn in self.sub_filenames(gdal_dataset):
            if ('GEOLOCATION_X_DATASET' in fn or 'longitude' in fn or
                    'GEOLOCATION_Y_DATASET' in fn or 'latitude' in fn):
                continue
            subds = gdal.Open(fn)
            for i in range(subds.RasterCount):
                band_num = i + 1
                band = subds.GetRasterBand(band_num)
                band_metadata = self._clean_band_metadata(band)
                # Keep only desired bands (given in "bands" list)
                try:
                    if bands:
                        if not band_metadata.has_key('standard_name'):
                            raise ContinueI
                        if not band_metadata['standard_name'] in bands:
                            raise ContinueI
                except ContinueI:
                    continue
                # Keep only desired slices following "netcdf_dim" dictionary
                try:
                    for key, val in netcdf_dim.iteritems():
                        match = [s for s in band_metadata if key in s]
                        if key=='time' and type(val)==np.datetime64:
                            # Select band directly from given timestamp, and
                            # break the for loop
                            band_num = np.argmin(np.abs(self.times(filename) -
                                val))
                            metadictlist.append(self._band_dict(fn,
                                band_num, subds))
                            raise BreakI
                        if not match or not band_metadata[match[0]]==val:
                            raise ContinueI
                        else:
                            band_metadata[key] = band_metadata.pop(match[0])
                except ContinueI:
                    continue
                except BreakI:
                    break

                # append band with src and dst dictionaries
                metadictlist.append(self._band_dict(fn, band_num, subds,
                    band=band, band_metadata=band_metadata))

        return metadictlist

    def _clean_band_metadata(self, band, remove = ['_Unsigned', 'ScaleRatio',
        'ScaleOffset', 'PixelFunctionType']):

        band_metadata = band.GetMetadata_Dict()
        for key in remove:
            if key in band_metadata:
                band_metadata.pop(key)

        return band_metadata

    def _band_dict(self, filename, band_num, subds, band=None, band_metadata=None):

        if not band:
            band = subds.GetRasterBand(band_num)
        if not band_metadata:
            band_metadata = self._clean_band_metadata(band)

        # Generate source metadata
        src = {'SourceFilename': filename, 'SourceBand': band_num}

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
            if (len(bandName) > 0 and band_metadata.has_key('origin')
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
        subfiles = self.sub_filenames_with_projection(gdal_dataset)
        if not subfiles:
            ''' In this case, gdal cannot find the projection of any
            subdatasets. We therefore assume a regular longitude/latitude grid,
            and set the projection to the Nansat Spatial Reference WKT
            [NSR().wkt], using the first subdataset as source
            '''
            fn = self.sub_filenames(gdal_dataset)
            if not fn:
                raise WrongMapperError
            sub = gdal.Open(fn[0])
            super(Mapper, self).__init__(
                    srcRasterXSize = sub.RasterXSize,
                    srcRasterYSize = sub.RasterYSize, 
                    srcGeoTransform = sub.GetGeoTransform(), 
                    srcProjection = NSR().wkt, 
                    srcMetadata = gdal_metadata)
        else:
            sub0 = gdal.Open(subfiles[0])
            super(Mapper, self).__init__(gdalDataset=sub0, srcMetadata=gdal_metadata)

    def _remove_strings_in_metadata_keys(self, gdal_metadata):
        if not gdal_metadata:
            raise WrongMapperError

        meta = gdal_metadata.copy()

        # These strings are added when datasets are exported in
        # nansat.nansat.Nansat.export...
        rm_strings = ['NC_GLOBAL#', 'NANSAT_', 'GDAL_']

        for rms in rm_strings:
            for key in meta.keys():
                newkey = key.replace(rms, '')
                meta[newkey] = meta.pop(key)

        return meta

    def sub_filenames(self, gdal_dataset):
        # Get filenames of subdatasets
        sub_datasets = gdal_dataset.GetSubDatasets()
        return [f[0] for f in sub_datasets]

    def sub_filenames_with_projection(self, gdal_dataset):
        # Get filenames of subdatasets containing projection
        sub_fnames = self.sub_filenames(gdal_dataset)
        with_proj = []
        for f in sub_fnames:
            if gdal.Open(f).GetProjection():
                with_proj.append(f)
        return with_proj

    #def ds_size(self, ds):
    #    return ds.RasterXSize, ds.RasterYSize

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

        ### GET start time from time variable
        if (time_coverage_start is None and
                 'time#standard_name' in gdal_metadata and
                 gdal_metadata['time#standard_name'] == 'time' and
                 'time#units' in gdal_metadata and
                 'time#calendar' in gdal_metadata):
            # get data from netcdf data
            ncFile = netcdf_file(inputFileName, 'r')
            timeLength = ncFile.variables['time'].shape[0]
            timeValueStart = ncFile.variables['time'][0]
            timeValueEnd = ncFile.variables['time'][-1]
            ncFile.close()
            try:
                timeDeltaStart = Units.conform(timeValueStart,
                                  Units(subMetadata['time#units'],
                                        calendar=subMetadata['time#calendar']),
                                  Units('days since 1950-01-01'))
            except ValueError:
                self.logger.error('calendar units are wrong: %s' %
                                  subMetadata['time#calendar'])
            else:
                time_coverage_start = (datetime.datetime(1950,1,1) +
                                   datetime.timedelta(float(timeDeltaStart)))

                if timeLength > 1:
                    timeDeltaEnd = Units.conform(timeValueStart,
                                          Units(subMetadata['time#units'],
                                                calendar=subMetadata['time#calendar']),
                                          Units('days since 1950-01-01'))
                else:
                    timeDeltaEnd = timeDeltaStart + 1
                time_coverage_end = (datetime.datetime(1950,1,1) +
                                     datetime.timedelta(float(timeDeltaEnd)))

        # set time_coverage_start if available
        if time_coverage_start is not None:
            self.dataset.SetMetadataItem('time_coverage_start',
                                    time_coverage_start.isoformat())
        # set time_coverage_end if available
        if time_coverage_end is not None:
            self.dataset.SetMetadataItem('time_coverage_end',
                                    time_coverage_end.isoformat())

