''' Nansat NetCDF-CF mapper

    Check CF-compliance of your files here:
    http://cfconventions.org/compliance-checker.html
'''
import warnings, os, datetime
import numpy as np
from dateutil.parser import parse

import gdal

from netCDF4 import Dataset

from nansat.vrt import VRT, GeolocationArray
from nansat.tools import WrongMapperError

class Mapper(VRT):

    def __init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs):

        if not filename.endswith('nc'):
            raise WrongMapperError
        
        gdal_metadata = self._remove_strings_in_metadata_keys(gdal_metadata)

        # Check conventions metadata
        if not gdal_metadata.has_key('Conventions') or not 'CF' in gdal_metadata['Conventions']:
            raise WrongMapperError

        # Create empty VRT dataset with geo-reference
        self._create_empty(gdal_dataset, gdal_metadata)
    
        # Add bands with metadata and corresponding values to the empty VRT
        self._create_bands(self._band_list(filename, gdal_metadata, *args, **kwargs))
        #xsize, ysize = self.ds_size(sub0)

        # Add GCMD/DIF compatible metadata in inheriting mappers

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
        tt = np.array([np.datetime64(t0 + datetime.timedelta(seconds=tn)) for
            tn in times])

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

    def _clean_band_metadata(self, band):
        band_metadata = band.GetMetadata_Dict()
        if 'PixelFunctionType' in band_metadata:
            band_metadata.pop('PixelFunctionType')
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
                bandName = bandMetadata.get('dods_variable', '')

            ## remove digits added by gdal when exporting to netcdf
            #if (len(bandName) > 0 and gdal_metadata.has_key('origin')
            #        and gdal_metadata['origin'].lower()=='nansat'):
            #    if bandName[-1:].isdigit():
            #        bandName = bandName[:-1]
            #    if bandName[-1:].isdigit():
            #        bandName = bandName[:-1]

        # if still no bandname, create one
        if len(bandName) == 0:
            bandName = 'band_%03d' % i

        dst['name'] = bandName

        return {'src': src, 'dst': dst}

    def _create_empty(self, gdal_dataset, gdal_metadata):
        subfiles = self.sub_filenames(gdal_dataset)
        if len(subfiles) == 0:
            raise WrongMapperError
        sub0 = gdal.Open(subfiles[0])

        super(Mapper, self).__init__(gdalDataset=sub0, srcMetadata=gdal_metadata)

        if not (self.dataset.GetGeoTransform() or self.geolocationArray.xVRT):
            # Set geolocation array from bands
            self.add_geolocationArray(
                    GeolocationArray(
                        [s for s in subfiles if 'longitude' in s or
                            'GEOLOCATION_X_DATASET' in s][0],
                        [s for s in subfiles if 'latitude' in s or
                            'GEOLOCATION_Y_DATASET' in s][0]))

        if not self.get_projection():
            # Get band projections
            projections = [gdal.Open(sub).GetProjection() for sub in subfiles if
                    gdal.Open(sub).GetProjection()]
            if not projections:
                raise WrongMapperError

            # Check that projection is the same for all bands
            assert all(proj==projections[0] for proj in projections)
            # Set projection
            self.dataset.SetProjection(projections[0])

    def _remove_strings_in_metadata_keys(self, gdal_metadata):
        if not gdal_metadata:
            raise WrongMapperError

        for key in gdal_metadata.keys():
            newkey = key.replace('NC_GLOBAL#', '')
            gdal_metadata[newkey] = gdal_metadata.pop(key)
        # Don't do this yet...
        #nans = 'NANSAT_'
        #for key in gdal_metadata.keys():
        #    if nans in key:
        #        gdal_metadata[key.replace(nans, '')] = gdal_metadata.pop(newkey)
        #        #gdal_metadata['origin'] = 'NANSAT'

        return gdal_metadata

    def sub_filenames(self, gdal_dataset):
        # Get filenames of subdatasets
        sub_datasets = gdal_dataset.GetSubDatasets()
        # add subdatasets with Projection only
        sub_fnames = []
        for f in sub_datasets:
            if gdal.Open(f[0]).GetProjection():
                sub_fnames.append(f[0])
        return sub_fnames

    #def ds_size(self, ds):
    #    return ds.RasterXSize, ds.RasterYSize
