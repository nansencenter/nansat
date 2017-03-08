''' Nansat NetCDF-CF mapper

    Check CF-compliance of your files here:
    http://cfconventions.org/compliance-checker.html
'''
import gdal

from nansat.vrt import VRT, GeolocationArray
from nansat.tools import WrongMapperError

class Mapper(VRT):

    def __init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs):

        gdal_metadata = self._remove_strings_in_metadata_keys(gdal_metadata)

        # Create empty VRT dataset with geo-reference
        self._create_empty(gdal_dataset, gdal_metadata)
    
        # Add bands with metadata and corresponding values to the empty VRT
        self._create_bands(self._band_list(filename, gdal_metadata, *args, **kwargs))
        #xsize, ysize = self.ds_size(sub0)


    def _band_list(self, filename, gdal_metadata, time=None, bands=[]):
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
                band_metadata = band.GetMetadata_Dict()
                if 'PixelFunctionType' in band_metadata:
                    band_metadata.pop('PixelFunctionType')

                # Generate source metadata
                src = {'SourceFilename': filename,
                       'SourceBand': band_num}

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

                    # remove digits added by gdal when exporting to netcdf
                    if (len(bandName) > 0 and gdal_metadata.has_key('origin')
                            and gdal_metadata['origin'].lower()=='nansat'):
                        if bandName[-1:].isdigit():
                            bandName = bandName[:-1]
                        if bandName[-1:].isdigit():
                            bandName = bandName[:-1]

                # if still no bandname, create one
                if len(bandName) == 0:
                    bandName = 'band_%03d' % i

                dst['name'] = bandName

                # append band with src and dst dictionaries
                metadictlist.append({'src': src, 'dst': dst})

        return metadictlist

    def _create_empty(self, gdal_dataset, metadata):
        subfiles = self.sub_filenames(gdal_dataset)
        sub0 = gdal.Open(subfiles[0])

        super(Mapper, self).__init__(self, sub0, srcMetadata=metadata)

        if not self.get_projection():
            # Set geolocation array from bands
            self.add_geolocationArray(
                    GeolocationArray(
                        [s for s in subfiles if 'longitude' in s or
                            'GEOLOCATION_X_DATASET' in s][0],
                        [s for s in subfiles if 'latitude' in s or
                            'GEOLOCATION_Y_DATASET' in s][0]))
            # Get band projections
            projections = [gdal.Open(sub).GetProjection() for sub in subfiles if
                    gdal.Open(sub).GetProjection()]
            # Check that projection is the same for all bands
            assert all(proj==projections[0] for proj in projections)
            # Set projection
            self.dataset.SetProjection(projections[0])

    def _remove_strings_in_metadata_keys(self, gdalMetadata):
        if not gdalMetadata:
            raise WrongMapperError

        for key in gdalMetadata.keys():
            newkey = key.replace('NC_GLOBAL#', '')
            gdalMetadata[newkey] = gdalMetadata.pop(key)
        nans = 'NANSAT_'
        for key in gdalMetadata.keys():
            if nans in key:
                gdalMetadata[key.replace(nans, '')] = gdalMetadata.pop(newkey)
                gdalMetadata['origin'] = 'NANSAT'

        return gdalMetadata

    def sub_filenames(self, gdal_dataset):
        # Get filenames of subdatasets
        sub_datasets = gdal_dataset.GetSubDatasets()
        return [f[0] for f in sub_datasets]

    #def ds_size(self, ds):
    #    return ds.RasterXSize, ds.RasterYSize
