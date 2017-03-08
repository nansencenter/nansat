''' Nansat NetCDF-CF mapper

    Check CF-compliance of your files here:
    http://cfconventions.org/compliance-checker.html
'''
import gdal

from nansat.vrt import VRT, GeolocationArray
from nansat.tools import WrongMapperError

class Mapper(VRT):

    def __init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs):

        gdal_metadata = self.replace_nc_global(gdal_metadata)

        # Create empty VRT dataset with geo-reference
        self._create_empty(gdal_dataset, gdal_metadata)
    
        import ipdb
        ipdb.set_trace()

        # Add bands with metadata and corresponding values to the empty VRT
        self._create_bands(self.band_list(*args, **kwargs))
        #xsize, ysize = self.ds_size(sub0)


        print 'hei'

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

    def band_list(self):
        return 0

    def replace_nc_global(self, gdalMetadata):
        if not gdalMetadata:
            raise WrongMapperError

        for key in gdalMetadata.keys():
             gdalMetadata[key.replace('NC_GLOBAL#', '')] = gdalMetadata.pop(
                                                                key)
        return gdalMetadata

    def sub_filenames(self, gdal_dataset):
        # Get filenames of subdatasets
        sub_datasets = gdal_dataset.GetSubDatasets()
        return [f[0] for f in sub_datasets]

    #def ds_size(self, ds):
    #    return ds.RasterXSize, ds.RasterYSize
