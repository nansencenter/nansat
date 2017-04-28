import os, glob
import numpy as np
from dateutil.parser import parse

import gdal

from nansat.vrt import VRT, GeolocationArray
from nansat.tools import WrongMapperError
from nansat.nsr import NSR
from nansat.mappers.mapper_netcdf_cf import Mapper as NCCF_mapper


class Mapper(VRT):
    '''
        Create VRT with mapping of Sentinel-1 OCN product
    '''

    def __init__(self, filename, gdal_dataset, gdal_metadata, product='RVL',
            subswath=0, *args, **kwargs):
        '''
        Parameters
        ----------
        product: string
            Sentinel-1 level-2 ocean product type/component, i.e. ocean swell
            spectra (OSW), ocean wind field (OWI), or radial surface velocity
            (RVL) (RVL is the default)
        #GCP_COUNT : int
        #    number of GCPs along each dimention
        '''
        if not product.lower() in ['rvl', 'owi', 'osw']:
            raise WrongMapperError

        mdir = os.path.join(filename, 'measurement')
        if not os.path.exists(mdir):
            raise WrongMapperError

        ncFile = glob.glob(mdir+'/*.nc')
        if not ncFile:
            raise WrongMapperError
        ncFile = ncFile[0] # TODO: refine... could perhaps be more than one(?)

        gdal_dataset = gdal.Open(ncFile)
        gdal_metadata = gdal_dataset.GetMetadata()
        # Check if it is Sentinel-1 (or ASAR) level-2 (in S1 data format)
        if not gdal_metadata:
            raise WrongMapperError

        gdal_metadata = self._remove_strings_in_metadata_keys(gdal_metadata)
        title = gdal_metadata['title']

        # Raise error if it is not Sentinel-1
        if not 'Sentinel-1' in title:
            raise WrongMapperError

        subfiles = self.sub_filenames(gdal_dataset)
        sub0 = gdal.Open(subfiles[0])

        lonfn = [s for s in subfiles if product.lower()+'Lon' in s][0]
        lon = gdal.Open(lonfn).ReadAsArray()[:,:,subswath]
        lon[lon==-999.]=np.nan
        latfn = [s for s in subfiles if product.lower()+'Lat' in s][0]
        lat = gdal.Open(latfn).ReadAsArray()[:,:,subswath]
        lat[lat==-999.]=np.nan
        super(Mapper, self).__init__(srcMetadata=gdal_metadata, lat=lat,
                lon=lon)

        bands = [b.split(':')[-1] for b in subfiles if 'rvl' in b]

        metaDict = []
        self.bandVRTs = {}
        for fn in subfiles:
            if ('Lon' in fn or 'Lat' in fn or not 'rvl' in fn or
                    'ZeroDopplerTime' in fn):
                continue
            band_name = fn.split(':')[-1]
            band_array = gdal.Open(fn).ReadAsArray()[:,:,subswath]
            band_array[band_array==-999.] = np.nan
            self.bandVRTs[band_name] = VRT(array=band_array, lat=lat, lon=lon)
            if band_name=='rvlIncidenceAngle':
                dst = {'wkv': 'angle_of_incidence', 'name': 'incidence_angle'}
            elif band_name=='rvlDcObs':
                dst = {'wkv':
                    'surface_backwards_doppler_centroid_frequency_shift_of_radar_wave'}
            elif band_name=='rvlNrcs':
                dst = {'wkv':
                        'surface_backwards_scattering_coefficient_of_radar_wave'}
            elif band_name=='rvlDcGeo':
                dst = {'wkv':
                    'predicted_surface_backwards_doppler_frequency_shift_of_radar_wave'}
            else:
                dst = {'name': band_name}
            metaDict.append({
                'src': {
                    'SourceFilename': self.bandVRTs[band_name].fileName,
                    'SourceBand': 1},
                'dst': dst
            })

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)
