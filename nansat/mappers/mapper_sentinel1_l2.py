# Name:     mapper_s1a_l2.py
# Purpose:      Mapping for Sentinel-1 level-2 data
#
# Author:       Morten Wergeland Hansen
# Modified: Morten Wergeland Hansen
#
# Created:  13.03.2014
# Last modified:17.07.2014 11:57
# Copyright:    (c) NERSC
# License:      This file is part of NANSAT. NANSAT is free software: you can
#               redistribute it and/or modify it under the terms of the GNU
#               General Public License as published by the Free Software
#               Foundation, version 3 of the License.
#               http://www.gnu.org/licenses/gpl-3.0.html This program is
#               distributed in the hope that it will be useful, but WITHOUT ANY
#               WARRANTY without even the implied warranty of MERCHANTABILITY
#               or FITNESS FOR A PARTICULAR PURPOSE.
import os
from dateutil.parser import parse

import gdal

from nansat.vrt import VRT
from nansat.exceptions import WrongMapperError
from nansat.nsr import NSR


class Mapper(VRT):
    '''
        Create VRT with mapping of Sentinel-1A stripmap mode (S1A_SM)
    '''

    def __init__(self, filename, gdalDataset, gdalMetadata,
                 product_type='RVL', GCP_COUNT=10, **kwargs):
        '''
        Parameters
        ----------
        product_type: string
            Sentinel-1 level-2 ocean product type/component, i.e. ocean swell
            spectra (OSW), ocean wind field (OWI), or radial surface velocity
            (RVL) (RVL is the default)
        GCP_COUNT : int
            number of GCPs along each dimention
        '''
        fPathName, fExt = os.path.splitext(filename)

        # List of Sentinel-1 level-2 components
        unwanted_product_components = ['osw', 'owi', 'rvl']
        # Remove requested 'product_type' from list of unwanted
        unwanted_product_components.pop(unwanted_product_components.index(
                                        product_type.lower()))

        # Check if it is Sentinel-1 (or ASAR) level-2 (in S1 data format)
        if not gdalMetadata or not 'NC_GLOBAL' in gdalMetadata.keys():
            raise WrongMapperError(filename)
        else:
            title = gdalMetadata['NC_GLOBAL#TITLE']

        # Raise error if it is not Sentinel-1 format
        if not 'Sentinel-1' or 'ASA' in title:
            raise WrongMapperError(filename)

        metadata = {}
        for key, val in gdalMetadata.iteritems():
            new_key = key.split('#')[-1]
            metadata[new_key] = val

        subDatasets = gdalDataset.GetSubDatasets()
        filenames = [f[0] for f in subDatasets]

        rm_bands = []
        # Find all data that is not relevant for the selected product type
        # and get bands of longitude, latitude and zero doppler time
        for i, f in enumerate(filenames):
            if f.split(':')[-1][:3] in unwanted_product_components:
                rm_bands.append(i)
            if 'Lon' in f.split(':')[-1]:
                lon_ds = gdal.Open(f)
                rm_bands.append(i)
            if 'Lat' in f.split(':')[-1]:
                lat_ds = gdal.Open(f)
                rm_bands.append(i)
            if 'ZeroDopplerTime' in f.split(':')[-1]:
                zdt_ds = gdal.Open(f)
                rm_bands.append(i)
        # Remove bands in rm_bands from the list of bands to add to the Nansat
        # object
        filenames = [f for i, f in enumerate(filenames) if not i in rm_bands]
           #     (
           # 'Lon' in f.split(':')[-1] or
           # 'Lat' in f.split(':')[-1] or
           # 'ZeroDopplerTime' in f.split(':')[-1] )]

        # create empty VRT dataset
        self._init_from_gdal_dataset(gdal.Open(subDatasets[0][0]), metadata=metadata)

        # The zero Doppler time grid is 3-dimensional - the last dimension is a
        # char array with the time as year, month, day, etc.
        # Will not bother with it yet...
        #for iBand in range(zdt_ds.RasterCount):
        #    subBand = zdt_ds.GetRasterBand(iBand+1)

        XSize = lon_ds.RasterXSize
        YSize = lon_ds.RasterYSize

        # get projection from the lon and lat datasets
        longitude = lon_ds.ReadAsArray()
        latitude = lat_ds.ReadAsArray()
        # estimate step of GCPs
        step0 = max(1, int(float(latitude.shape[0]) / GCP_COUNT))
        step1 = max(1, int(float(latitude.shape[1]) / GCP_COUNT))
        self.logger.debug('gcpCount: >%s<, %d %d %f %d %d',
                          title,
                          latitude.shape[0], latitude.shape[1],
                          GCP_COUNT, step0, step1)

        # estimate pixel/line step of the geolocation arrays
        pixelStep = 1
        lineStep = 1
        self.logger.debug('pixel/lineStep %f %f' % (pixelStep, lineStep))
        # generate list of GCPs
        dx = .5
        dy = .5
        gcps = []
        k = 0
        for i0 in range(0, latitude.shape[0], step0):
            for i1 in range(0, latitude.shape[1], step1):
                # create GCP with X,Y,pixel,line from lat/lon matrices
                lon = float(longitude[i0, i1])
                lat = float(latitude[i0, i1])
                if (lon >= -180 and lon <= 180 and lat >= -90 and lat <= 90):
                    gcp = gdal.GCP(lon, lat, 0, i1 * pixelStep + dx,
                                   i0 * lineStep + dy)
                    self.logger.debug('%d %d %d %f %f',
                                      k, gcp.GCPPixel, gcp.GCPLine,
                                      gcp.GCPX, gcp.GCPY)
                    gcps.append(gcp)
                    k += 1

        # append GCPs and lat/lon projection to the vsiDataset
        self.dataset.SetGCPs(gcps, NSR().wkt)

        # define band specific parameters
        metaDict = []
        geoFileDict = {}
        xDatasetSource = ''
        yDatasetSource = ''
        for i, filename in enumerate(filenames):
            band = gdal.Open(filename)

            # check that the band size is the same size as the latitude and
            # longitude grids
            if (band.RasterXSize != XSize or
                    band.RasterYSize != YSize):
                raise IndexError(('Size of sub-dataset is different from size '
                                  'of longitude and latitude grids'))

            bandMetadata = band.GetMetadata()
            # generate src metadata
            src = {'SourceFilename': filename,
                   'SourceBand': 1
                   }

            # Generate dst metadata
            short_name = filename.split(':')[-1]
            dst = {'name': short_name,
                   'short_name': short_name,
                   'long_name': bandMetadata[short_name+'#long_name'],
                   'units': bandMetadata[short_name+'#units'],
                   #'wkv': ,
                   }

            # append band with src and dst dictionaries
            metaDict.append({'src': src, 'dst': dst})

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        metaDict = []
        for i in range(self.dataset.RasterCount):
            if 'Nrcs' in self.dataset.GetRasterBand(i+1).GetMetadata()['name']:
                metaDict.append({
                    'src': {'SourceFilename': (self.dataset.GetRasterBand(i+1).
                                GetMetadata()['SourceFilename']),
                            'SourceBand': 1
                            },
                    'dst': {
                        'short_name': 'sigma0',
                        'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                        'PixelFunctionType':  'dB2pow',
                        'polarization': (self.dataset.
                                         GetMetadata()['POLARISATION']),
                        'suffix': self.dataset.GetMetadata()['POLARISATION'],
                        'dataType': 6,
                    }
                })

        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        # set time
        self.dataset.SetMetadataItem('time_coverage_start',
            parse(self.dataset.GetMetadata()['SOURCE_ACQUISITION_UTC_TIME']).isoformat())
