# Name:         mapper_opendap_ostia.py
# Purpose:      Nansat mapping for  GHRSST Level 4 OSTIA Global Foundation Sea Surface
#               Temperature Analysis
# Author:       Artem Moiseev
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from nansat.mappers.opendap import Opendap
from nansat.nsr import NSR
import pythesint as pti
import os
from datetime import datetime
import numpy as np
import json
from netCDF4 import Dataset


class Mapper(Opendap):

    baseURLs = [
        'https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/ghrsst/data/L4/GLOB/UKMO/OSTIA',
        'https://opendap.jpl.nasa.gov:443/opendap/OceanTemperature/ghrsst/data/L4/GLOB/UKMO/OSTIA'
    ]

    timeVarName = 'time'
    xName = 'lon'
    yName = 'lat'
    timeCalendarStart = '1981-01-01'
    srcDSProjection = NSR().wkt

    def __init__(self, filename, gdal_dataset, gdal_metadata, date=None,
                 ds=None, bands=None, cachedir=None, *args, **kwargs):

        self.test_mapper(filename)
        timestamp = date if date else self.get_date(filename)
        ds = Dataset(filename)
        self.create_vrt(filename, gdal_dataset, gdal_metadata, timestamp, ds, bands, cachedir)
        self.dataset.SetMetadataItem('entry_title', str(ds.getncattr('title')))
        self.dataset.SetMetadataItem('data_center', json.dumps(pti.get_gcmd_provider('UK/MOD/MET')))
        self.dataset.SetMetadataItem('ISO_topic_category',
                pti.get_iso19115_topic_category('oceans')['iso_topic_category'])
        self.dataset.SetMetadataItem('gcmd_location', json.dumps(pti.get_gcmd_location('sea surface')))

        #mm = pti.get_gcmd_instrument('amsr-e')
        #ee = pti.get_gcmd_platform('aqua')
        #self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        #self.dataset.SetMetadataItem('platform', json.dumps(ee))
        self.dataset.SetMetadataItem('platform/instrument',
                json.dumps(self.get_platform_and_instrument_list(ds)))

    @staticmethod
    def get_date(filename):
        """Extract date and time parameters from filename and return
        it as a formatted (isoformat) string

        Parameters
        ----------

        filename: str
            nn

        Returns
        -------
            str, YYYY-mm-ddThh:MMZ

        """
        _, filename = os.path.split(filename)
        t = datetime.strptime(filename.split('-')[0], '%Y%m%d')
        return datetime.strftime(t, '%Y-%m-%dT%H:%M:00Z')

    def convert_dstime_datetimes(self, ds_time):
        """Convert time variable to np.datetime64"""
        ds_datetimes = np.array(
            [(np.datetime64(self.timeCalendarStart).astype('M8[s]')
              + np.timedelta64(int(sec), 's').astype('m8[s]')) for sec in ds_time]).astype('M8[s]')
        return ds_datetimes

    @staticmethod
    def get_platform_and_instrument_list(ds):
        """ This method uses the source_data in the OPeNDAP dataset to select the platforms
        and instruments. It checks the hardcoded dictionary, pi, to find the platform and instrument
        for products given in the source_data. The reason for that is that the items in source_data
        are not in the GCMD keywords but rather refer to products of certain instruments. If you
        search the internet for the dictionary keys, you'll find a dataset description which
        includes the listed platforms and instruments.
        """
        pi = {
                'AVHRR18_G-NAVO-L2P-V1.0': [pti.get_gcmd_platform('noaa-18'),
                    pti.get_gcmd_instrument('avhrr-3')],
                'AVHRR19_G-NAVO-L2P-V1.0': [pti.get_gcmd_platform('noaa-19'),
                    pti.get_gcmd_instrument('avhrr')],
                'AVHRR_SST_METOP_B-OSISAF-L2P-V1.0': [pti.get_gcmd_platform('metop-b'),
                    pti.get_gcmd_instrument('avhrr')],
                'VIIRS_NPP-OSPO-L2P-V2.3': [pti.get_gcmd_platform('suomi-npp'),
                    pti.get_gcmd_instrument('viirs')],
                'AMSR2-REMSS-L2P-V07.2': [pti.get_gcmd_platform('gcom-w1'),
                    pti.get_gcmd_instrument('amsr2')],
                'GOES13-OSISAF-L3C-V1.0': [pti.get_gcmd_platform('goes-16'),
                    pti.get_gcmd_instrument('abi')],
                'SEVIRI_SST-OSISAF-L3C-V1.0': [pti.get_gcmd_platform('msg'),
                    pti.get_gcmd_instrument('seviri')],
                'OSISAF_ICE': [pti.get_gcmd_platform('earth observation satellites'),
                    pti.get_gcmd_instrument('Imaging Spectrometers/Radiometers')],
                'NCEP_ICE': [pti.get_gcmd_platform('ncep-gfs'),
                    pti.get_gcmd_instrument('computer')],
                'AMSRE': [pti.get_gcmd_platform('aqua'), pti.get_gcmd_instrument('amsr-e')],
                'ATS_NR_2P': [pti.get_gcmd_platform('envisat'), pti.get_gcmd_instrument('aatsr')],
                'AVHRR18_G': [pti.get_gcmd_platform('noaa-18'), pti.get_gcmd_instrument('avhrr-3')],
                'AVHRR17_NAR': [pti.get_gcmd_platform('noaa-17'), pti.get_gcmd_instrument('avhrr')],
                'AVHRR18_NAR': [pti.get_gcmd_platform('noaa-18'), pti.get_gcmd_instrument('avhrr-3')],
                'SEVIRI': [pti.get_gcmd_platform('msg'), pti.get_gcmd_instrument('seviri')],
                'TMI': [pti.get_gcmd_platform('trmm'), pti.get_gcmd_instrument('tmi')],
        }
        pi_list = []
        # Here, we may have a dilemma: for example, the dataset at
        # https://opendap.jpl.nasa.gov:443/opendap/OceanTemperature/ghrsst/data/L4/GLOB/UKMO/OSTIA/2008/002/20080102-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc.bz2
        # has source data AVHRR18_G and AVHRR18_NAR. I don't know the difference between them but
        # presume both are noaa18/AVHRR-3. It means duplication. Might not be a problem in Nansat,
        # though, and it could be easy to solve in django-geo-spaas...
        for source_data in ds.source_data.split(','):
            pi_list.append(pi[source_data.strip()])
        return pi_list

