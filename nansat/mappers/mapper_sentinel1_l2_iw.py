# ------------------------------------------------------------------------------
# Name:     mapper_sentinel1_l2_iw.py
# Purpose:  Mapper for Sentinel-1 IW L2 OCN products with support of the RLV and OWI
#
# Author:       Artem Moiseev
#
# Created:  28.05.2021
# Copyright:    (c) NERSC
# License: GPL V3
# ------------------------------------------------------------------------------


import os
import re   
from datetime import datetime

from nansat.vrt import VRT
from nansat import Domain, NSR
from nansat.exceptions import WrongMapperError

import json
import pythesint as pti
from netCDF4 import Dataset, date2num
import numpy as np


class Mapper(VRT):
    """"""

    def __init__(self, file_src, gdal_dataset, gdal_metadata, *args, **kwargs):
        # Check if an input file can be read with this mapper, else raise error
        # Check the filename
        Mapper.check_input(file_src)

        # Check if the subswath number [1, 2, 3] is provided
        # Subswath number is required for mapping the RVL products
        if 'subswath_n' in kwargs.keys():
            subswath_id = kwargs['subswath_n'] - 1
        else:
            subswath_id = 0
        # Check if the product type is provided: rvl 
        if 'product_type' in kwargs.keys():
            product = kwargs['product_type']
        else:
            product = 'rvl'

        if 'dop_time' in kwargs.keys():
            dop_time = kwargs['dop_time']
        else:
            dop_time = False
        # Read the NetCDF dataset
        nc_dataset = Dataset(file_src)
        # Generate a domain for the required subswath
        # (1) Import lon lat arrays for the required subswath from the file
        if product == 'rvl':
            lon = nc_dataset.variables[product + 'Lon'][:, :, subswath_id]
            lat = nc_dataset.variables[product + 'Lat'][:, :, subswath_id]
        elif product == 'owi':
            lon = nc_dataset.variables[product + 'Lon'][:, :]
            lat = nc_dataset.variables[product + 'Lat'][:, :]
        else:
            raise ValueError

        # (2) Filter part not covered with data following the data mask from the product
        lon = self.filter_data(lon)
        lat = self.filter_data(lat)
        # (3) Generale the subswath domain from the lon lat arrays
        self._init_from_lonlat(lon=lon, lat=lat)
        # Get GCPs from the created domain and set them for the input file
        dom = Domain.from_lonlat(lon, lat)
        self.dataset.SetGCPs(dom.vrt.dataset.GetGCPs(), NSR().wkt)

        # Add bands
        for subdataset_name in nc_dataset.variables.keys():
            # The ODL Sentinel-1 data contain several products. eact var name 
            # starts with the name of the generic product e.g.: rvlDoppler or oswAngularBinSize
            # Use only required type of products e.g., rvl
            if not str.startswith(subdataset_name, product):
                continue
            # Import band data
            band = nc_dataset.variables[subdataset_name]

            if product == 'rvl':
                if band.dimensions[:2] != ('rvlAzSize', 'rvlRaSize'):
                    continue
                elif subdataset_name == 'rvlZeroDopplerTime':
                    if dop_time:
                        band_data, band_metadata = self.process_zero_doppler_time_band(
                            nc_dataset, band, subswath_id)
                else:
                    band_data = self.filter_data(nc_dataset.variables[subdataset_name][:, :, subswath_id])
                    band_metadata = {}
            
            elif product == 'owi':
                # Process only part of the OWI product which prowided on the OWI grid 
                # Some part of the OWI product is provided on the OSW grid and will be ignored
                if len(band.dimensions) < 2:
                    continue 
                elif band.dimensions[:2] != ('owiAzSize', 'owiRaSize'):
                    continue
                elif len(band.dimensions) == 2:
                    band_data = self.filter_data(nc_dataset.variables[subdataset_name][:, :])
                    band_metadata = {}
                else:
                    continue
            else:
                raise ValueError

            # Create band VRT
            # NOTE: The try/except routine was added to fix the issue with the mask in the UssX/UssY bands
            # TODO: Fix the issue. Use a SAR mask while filtering the data instead of using the mask from each NetCDF band
            try:
                band_vrt = VRT.from_array(band_data)
            except RuntimeError:
                continue
            # Add VRT to the list of all dataset vrts
            self.band_vrts[subdataset_name + 'VRT'] = band_vrt
            src = {'SourceFilename': self.band_vrts[subdataset_name + 'VRT'].filename, 
                   'SourceBand': 1}
            # Add band specific metadata
            dst = {'name': subdataset_name}
            if len(band_metadata.keys()) != 0:
                for key, value in band_metadata.items():
                    dst[key] = value

            for subdataset_att in nc_dataset.variables[subdataset_name].ncattrs():
                dst[subdataset_att] = nc_dataset.variables[subdataset_name].getncattr(subdataset_att)
            # Create band
            self.create_band(src, dst)
            self.dataset.FlushCache()

        # Set generic metadata
        try:
            platform = nc_dataset.getncattr('missionName').replace('S', 'Sentinel-')
            self.dataset.SetMetadataItem('platform', json.dumps(pti.get_gcmd_platform(platform)))
            self.dataset.SetMetadataItem('instrument', json.dumps(pti.get_gcmd_instrument('SAR')))
            self.dataset.SetMetadataItem('data_center', json.dumps(
                pti.get_gcmd_provider('FR/IFREMER/CERSAT')))
            self.dataset.SetMetadataItem('ISO_topic_category',
                                         json.dumps(pti.get_iso19115_topic_category('Oceans')))
            self.dataset.SetMetadataItem('summary', '')
            self.dataset.SetMetadataItem('entry_title', '')
        except json.JSONDecodeError:
            pass
        # Add metadata from the product
        for key, value in nc_dataset.__dict__.items():
            if isinstance(value, str):
                self.dataset.SetMetadataItem(key, value)
            elif isinstance(value, np.ndarray):
                self.dataset.SetMetadataItem(key, ','.join([str(a) for a in value]))
            else:
                continue

        # Set time coverage metadata
        time_coverage_start = datetime.strptime(nc_dataset.getncattr('firstMeasurementTime')[:-1], '%Y-%m-%dT%H:%M:%S')
        time_coverage_end = datetime.strptime(nc_dataset.getncattr('lastMeasurementTime')[:-1], '%Y-%m-%dT%H:%M:%S')
        self.dataset.SetMetadataItem('time_coverage_start', time_coverage_start.isoformat())
        self.dataset.SetMetadataItem('time_coverage_end', time_coverage_end.isoformat())

    def process_zero_doppler_time_band(self, nc_ds, band, subswath_id, calendar_start='seconds since 2010-01-01'):
        zero_dop_time_raw = np.char.decode(self.filter_data(band[:, :, :, subswath_id]))
        timeseries = zero_dop_time_raw.reshape(
            (zero_dop_time_raw.shape[0] * zero_dop_time_raw.shape[1], -1))
        # Assemble raw timestamp strings and convert it to datetime objects
        zero_dop_time_datetime = np.array(list(map(self._join_test, timeseries)))

        zero_dop_time_datetime_2d = zero_dop_time_datetime.reshape(zero_dop_time_raw.shape[:2])

        band_data = date2num(zero_dop_time_datetime_2d, calendar_start, calendar='standard')
        band_metadata = {'units': calendar_start,
                         'long_name': band.long_name}
        return band_data, band_metadata

    def _convert_time(self, time_str):
        try:
            time_stamp = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S.%f')
        except:
            time_stamp = datetime(1999, 1, 1, 0, 0, 0)
        return time_stamp

    def _join_test(self, arr):
        """Assemble the time string"""
        return self._convert_time(''.join(arr)[:-3])

    def filter_data(self, arr):
        # Remove data covered by mask from the array
        if isinstance(arr.mask, np.ndarray):
            unmasked_data = arr.data[:, arr.mask[0] == False]
        else:
            unmasked_data = arr
        return unmasked_data

    @staticmethod
    def check_input(file_src):
        """Check that input file (according to file_src) can be read by the mapper
        Parameters
        ----------
        file_src : str
            path of input Sentinel-1 IW file
        """
        # Extract filename and file format for the src
        filename_base, file_format = Mapper.parse_filename(file_src.strip())

        if file_format != '.nc':
            raise WrongMapperError('Only netCDF format is supported')
        elif not re.match(r's1.-iw-ocn-vv*', filename_base):
            raise WrongMapperError()

    @staticmethod
    def parse_filename(file_src):
        """Parse input file src to extract file name and format
        Parameters
        ----------
        file_src : str
            path of input Sentinel-1 IW file
        Returns
        -------
        filename_base: str
            filename without format (i.e. without '.nc')
        file_format: str
            format of the file (i.e '.nc')
        """
        _, filename = os.path.split(file_src)
        filename_base, file_format = os.path.splitext(filename)
        return filename_base, file_format
