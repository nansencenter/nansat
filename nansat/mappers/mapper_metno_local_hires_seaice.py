# Name:         mapper_metno_hires_seaice.py
# Purpose:      Nansat mapping for high resolution sea ice
#               from met.no Thredds server
# Authors:      Knut-Frode Dagestad
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

# High resolution (1 km) manual ice concentration, based on SAR imagery
#
# Mapper is called with keyword (fake filename):
#   'metno_hires_seaice_20140109'
# and [optional] keyword: iceFolder = <local folder
# where netCDF files with ice are stored>
#
# The closest available data within +/- 3 days is returned
from __future__ import print_function, absolute_import, unicode_literals
import sys
import os
from datetime import datetime, timedelta

from nansat.tools import gdal, ogr
from nansat.exceptions import WrongMapperError
from nansat.vrt import VRT
import mapper_generic as mg


class Mapper(mg.Mapper):
    """Create VRT with mapping of WKV for Met.no seaice"""

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):
        """Create VRT"""

        try:
            ice_folder_name = kwargs['iceFolder']
        except:
            #iceFolderName = '/vol/istjenesten/data/metnoCharts/'
            ice_folder_name = '/vol/data/metnoCharts/'

        keyword_base = 'metno_local_hires_seaice'

        if filename[0:len(keyword_base)] != keyword_base:
            raise WrongMapperError

        keyword_time = filename[len(keyword_base)+1:]
        requested_time = datetime.strptime(keyword_time, '%Y%m%d')
        # Search for nearest available file, within the closest 3 days
        found_dataset = False
        for delta_day in [0, -1, 1, -2, 2, -3, 3]:
            valid_time = (requested_time + timedelta(days=delta_day) +
                          timedelta(hours=15))
            filename = (ice_folder_name + 'ice_conc_svalbard_' +
                        valid_time.strftime('%Y%m%d1500.nc'))
            if os.path.exists(filename):
                print('Found file:')
                print(filename)
                gdal_dataset = gdal.Open(filename)
                gdal_metadata = gdalDataset.GetMetadata()
                mg.Mapper.__init__(self, filename, gdal_dataset, gdal_metadata)
                found_dataset = True
                # Modify GeoTransform from netCDF file
                # - otherwise a shift is seen!
                self.dataset.SetGeoTransform(
                    (-1243508 - 1000, 1000, 0, -210526 - 7000, 0, -1000))
                break  # Data is found for this day

        if found_dataset is False:
            AttributeError("No local Svalbard-ice files available")
            sys.exit()
