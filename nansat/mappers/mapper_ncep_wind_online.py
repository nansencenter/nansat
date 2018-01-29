# Name:         mapper_ncep_wind_online.py
# Purpose:      Nansat mapping for NCEP GFS model data, stored online
# Author:       Knut-Frode Dagestad, Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#
# Mapper searches two online archives for NCEP GFS grib files
# covering the requested time, and downloads using curl, if found:
#    1. ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/ (~last month)
#    2. http://nomads.ncdc.noaa.gov/data/gfs4/ (back to June 2012, with holes)
#
# Usage:
#    w = Nansat('ncep_wind_online:YYYYMMDDHHMM')
# Example:
#    w = Nansat('ncep_wind_online:201405011000')

from __future__ import absolute_import, print_function, division, unicode_literals

import os
import sys
from datetime import datetime, timedelta

from nansat.vrt import VRT
from nansat.exceptions import WrongMapperError
from nansat.nansat import Nansat

# Place to store downloads - this can be changed via the "outFolder" argument
# to Mapper.__init__
downloads = os.path.join(os.path.expanduser('~'), 'ncep_gfs_downloads')


class Mapper(VRT, object):
    """VRT with mapping of WKV for NCEP GFS"""

    def __init__(self, filename, gdalDataset, gdalMetadata,
                 outFolder=downloads, **kwargs):
        """Create NCEP VRT"""

        if not os.path.exists(outFolder):
            os.mkdir(outFolder)

        ##############
        # Get time
        ##############
        keyword_base = 'ncep_wind_online'
        if filename[0:len(keyword_base)] != keyword_base:
            raise WrongMapperError

        time_str = filename[len(keyword_base)+1::]
        time = datetime.strptime(time_str, '%Y%m%d%H%M')
        print(time)

        ########################################
        # Find and download online grib file
        ########################################
        # Find closest 6 hourly modelrun and forecast hour
        model_run_hour = round((time.hour + time.minute/60.)/6)*6
        nearest_model_run = (datetime(time.year, time.month, time.day)
                             + timedelta(hours=model_run_hour))
        if sys.version_info < (2, 7):
            td = (time - nearest_model_run)
            forecast_hour = (td.microseconds +
                             (td.seconds + td.days * 24 * 3600)
                             * 10 ** 6) / 10 ** 6 / 3600.
        else:
            forecast_hour = (time - nearest_model_run).total_seconds() / 3600.
        if model_run_hour == 24:
            model_run_hour = 0
        if forecast_hour < 1.5:
            forecast_hour = 0
        else:
            forecast_hour = 3

        #########################################################
        # Try first to get NRT data from
        # ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/
        # - avaliable approximately the latest month
        #########################################################
        url = ('ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/' + 'gfs.'
               + nearest_model_run.strftime('%Y%m%d') + '%.2d' % model_run_hour + '/gfs.t' + '%.2d'
               % model_run_hour + 'z.master.grbf' +  '%.3d' % forecast_hour + '.10m.uv.grib2')
        out_filename = os.path.join(outFolder,
                                   ('ncep_gfs_' +
                                    nearest_model_run.strftime('%Y%m%d_%HH_') +
                                    '%.2d' % forecast_hour +
                                    '.10m.uv.grib2'))
        if os.path.exists(out_filename):
            print('NCEP wind is already downloaded: ' + out_filename)
        else:
            os.system('curl -so ' + out_filename + ' ' + url)
            if os.path.exists(out_filename):
                print('Downloaded ' + out_filename)
            else:
                print('NRT GRIB file not available: ' + url)
                #########################################################
                # If NRT file not available, search in long term archive
                #########################################################
                url = ('http://nomads.ncdc.noaa.gov/data/gfs4/' +
                       nearest_model_run.strftime('%Y%m/%Y%m%d/'))
                basename = ('gfs_4_' + nearest_model_run.strftime('%Y%m%d_') +
                            nearest_model_run.strftime('%H%M_') +
                            '%.3d' % forecast_hour)
                filename = basename + '.grb2'
                out_filename = os.path.join(outFolder, filename)
                print('Downloading ' + url + filename)

                # Download subset of grib file
                mapper_dir = os.path.dirname(os.path.abspath(__file__))
                get_inv = os.path.join(mapper_dir, 'get_inv.pl')
                if not os.path.isfile(get_inv):
                    raise IOError('%s: File not found' % get_inv)
                get_grib = os.path.join(mapper_dir, 'get_grib.pl')

                if not os.path.isfile(get_grib):
                    raise IOError('%s: File not found' % get_grib)

                if not os.path.isfile(out_filename):
                    command = (get_inv + ' ' + url + basename +
                               '.inv | egrep "(:UGRD:10 m |:VGRD:10 m )" | ' +
                               get_grib + ' ' + url + filename + ' ' + out_filename)
                    os.system(command)
                    if os.path.isfile(out_filename):
                        print('Downloaded ' + filename + ' to ' + outFolder)
                else:
                    print('Already downloaded %s' % out_filename)

                if not os.path.isfile(out_filename):
                    sys.exit('No NCEP wind files found for requested time')

        ######################################################
        # Open downloaded grib file with a(ny) Nansat mapper
        ######################################################
        w = Nansat(out_filename)
        self._copy_from_dataset(w.vrt.dataset)

        return

    #def __del__(self):
    #    super(Mapper,self).__del__()
    #    # Delete the downloaded ncep data
    #    try:
    #        os.unlink(self.dataset.GetFileList()[0])
    #    except:
    #        pass
