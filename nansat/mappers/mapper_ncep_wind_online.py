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

import os
import sys
from datetime import datetime, timedelta
import inspect

from nansat.vrt import VRT

# Could be input as argument, or defined as environment variable
outFolder = './'

class Mapper(VRT, object):
    ''' VRT with mapping of WKV for NCEP GFS '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        ''' Create NCEP VRT '''

        ##############
        # Get time
        ##############
        keywordBase = 'ncep_wind_online'
        if fileName[0:len(keywordBase)] != keywordBase:
            raise AttributeError("Wrong mapper")

        timestr = fileName[len(keywordBase)+1::]
        time = datetime.strptime(timestr, '%Y%m%d%H%M')
        print time

        ########################################
        # Find and download online grib file
        ########################################
        # Find closest 6 hourly modelrun and forecast hour
        modelRunHour = round((time.hour + time.minute/60.)/6)*6
        nearestModelRun = datetime(time.year, time.month, time.day) \
            + timedelta(hours=modelRunHour)
        forecastHour = (time - nearestModelRun).total_seconds()/3600.
        if modelRunHour == 24:
            modelRunHour = 0
        if forecastHour < 1.5:
            forecastHour = 0
        else:
            forecastHour = 3

        #########################################################
        # Try first to get NRT data from
        # ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/
        # - avaliable approximately the latest month
        #########################################################
        url = 'ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/' \
                + 'gfs.' + nearestModelRun.strftime('%Y%m%d') \
                + '%.2d' % modelRunHour \
                + '/gfs.t' + '%.2d' % modelRunHour + 'z.master.grbf' \
                + '%.2d' % forecastHour + '.10m.uv.grib2'
        outFileName = outFolder + 'ncep_gfs_' + nearestModelRun.strftime(
                    '%Y%m%d_%HH_') + '%.2d' % forecastHour + '.10m.uv.grib2'
        if os.path.exists(outFileName):
            print 'NCEP wind is already downloaded: ' + outFileName
        else:
            os.system('curl -so ' + outFileName + ' ' + url)
            if os.path.exists(outFileName):
                print 'Downloaded ' + outFileName
            else:
                print 'NRT GRIB file not available: ' + url
                #########################################################
                # If NRT file not available, search in long term archive
                #########################################################
                url = 'http://nomads.ncdc.noaa.gov/data/gfs4/' + \
                    nearestModelRun.strftime('%Y%m/%Y%m%d/')
                baseName = 'gfs_4_' + nearestModelRun.strftime('%Y%m%d_') \
                        + nearestModelRun.strftime('%H%M_') \
                        + '%.3d' % forecastHour
                fileName = baseName + '.grb2'
                outFileName = outFolder + fileName
                print 'Downloading ' + url + fileName

                # Download subset of grib file
                mapperDir = os.path.dirname(os.path.abspath(
                                inspect.getfile(inspect.currentframe())))
                get_inv = mapperDir + '/get_inv.pl '
                get_grib = mapperDir + '/get_grib.pl '
                if not os.path.exists(fileName):
                    command = get_inv + url + baseName \
                          + '.inv | egrep "(:UGRD:10 m |:VGRD:10 m )" | ' \
                          + get_grib + url + fileName + ' ' + outFileName
                    os.system(command)
                else:
                    print 'Already downloaded'
                if not os.path.exists(outFileName):
                    sys.exit('No NCEP wind files found for requested time')


        ######################################################
        # Open downloaded grib file with a(ny) Nansat mapper
        ######################################################
        # baaaad solution!
        from nansat import Nansat
        w = Nansat(outFileName)
        VRT.__init__(self, vrtDataset=w.vrt.dataset)

        return

    #def __del__(self):
    #    super(Mapper,self).__del__()
    #    # Delete the downloaded ncep data
    #    try:
    #        os.unlink(self.dataset.GetFileList()[0])
    #    except:
    #        pass
