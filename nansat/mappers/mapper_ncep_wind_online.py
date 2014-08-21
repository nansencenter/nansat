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
from nansat.tools import WrongMapperError

# Place to store downloads - this can be changed via the "outFolder" argument
# to Mapper.__init__
downloads = os.path.join(os.path.expanduser('~'), 'ncep_gfs_downloads')

class Mapper(VRT, object):
    ''' VRT with mapping of WKV for NCEP GFS '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, outFolder=downloads, **kwargs):
        ''' Create NCEP VRT '''

        if not os.path.exists(outFolder):
            os.mkdir(outFolder)

        ##############
        # Get time
        ##############
        keywordBase = 'ncep_wind_online'
        if fileName[0:len(keywordBase)] != keywordBase:
            raise WrongMapperError(__file__, "Wrong mapper")

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
        outFileName = os.path.join(outFolder, 'ncep_gfs_' + nearestModelRun.strftime(
                    '%Y%m%d_%HH_') + '%.2d' % forecastHour + '.10m.uv.grib2' )
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
                outFileName = os.path.join(outFolder, fileName)
                print 'Downloading ' + url + fileName

                # Download subset of grib file
                mapperDir = os.path.dirname(os.path.abspath(__file__))
                get_inv = os.path.join(mapperDir, 'get_inv.pl')
                if not os.path.isfile(get_inv):
                    raise IOError('%s: File not found' %get_inv)
                get_grib = os.path.join(mapperDir, 'get_grib.pl')
                if not os.path.isfile(get_grib):
                    raise IOError('%s: File not found' %get_grib)
                if not os.path.isfile(outFileName):
                    command = get_inv + ' ' + url + baseName \
                          + '.inv | egrep "(:UGRD:10 m |:VGRD:10 m )" | ' \
                          + get_grib + ' ' + url + fileName + ' ' + outFileName
                    os.system(command)
                    if os.path.isfile(outFileName):
                        print 'Downloaded ' + fileName + ' to ' + outFolder
                else:
                    print 'Already downloaded %s' %outFileName
                if not os.path.isfile(outFileName):
                    sys.exit('No NCEP wind files found for requested time')


        ######################################################
        # Open downloaded grib file with a(ny) Nansat mapper
        ######################################################
        from nansat.nansat import Nansat
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
