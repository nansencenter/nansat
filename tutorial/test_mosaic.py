#!/usr/bin/env python
# Name:    mosaic_test.py
# Purpose: Tutorial for nansat.Mosaic
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      23.08.2013
# Copyright:    (c) NERSC 2011 - 2013
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

import os
import datetime as dt
from nansat import Domain, Mosaic
# input and output file names
from testio import testio

def main():
    iPath, oPath = testio()
    iFileNames = [os.path.join(iPath, 'gcps.tif'),
                  os.path.join(iPath, 'stere.tif')]
    print 'Input files: ', iFileNames
    oFileName = os.path.join(oPath, 'output_mosaic_')
    print 'Output file:', oFileName
    
    ''' Mosaic class includes mosaicing methods
    
    mosaicing methods:
    1. average
    2. median
    3. latest (add the latest image on top)
    
    '''
    
    # Create target domain
    domain = Domain(4326, '-lle 27 70 31 72 -ts 700 650')
    
    
    print '\n\n\n ============== AVERAGE =============\n\n\n'
    # A. Perform averaging of several files
    # 1. Create destination Nansat object with desired projection
    nMosaic = Mosaic(domain=domain, logLevel=10)
    # 2. Perfom averaging
    nMosaic.average(iFileNames, bands=['L_645', 'L_555', 'L_469'])
    # 3. Get mask of valid pixels
    mask = nMosaic['mask']
    # 4. Output averaged data using the mask
    nMosaic.write_figure(fileName=oFileName + '.png',
                         bands=['L_645', 'L_555', 'L_469'],
                         clim='hist',
                         mask_array=mask,
                         mask_lut={0: [128,128,128]})
    # 5. Get values of standard deviation from averaging of input files
    L_469_std = nMosaic['L_469_std']
    
    print '\n\n\n ============== MEDIAN =============\n\n\n'
    # B. calculate median from the first band (very slow thus comented)
    nMosaic.median(iFileNames)
    
    print '\n\n\n ============== LATEST =============\n\n\n'
    # C. fill the result with the latest image
    nMosaic.latest(iFileNames)
    
    # D. Average only files that fall within given period
    # create new Mosaic
    nMosaic = Mosaic(domain=domain, logLevel=10)
    # define period
    period = [dt.datetime(2011, 8, 15, 0, 0), dt.datetime(2011, 8, 15, 23, 59)]
    # run averaging of only files within period.
    # Only gcps.tif will be averaged since stere.tif has no time information
    print '\n\n\n ============== AVERAGE WITHIN PERIOD =============\n\n\n'
    nMosaic.average(iFileNames, period=period)
    
    print '\n*** mosaic_test completed successfully. Output files are found here:' + oFileName
    
    
if __name__ == "__main__":
    main()
