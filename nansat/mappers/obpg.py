# Name:        obpg
# Purpose:     Base class for mapping for L2 data from the OBPG web-site
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from nansat.vrt import VRT


class OBPGL2BaseClass(VRT):
    ''' Base Class for Mappers for SeaWIFS/MODIS/MERIS/VIIRS L2 data from OBPG
    '''

    titles = ['HMODISA Level-2 Data',
              'MODISA Level-2 Data',
              'HMODIST Level-2 Data',
              'MERIS Level-2 Data',
              'GOCI Level-2 Data',
              'VIIRSN Level-2 Data',
              'SeaWiFS Level-2 Data']
