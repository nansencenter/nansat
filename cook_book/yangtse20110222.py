#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:		yangtse20110222.py
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	06.11.2014
# Last modified:06.11.2014 13:00
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import sys, os
home = os.path.expanduser("~")

import ipdb
ipdb.set_trace()

import numpy as np
import matplotlib.pyplot as plt

from nansat.nansatmap import Nansatmap
from nansat.nansat import Nansat, Domain


n = Nansat(os.path.join(home,
        'conferences/ESABigData2014/demo_data/MER_RR__1PRACR20110222_020119_000025973099_00391_46957_0000.N1'))

d = Domain(4326, '-lle 118 28 132 40 -ts 1000 800')

n.reproject(d)

w = Nansat(os.path.join(home,
        'conferences/ESABigData2014/demo_data/gfs20110222/gfs.t00z.master.grbf03'))

w.reproject(d)

L_560 = n['L_560']

L_560[L_560>90] = np.nan

u = w['U']

v = w['V']

nMap = Nansatmap(n, resolution='f')

nMap.pcolormesh(L_560, vmin=20, vmax=90)

nMap.quiver(u,v,quivectors=20)

nMap.add_colorbar()

nMap.drawgrid()

nMap.quiver(u,v,quivectors=20)

plt.suptitle('TOA radiance from MERIS image and wind direction, 2011.02.22')

nMap.save('20110222_L_560.png')
