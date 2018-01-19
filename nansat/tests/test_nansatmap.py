#------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      Test the Nansat class
#
# Author:       Morten Wergeland Hansen, Anton Korosov, Asuka Yamakawa
# Modified: Morten Wergeland Hansen
#
# Created:      18.06.2014
# Last modified:18.11.2014 11:48
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------
import unittest
import os

import matplotlib.pyplot as plt
from nansat import Nansat
try:
    from nansat import Nansatmap
except ImportError:
    NANSATMAP_EXISTS = False
else:
    NANSATMAP_EXISTS = True
    

import nansat_test_data as ntd

IS_CONDA = 'conda' in os.environ['PATH']


class NansatmapTest(unittest.TestCase):
    def setUp(self):
        self.test_file_stere = os.path.join(ntd.test_data_path, 'stere.tif')
        plt.switch_backend('Agg')

        if not os.path.exists(self.test_file_stere):
            raise ValueError('No test data available')

    @unittest.skipUnless(NANSATMAP_EXISTS, 'Nansatmap is required')
    def test_create_map(self):
        ''' should simply create a Nansatmap instance '''
        n = Nansat(self.test_file_stere, logLevel=40)
        nmap = Nansatmap(n)

        self.assertEqual(type(nmap), Nansatmap)

    @unittest.skipUnless(NANSATMAP_EXISTS, 'Nansatmap is required')
    def test_imshow(self):
        ''' Should use Nansatmap.imshow '''
        n = Nansat(self.test_file_stere, logLevel=40)
        b1 = n[1]
        nmap = Nansatmap(n)
        nmap.imshow(b1, cmap='ak01')
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansatmap_imshow.png')
        nmap.save(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    @unittest.skipUnless(NANSATMAP_EXISTS, 'Nansatmap is required')
    def test_imshow_random(self):
        ''' Should use Nansatmap.imshow '''
        n = Nansat(self.test_file_stere, logLevel=40)
        b1 = n[1]
        nmap = Nansatmap(n)
        nmap.imshow(b1/5, cmap='random')
        nmap.add_colorbar()
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansatmap_imshow_random.png')
        nmap.save(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

    @unittest.skipUnless(NANSATMAP_EXISTS, 'Nansatmap is required')
    def test_pcolormesh(self):
        ''' Should use Nansatmap.pcolormesh '''
        n = Nansat(self.test_file_stere, logLevel=40)
        b1 = n[1]
        nmap = Nansatmap(n)
        nmap.pcolormesh(b1)
        tmpfilename = os.path.join(ntd.tmp_data_path, 'nansatmap_pcolormesh.png')
        nmap.save(tmpfilename)

        self.assertTrue(os.path.exists(tmpfilename))

if __name__ == "__main__":
    unittest.main()
