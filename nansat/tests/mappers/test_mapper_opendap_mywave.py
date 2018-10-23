import unittest
from nansat.mappers.mapper_opendap_mywave import Mapper


class MyWaveOpenDAPTests(unittest.TestCase):

    def setUp(self):
        self.src = 'http://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive' \
                   '/2017/10/29/MyWave_wam4_WAVE_20171029T18Z.nc'

    def test_get_date(self):
        res = Mapper.get_date(self.src)
        self.assertIsInstance(res, str)
        self.assertEqual(res, '2017-10-29T18:00:00Z')
