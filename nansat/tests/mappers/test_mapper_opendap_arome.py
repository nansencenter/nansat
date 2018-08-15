import unittest
from nansat.mappers.mapper_opendap_arome import Mapper
from nansat import Nansat


class AROMEOpenDAPTests(unittest.TestCase):

    def test_get_date(self):
        res = Mapper.get_date('https://thredds.met.no/thredds/dodsC/aromearcticarchive/2017/'
                              '10/30/arome_arctic_full_2_5km_20171030T21Z.nc')
        self.assertIsInstance(res, str)
        self.assertEqual(res, '2017-10-30T21:00Z')
