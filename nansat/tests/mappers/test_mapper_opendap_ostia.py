import unittest
from nansat.mappers.mapper_opendap_ostia import Mapper


class MyWaveOpenDAPTests(unittest.TestCase):

    def setUp(self):
        self.src = 'https://podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/L4/GLOB/' \
                   'UKMO/OSTIA/2016/006/20160106-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc.bz2'

    def test_get_date(self):
        res = Mapper.get_date(self.src)
        self.assertIsInstance(res, str)
        self.assertEqual(res, '2016-01-06T00:00:00Z')