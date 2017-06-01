import unittest

from nansat.nansat import Nansat
from nansat.domain import Domain
from nansat.nsr import NSR
try:
    from sardoppler.sardoppler import Doppler
except Exception as e:
    print e.message

class TestOpenIssues(unittest.TestCase):

    def test_issue_189(self):
        fn = '/mnt/10.11.12.232/sat_downloads_asar/level-0/2010-01/descending/VV/gsar_rvl/RVL_ASA_WS_20100110211812087.gsar'
        n = Doppler(fn)
        xlon, xlat = n.get_corners()
        d = Domain(NSR(3857),
                   '-lle %f %f %f %f -tr 1000 1000' % (
                    xlon.min(), xlat.min(), xlon.max(), xlat.max()))
        n.reproject(d, eResampleAlg=1, tps=True)
        inci = n['incidence_angle']

    def test_issue_193(self):
        fn = [
            '/vagrant/shared/test_data/cmems/GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS-x10-X30-y55-Y73-201705181200-201705271200.nc',
            '/vagrant/shared/test_data/cmems/ARC-METNO-ARC-TOPAZ4_2_PHYS-FOR-TDS-x10-X30-y55-Y73-20170518-20170526.nc',
            '/vagrant/shared/test_data/cmems/GLOBAL_ANALYSIS_FORECAST_BIO_001_014-TDS-x-180-X179.5-y-89-Y90-20170520-20170527.nc',
        ]
        for f in fn:
            n = Nansat(f)
            self.assertTrue(n.get_metadata().has_key('time_coverage_start'))
            self.assertTrue(n.get_metadata().has_key('time_coverage_end'))
            self.assertTrue(n.get_metadata().has_key('instrument'))
            self.assertTrue(n.get_metadata().has_key('platform'))
            self.assertEqual(n.mapper, 'cmems')
