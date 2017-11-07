import unittest

from nansat.nansat import Nansat
from nansat.domain import Domain
from nansat.nsr import NSR

doppler_installed = True
try:
    from sardoppler.sardoppler import Doppler
except Exception as e:
    print e.message
    doppler_installed = False

class TestOpenIssues(unittest.TestCase):

    def test_issue_189(self):
        fn = '/mnt/10.11.12.232/sat_downloads_asar/level-0/2010-01/descending/VV/gsar_rvl/RVL_ASA_WS_20100110211812087.gsar'
        if doppler_installed:
            n = Doppler(fn)
            xlon, xlat = n.get_corners()
            d = Domain(NSR(3857),
                    '-lle %f %f %f %f -tr 1000 1000' % (
                        xlon.min(), xlat.min(), xlon.max(), xlat.max()))
            n.reproject(d, eResampleAlg=1, tps=True)
            inci = n['incidence_angle']

