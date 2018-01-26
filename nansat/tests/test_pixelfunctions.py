import sys
import importlib
import unittest
pixfun_module_name = 'nansat._pixfun_py{0}'.format(sys.version_info[0])

class TestPixelFunctions(unittest.TestCase):
    def test_import_pixel_functions(self):
        try:
            pixfun = importlib.import_module(pixfun_module_name)
            pixfun.registerPixelFunctions()
        except ImportError:
            self.fail('Cannot import pixel functions')
