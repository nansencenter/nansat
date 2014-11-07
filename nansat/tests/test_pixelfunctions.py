import unittest


class TestPixelFunctions(unittest.TestCase):
    def test_import_pixel_functions(self):
        try:
            from nansat._pixfun import registerPixelFunctions
        except:
            self.fail('Cannot import pixel functions')
