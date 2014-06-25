#-------------------------------------------------------------------------------
# Name:         test_nansat.py
# Purpose:      To test nansat
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	18.06.2014
# Last modified:18.06.2014 14:53
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import os, sys, glob
from types import *
import unittest

from nansat import *

class NansatTest(unittest.TestCase):

    def test_mapper_imports(self):
        for folder in sys.path:
            for mapper in glob.glob(folder + '/mapper_*.py'):
                mm = os.path.basename(mapper.replace('.py',''))
                assert type(__import__(mm)) is ModuleType
