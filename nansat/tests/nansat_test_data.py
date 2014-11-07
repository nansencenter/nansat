#------------------------------------------------------------------------------
# Name:         nansat_test_data.py
# Purpose:      Get/Create directories to store test data and results of tests
#
# Author:       Anton Korosov
#
# Created:      29.09.2014
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
#------------------------------------------------------------------------------
import os

tests_path = os.path.dirname(os.path.abspath(__file__))
test_data_path = os.path.join(tests_path, 'data')
tmp_data_path = os.path.join(tests_path, 'data', 'test_data')

if not os.path.exists(tmp_data_path):
    os.mkdir(tmp_data_path)
