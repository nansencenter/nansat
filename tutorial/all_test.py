#!/usr/bin/env python
# Name:    all_test.py
# Purpose: Test all class of Nansat
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      27.08.2013
# Copyright:    (c) NERSC 2011 - 2013
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# import all test files of Nansat
import domain_test
import nansat_test
import figure_test
import nansatmap_test
import nansatshape_test
import mosaic_test
import pointbrowser_test

def main():
    # Test domain_test
    domain_test.main()

    # Test nansat_test
    nansat_test.main()

    # Test figure_test
    figure_test.main()

    # Test nansatmap_test
    nansatmap_test.main()

    # Test nansatshape_test
    nansatshape_test.main()

    # Test mosaic_test
    mosaic_test.main()

    # Test pointbrowser_test
    pointbrowser_test.main()



