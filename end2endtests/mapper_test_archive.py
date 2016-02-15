# ------------------------------------------------------------------------------
# Name:         mapper_test_archive.py
# Purpose:      To discover test data for the end2endtests
#
# Author:       Anton Korosov, Morten Wergeland Hansen, Asuka Yamakawa
# Modified:     Morten Wergeland Hansen, Aleksander Vines
#
# Created:      2014-06-18
# Last modified:2015-12-28 13:36
# Copyright:    (c) NERSC
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------------------
import os
import glob


class DataForTestingMappers(object):
    def __init__(self):
        ''' Find test files and corresponding mapper names '''
        existingTestFiles = self.find_existing_files()
        filesAndMappers = self.identify_mappers(existingTestFiles)

        # FIXME: Why are we operating with two different names on this var?
        self.mapperData = filesAndMappers

    def find_existing_files(self):
        ''' Find all files for testing inside MAPPER_TEST_DATA_DIR'''
        testFiles = []

        testDataEnv = os.getenv('MAPPER_TEST_DATA_DIR')
        # FIXME: Give a warning if this is None?
        if testDataEnv is not None:
            testDataDirs = testDataEnv.split(':')
            for testDataDir in testDataDirs:
                if os.path.isdir(testDataDir):
                    testFiles += glob.glob(os.path.join(testDataDir, '*', '*'))

        testFiles = [f for f in testFiles if self.readable(f)]

        return testFiles

    def identify_mappers(self, testFiles):
        ''' Get the name of the mapper from the sub-directory name '''

        mapperNames = [os.path.split(os.path.split(testFile)[0])[1]
                       for testFile in testFiles]
        return zip(testFiles, mapperNames)

    def readable(self, testFile):
        ''' Test if file is readable at OS level '''
        if not os.path.exists(testFile):
            return False
        if not os.access(testFile, os.R_OK):
            return False
        if os.stat(testFile).st_size == 0:
            return False

        return True
