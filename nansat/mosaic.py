# Name:    nansat.py
# Name:  nansat.py
# Purpose: Container of Nansat class
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
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
from __future__ import absolute_import
import multiprocessing as mp
import datetime

import numpy as np
import scipy.stats as st

from nansat.nansat import Nansat


class Mosaic(Nansat):
    '''Container for mosaicing methods

    Mosaic inherits everything from Nansat
    '''

    # default parameters
    nClass = Nansat
    eResampleAlg = 0
    period = None, None
    threads = 1
    maskName = 'mask'
    doReproject = True
    bandIDs = [1]
    mapper='mosaic'

    def _set_defaults(self, idict):
        '''Check input params and set defaut values

        Look throught default parameters (in self) and given parameters (dict)
        and paste value from input if the key matches

        Parameters
        ----------
        idict : dictionary
            parameter names and values

        Modifies
        ---------
        self

        '''
        for key in idict:
            if hasattr(self, key):
                setattr(self, key, idict[key])

    def _get_layer_image(self, f):
        '''Get nansat object from the specifed file

        Open file with Nansat
        Return, if it is within a given period,

        Parameters:
        -----------
        f : string
            name of the file to be opened
        Returns:
        --------
            Nansat object or None
        '''
        # open file using Nansat or its child class
        # the line below is for debugging
        #n = self.nClass(f, logLevel=self.logger.level)
        self.logger.info('Try to open %s' % f)
        #n = self.nClass(f, logLevel=self.logger.level)
        try:
            n = self.nClass(f, logLevel=self.logger.level)
        except:
            self.logger.error('Unable to open %s' % f)
            return None

        # check if image is out of period
        self.logger.info('Try to get time from %s' % f)
        if n is not None:
            ntime = n.get_time()

            if (ntime[0] is None and any(self.period)):
                self.logger.error('%s has no time' % f)
                return None

            if (self.period[0] is not None and
                    ntime[0] < self.period[0]):
                self.logger.info('%s is taken before the period' % f)
                return None

            if (self.period[1] is not None and
                    ntime[0] > self.period[1]):
                self.logger.info('%s is taken after the period' % f)
                return None

        return n

    def _get_layer_mask(self, n):
        '''Get mask from input Nansat object

        Open files, reproject, get mask and metadata

        Parameters:
        -----------
        n : Nansat
            input object
        doReproject : boolean
            Should we reproject input files?
        maskName : string
            Name of the mask in the input file

        Returns:
        --------
        mask : Numpy array with L2-mask
        '''
        mask = 64 * np.ones(self.shape()).astype('int8')
        # add mask band [0: nodata, 1: cloud, 2: land, 64: data]
        self.logger.info('Try to get raw mask')
        try:
            mask = n[self.maskName]
        except:
            self.logger.error('Cannot get mask from %s' % n.fileName)
            n.add_band(array=mask, parameters={'name': self.maskName})
        self.logger.debug('Got raw mask - OK')

        if self.doReproject:
            # reproject image and get reprojected mask
            self.logger.debug('Try to get reprojected mask')
            n.reproject(self, eResampleAlg=self.eResampleAlg)
            try:
                mask = n[self.maskName]
            except:
                self.logger.error('Unable to get reprojected mask!')
            self.logger.debug('Get reprojected mask - OK')

        return mask

    def _get_layer(self, f):
        '''Get nansat object and mask from input file

        Parameters:
        -----------
        f : string
            input filename
        doReproject : boolean
            Should we reproject input files?
        maskName : string
            Name of the mask in the input file

        Returns:
        --------
        n : Nansat object of input file
        mask : Numpy array with array
        '''
        mask = None
        n = self._get_layer_image(f)
        if n is not None:
            mask = self._get_layer_mask(n)

        return n, mask

    def _get_cube(self, files, band):
        '''Make cube with data from one band of input files

        Open files, reproject, get band, insert into cube

        Parameter:
        ----------
        files : list of strings
            input filenames
        band : int or string
            ID of the band
        doReproject : boolean
            Should we reproject input files?
        maskName : string
            Name of the mask in the input file

        Returns:
        --------
        dataCube : Numpy 3D array with bands
        mask : Numpy array with L2-mask
        '''
        # preallocate 3D cube and mask
        self.logger.debug('Allocating 3D cube')
        dataCube = np.zeros((len(files), self.shape()[0], self.shape()[1]))
        maskMat = np.zeros((2, self.shape()[0], self.shape()[1]), 'int8')

        # for all input files
        for i, f in enumerate(files):
            self.logger.info('Processing %s' % f)

            # get image and mask
            n, mask = self._get_layer(f)
            if n is None:
                continue
            # get band from input image
            a = None
            try:
                a = n[band].astype('float32')
            except:
                self.logger.error('%s is not in %s' % (band, n.fileName))
            if a is not None:
                # mask invalid data
                a[mask <= 2] = np.nan

            # add band to the cube
            dataCube[i, :, :] = a

            # add data to mask matrix (maximum of 0, 1, 2, 64)
            maskMat[0, :, :] = mask
            maskMat[1, :, :] = maskMat.max(0)

        return dataCube, maskMat.max(0)

    def average(self, files=[], bands=[1], doReproject=True, maskName='mask',
                threads=1, **kwargs):
        '''Memory-friendly, multithreaded mosaicing(averaging) of input files

        Convert all input files into Nansat objects, reproject onto the
        Domain of the current object, get bands, from each object,
        calculate average and STD, add averaged bands (and STD) to the current
        object.

        average() tries to get band 'mask' from the input files. The mask
        should have the following coding:
            0 : nodata
            1 : clouds
            2 : land
            64 : valid pixel
        If it gets that band (which can be provided by some mappers or Nansat
        childs, e.g.  ModisL2Image) it uses it to select averagable pixels
        (i.e. where mask == 64).
        If it cannot locate the band 'mask' is assumes that all pixels are
        averagebale except for thouse out of swath after reprojection.

        average() adds bands to the object, so it works only with empty, or
        non-projected objects

        Parameters
        -----------
        files : list
            list of input files
        bands : list
            list of names/band_numbers to be processed
        doReproject : boolean, [True]
            reproject input files?
        maskName : str, ['mask']
            name of the mask in input files
        nClass : child of Nansat, [Nansat]
            This class is used to read input files
        threads : int
            number of parallel processes to use
        eResampleAlg : int, [0]
            agorithm for reprojection, see Nansat.reproject()
        period : [datetime0, datetime1]
            Start and stop datetime objects from pyhon datetime.

        '''
        # check inputs
        if len(files) == 0:
            self.logger.error('No input files given!')
            return

        # modify default values
        self.bandIDs = bands
        self.doReproject = doReproject
        self.maskName = maskName
        self.threads = threads
        self._set_defaults(kwargs)

        # get desired shape
        dstShape = self.shape()
        self.logger.debug('dstShape: %s' % str(dstShape))

        # preallocate 2D matrices:
        # sum, sum of squares, count of products and mask
        self.logger.debug('Allocating 2D matrices')
        avgMat = np.zeros((len(bands), dstShape[0], dstShape[1]))
        stdMat = np.zeros((len(bands), dstShape[0], dstShape[1]))
        cntMat = np.zeros((dstShape[0], dstShape[1]), 'float16')
        maskMat = np.zeros((2, dstShape[0], dstShape[1]), 'int8')

        # put 2D matrices into result queue (variable shared by sub-processes)
        matQueue = mp.Queue()
        matQueue.put((cntMat, maskMat, avgMat, stdMat, files[0]))

        # create task queue with file names
        fQueue = mp.JoinableQueue()

        # generate sub-processes
        procs = []
        for i in range(threads):
            procs.append(mp.Process(target=self._average_one_file,
                                    args=(fQueue, matQueue)))

        # start sub-processes
        for i in range(threads):
            procs[i].start()

        # put file names into task queue
        for f in files:
            fQueue.put(f)
        # add poison pill to task queue
        for i in range(threads):
            fQueue.put(None)

        # wait until sub-processes get all tasks from the task queue
        fQueue.join()

        # get data from result queue
        cntMat, maskMat, avgMat, stdMat, fName = matQueue.get()

        # average products
        cntMat[cntMat == 0] = np.nan
        for bi, b in enumerate(bands):
            self.logger.debug('    Averaging %s' % b)
            # get average
            avg = avgMat[bi] / cntMat
            # calculate STD
            # STD = sqrt(sum((x-M)^2)/n) = (sqrt((sum(x^2) -
            #                                2*mean(x)*sum(x) +
            #                                sum(mean(x)^2))/n))
            stdMat[bi] = np.sqrt((stdMat[bi] - 2.0 * avg * avgMat[bi] +
                                np.square(avg) * cntMat) / cntMat)
            # set mean
            avgMat[bi] = avg

        # calculate mask (max of 0, 1, 2, 4, 64)
        maskMat = maskMat.max(0)
        # if old 'valid' mask was applied in files, replace with new mask
        maskMat[maskMat == 128] = 64

        self.logger.debug('Adding bands')
        # add mask band
        self.logger.debug('    mask')
        self.add_band(array=maskMat, parameters={'name': maskName,
                                                 'long_name': 'L2-mask',
                                                 'standard_name': 'mask'})
        firstN = self._get_layer_image(fName)

        # add averaged bands with metadata
        for bi, b in enumerate(bands):
            self.logger.debug('    %s' % b)
            # get metadata of this band from the first image
            parameters = firstN.get_metadata(bandID=b)
            parameters.pop('dataType')
            parameters.pop('SourceBand')
            parameters.pop('SourceFilename')
            # add band and std with metadata
            self.add_band(array=avgMat[bi], parameters=parameters)
            parameters['name'] = parameters['name'] + '_std'
            self.add_band(array=stdMat[bi], parameters=parameters)

    def _average_one_file(self, fQueue, matQueue):
        ''' Parallel processing of one file

        In infinite loop wait for tasks in the task queue
        If the task is available, get it and proceed
        If task is None (poison pill) quit the infinite loop
        If task is filename:
            open the file
            reproject
            get data from file,
            get intermediate result from the result queue
            add data from file into the result
            put the intermedieate result back to the queue

        Parameters
        ----------
            fQueue : multiprocessing.JoinableQueue
                task queue with file names
            matQueue : multiprocessing.Queue
                result queue with cntMat, avgMat, stdMat and maskMat

        Modifies
        --------
            fQueue : get results from the task queue
            matQueue : get and put results from into the result queue
        '''
        # start infinite loop
        while True:
            # get task from the queue
            f = fQueue.get()

            if f is None:
                # if poison pill received, quit infinite loop
                fQueue.task_done()
                break

            # otherwise start processing of task
            self.logger.info('Processing %s' % f)

            dstShape = self.shape()

            # get image and mask
            self.logger.info('Open %s and get mask' % f)
            n, mask = self._get_layer(f)

            # skip processing of invalid image
            if n is None:
                self.logger.error('%s invalid file!' % f)
                fQueue.task_done()
                continue

            # create temporary matrices to store results
            cntMatTmp = np.zeros((dstShape[0], dstShape[1]), 'float16')
            cntMatTmp[mask == 64] = 1
            avgMatTmp = np.zeros((len(self.bandIDs), dstShape[0], dstShape[1]), 'float16')
            stdMatTmp = np.zeros((len(self.bandIDs), dstShape[0], dstShape[1]), 'float16')

            # add data to summation matrices
            for bi, b in enumerate(self.bandIDs):
                self.logger.info('    Adding %s to sum' % b)
                # get projected data from Nansat object
                a = None
                try:
                    a = n[b]
                except:
                    self.logger.error('%s is not in %s' % (b, n.fileName))
                if a is not None:
                    # mask invalid data
                    a[mask < 64] = 0
                    # sum of valid values and squares
                    avgMatTmp[bi] += a
                    stdMatTmp[bi] += np.square(a)
            # destroy Nansat image
            n = None

            # get intermediate results from queue
            cntMat, maskMat, avgMat, stdMat, fName  = matQueue.get()

            # add data to the counting matrix
            cntMat += cntMatTmp

            # add data to the mask matrix (maximum of 0, 1, 2, 64)
            maskMat[0, :, :] = mask
            maskMat[1, :, :] = maskMat.max(0)

            # add data to sum and square_sum matrix
            avgMat += avgMatTmp
            stdMat += stdMatTmp

            # remember file name
            fName = f

            # update intermediate results into queue
            matQueue.put((cntMat, maskMat, avgMat, stdMat, fName))

            # tell the queue that task is done
            fQueue.task_done()

    def median(self, files=[], bands=[1], doReproject=True, maskName='mask',
               **kwargs):
        '''Calculate median of input bands

        Memory and CPU greedy method. Generates 3D cube from bands of
        all input images and calculates median. Adds median bands to self

        Parameters
        -----------
        files : list
            list of input files
        bands : list
            list of names/band_numbers to be processed
        doReproject : boolean, [True]
            reproject input files?
        maskName : str, ['mask']
            name of the mask in input files
        nClass : child of Nansat, [Nansat]
            This class is used to read input files
        eResampleAlg : int, [0]
            agorithm for reprojection, see Nansat.reproject()
        period : [datetime0, datetime1]
            Start and stop datetime objects from pyhon datetime.

        '''
        # check inputs
        if len(files) == 0:
            self.logger.error('No input files given!')
            return

        # modify default values
        self.bandIDs = bands
        self.doReproject = doReproject
        self.maskName = maskName
        self._set_defaults(kwargs)

        lastN = self._get_layer_image(files[-1])
        # add medians of all bands
        for band in bands:
            bandCube, mask = self._get_cube(files, band)
            bandMedian = st.nanmedian(bandCube, axis=0)

            # get metadata of this band from the last image
            parameters = lastN.get_metadata(bandID=band)
            # add band and std with metadata
            self.add_band(array=bandMedian, parameters=parameters)

        self.add_band(array=mask, parameters={'name': 'mask'})

    def latest(self, files=[], bands=[1], doReproject=True, maskName='mask',
               **kwargs):
        '''Mosaic by adding the latest image on top without averaging

        Uses Nansat.get_time() to estimate time of each input file;
        Sorts images by aquisition time;
        Creates date_index band - with mask of coverage of each frame;
        Uses date_index to fill bands of self only with the latest data

        Parameters
        -----------
        files : list
            list of input files
        bands : list
            list of names/band_numbers to be processed
        doReproject : boolean, [True]
            reproject input files?
        maskName : str, ['mask']
            name of the mask in input files
        nClass : child of Nansat, [Nansat]
            This class is used to read input files
        eResampleAlg : int, [0]
            agorithm for reprojection, see Nansat.reproject()
        period : [datetime0, datetime1]
            Start and stop datetime objects from pyhon datetime.

        '''
        # check inputs
        if len(files) == 0:
            self.logger.error('No input files given!')
            return

        # modify default values
        self.bandIDs = bands
        self.doReproject = doReproject
        self.maskName = maskName
        self._set_defaults(kwargs)

        # collect ordinals of times of each input file
        itimes = np.zeros(len(files))
        for i in range(len(files)):
            n = self._get_layer_image(files[i])
            nstime = n.get_time()[0]
            if nstime is None:
                nstime = 693596  # 1900-01-01
            else:
                nstime = nstime.toordinal()
            itimes[i] = nstime

        # sort times
        ars = np.argsort(itimes)

        # maxIndex keeps mask of coverae of each frame
        maxIndex = np.zeros((2, self.shape()[0], self.shape()[1]))
        for i in range(len(files)):
            # open file and get mask
            n, mask = self._get_layer(files[ars[i]])
            # fill matrix with serial number of the file
            maskIndex = (np.zeros(mask.shape) + i + 1).astype('uint16')
            # erase non-valid values
            maskIndex[mask != 64] = 0
            # first layer of maxIndex keeps serial number of this file
            maxIndex[0, :, :] = maskIndex
            # second layer of maxIndex keeps maximum serial number
            # or serial number of the latest image
            maxIndex[1, :, :] = maxIndex.max(0)
        maxIndex = maxIndex.max(0)

        # preallocate 2D matrices for mosaiced data and mask
        self.logger.debug('Allocating 2D matrices')
        avgMat = {}
        for b in bands:
            avgMat[b] = np.zeros((maxIndex.shape[0], maxIndex.shape[1]))
        maskMat = np.zeros((maxIndex.shape[0], maxIndex.shape[1]))

        for i in range(len(files)):
            f = files[ars[i]]
            self.logger.info('Processing %s' % f)

            # get image and mask
            n, mask = self._get_layer(f)
            if n is None:
                continue
            # insert mask into result only for pixels masked
            # by the serial number of the input file
            maskMat[maxIndex == (i + 1)] = mask[maxIndex == (i + 1)]

            # insert data into mosaic matrix
            for b in bands:
                self.logger.debug('    Inserting %s to latest' % b)
                # get projected data from Nansat object
                a = None
                try:
                    a = n[b].astype('float32')
                except:
                    self.logger.error('%s is not in %s' % (b, n.fileName))
                if a is not None:
                    # insert data into result only for pixels masked
                    # by the serial number of the input file
                    avgMat[b][maxIndex == (i + 1)] = a[maxIndex == (i + 1)]

            # destroy input nansat
            n = None
        # keep last image opened
        lastN = self._get_layer_image(f)

        self.logger.debug('Adding bands')
        # add mask band
        self.logger.debug('    mask')
        self.add_band(array=maskMat, parameters={'name': maskName,
                                                 'long_name': 'L2-mask',
                                                 'standard_name': 'mask'})
        # add mosaiced bands with metadata
        for b in bands:
            self.logger.debug('    %s' % b)

            # get metadata of this band from the last image
            parameters = lastN.get_metadata(bandID=b)
            # add band with metadata
            self.add_band(array=avgMat[b], parameters=parameters)

        # compose list of dates of input images
        timeString = ''
        dt = datetime.datetime(1, 1, 1)
        for i in range(len(itimes)):
            timeString += dt.fromordinal(int(itimes[ars[i]])).strftime('%Y-%m-%dZ%H:%M ')
        # add band with mask of coverage of each frame
        self.add_band(array=maxIndex, parameters={'name': 'date_index',
                                                  'values': timeString})
