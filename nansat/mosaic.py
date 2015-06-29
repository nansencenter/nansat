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
import warnings
import multiprocessing as mp
import ctypes
import datetime

import numpy as np
import scipy.stats as st

from nansat.nansat import Nansat

## shared arrays for mean, squared mean and count
sharedArray = None
domain = None

def mparray2ndarray(sharedArray, shape, dtype='float32'):
    ''' convert shared multiprocessing Array to numpy ndarray '''
    # get access to shared array and convert to numpy ndarray
    sharedNDArray = np.frombuffer(sharedArray.get_obj(), dtype=dtype)
    # change shape to match bands
    sharedNDArray.shape = shape

    return sharedNDArray

def sumup(layer):
    ''' Sum up bands from input images in multiple threads'''
    global sharedArray
    global domain

    # get nansat from the input Layer
    layer.make_nansat_object(domain)
    # if not in the period, quit
    if not layer.within_period():
        return 1
    # get mask
    mask = layer.get_mask_array()
    # get arrays with data
    bandArrays = [layer.n[band] for band in layer.bands]
    bandArrays = np.array(bandArrays)
    finiteMask = np.isfinite(bandArrays.sum(axis=0))
    # get metadata
    bandMetadata = [layer.n.get_metadata(bandID=band) for band in layer.bands]

    with sharedArray.get_lock(): # synchronize access
        sharedNDArray = mparray2ndarray(sharedArray,
                                        (2+len(layer.bands)*2,
                                         mask.shape[0],
                                         mask.shape[1]),
                                        'float32')
        gpi = finiteMask * (mask == 64)

        # update counter
        sharedNDArray[0][gpi] += 1

        # update mask with max
        sharedNDArray[1] = np.max([sharedNDArray[1], mask], axis=0)

        # update sum for each band
        for i, bandArray in enumerate(bandArrays):
            sharedNDArray[2+i][gpi] += bandArray[gpi]

        # update squared sum for each band
        for i, bandArray in enumerate(bandArrays):
            sharedNDArray[2+len(layer.bands)+i][gpi] += bandArray[gpi]

    # release layer
    layer = None
    return bandMetadata

class Layer:
    ''' Small class to get mask and arrays from many bands '''
    def __init__(self, fileName, bands=[1],
                        opener=Nansat, maskName='mask',
                        doReproject=True, eResampleAlg=0,
                        period=(None,None),
                        logLevel=30):
        # Set parameters of processing
        self.fileName     = fileName
        self.bands        = bands
        self.opener       = opener
        self.maskName     = maskName
        self.doReproject  = doReproject
        self.period       = period
        self.eResampleAlg = eResampleAlg
        self.logLevel     = logLevel

    def make_nansat_object(self, domain):
        # Open self.fileName with self.opener
        print 'Layer', self.fileName, self.logLevel
        self.n = self.opener(self.fileName, logLevel=self.logLevel)
        if self.doReproject:
            self.n.reproject(domain, eResampleAlg=self.eResampleAlg)

    def within_period(self):
        ''' Test if given file is within period of time '''
        withinPeriod = True
        ntime = self.n.get_time()
        if (ntime[0] is None and any(self.period)):
            withinPeriod = False

        if (self.period[0] is not None and ntime[0] < self.period[0]):
            withinPeriod = False

        if (self.period[1] is not None and ntime[0] > self.period[1]):
            withinPeriod = False

        return withinPeriod

    def get_mask_array(self):
        ''' Get array with mask values '''
        if self.n.has_band(self.maskName):
            mask = self.n[self.maskName]
        elif self.doReproject:
            mask = self.n['swathmask'] * 64
        else:
            mask = np.ones(self.n.shape()) * 64

        return mask


class Mosaic(Nansat):
    '''Container for mosaicing methods

    Mosaic inherits everything from Nansat
    '''

    def average(self, files=[], bands=[1], doReproject=True, maskName='mask',
                opener=Nansat, threads=1, eResampleAlg=0, period=(None,None)):
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
        opener : child of Nansat, [Nansat]
            This class is used to read input files
        threads : int
            number of parallel processes to use
        eResampleAlg : int, [0]
            agorithm for reprojection, see Nansat.reproject()
        period : [datetime0, datetime1]
            Start and stop datetime objects from pyhon datetime.

        '''
        # shared array for multiple threads
        global sharedArray
        global domain

        # check inputs
        if len(files) == 0:
            self.logger.error('No input files given!')
            return

        # get desired shape
        dstShape = self.shape()
        # preallocate shared mem array
        sharedArray = mp.Array(ctypes.c_float, [0]*(2+len(bands)+len(bands)) * dstShape[0] * dstShape[1])

        # create list of layers
        domain = Nansat(domain=self)
        layers = [Layer(ifile, bands, opener, maskName, doReproject,
                        eResampleAlg, period, self.logger.level)
                        for ifile in files]

        # test in debug
        #sumup(layers[0])

        # prepare pool of processors
        pool = mp.Pool(threads)

        # run reprojection and summing up
        metadata = pool.map(sumup, layers)

        # get band metadata from the first valid file
        for bandsMeta in metadata:
            if type(bandsMeta) is list:
                break

        # average products
        sharedNDArray = mparray2ndarray(sharedArray,
                                     (2+len(bands)*2, dstShape[0], dstShape[1]),
                                     'float32')

        # cleanup
        pool.terminate()
        pool = None
        layers = None
        metadata = None
        sharedArray = None

        cntMat = sharedNDArray[0]
        maskMat = sharedNDArray[1]
        avgMat = sharedNDArray[2:2+len(bands)]
        stdMat = sharedNDArray[2+len(bands):]

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

        self.logger.debug('Adding bands')
        # add mask band
        self.logger.debug('    mask')
        self.add_band(array=maskMat, parameters={'name': maskName,
                                                 'long_name': 'L2-mask',
                                                 'standard_name': 'mask'})

        # add averaged bands with metadata
        for bi, b in enumerate(bands):
            self.logger.debug('    %s' % b)
            # add band and std with metadata
            self.add_band(array=avgMat[bi], parameters=bandsMeta[bi])
            bandsMeta[bi]['name'] = bandsMeta[bi]['name'] + '_std'
            self.add_band(array=stdMat[bi], parameters=bandsMeta[bi])


    def _get_cube(self, files, band, doReproject, maskName, opener,
                                                    eResampleAlg,
                                                    period,
                                                    vmin=-np.inf,
                                                    vmax=np.inf):
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
        opener : class
            Nansat or any Nansat child to open input image
        eResampleAlg : int
            parameter for Nansat.reproject()
        period : tuple
            valid (start_date, end_date) or (None, None)

        Returns:
        --------
            dataCube : Numpy 3D array with bands
            mask : Numpy array with L2-mask
            metadata : dict with band metadata
        '''
        # preallocate 3D cube and mask
        self.logger.debug('Allocating 3D cube')
        dataCube = np.zeros((len(files), self.shape()[0], self.shape()[1]))
        maskMat = np.zeros((2, self.shape()[0], self.shape()[1]), 'int8')

        # for all input files
        for i, f in enumerate(files):
            self.logger.info('Processing %s' % f)
            layer = Layer(f, [band], opener, maskName, doReproject,
                            eResampleAlg, period, logLevel=self.logger.level)
            # get nansat from the input Layer
            layer.make_nansat_object(domain)

            # if not in the period, quit
            if not layer.within_period():
                continue
            # get mask
            mask = layer.get_mask_array()
            # get arrays with data
            bandArray = layer.n[band].astype('float32')
            # remove invalid data
            bandArray[mask < 64] = np.nan
            bandArray[bandArray < vmin] = np.nan
            bandArray[bandArray > vmax] = np.nan

            # get metadata
            bandMetadata = layer.n.get_metadata(bandID=band)

            # add band to the cube
            dataCube[i, :, :] = bandArray

            # add data to mask matrix (maximum of 0, 1, 2, 64)
            maskMat[0, :, :] = mask
            maskMat[1, :, :] = maskMat.max(0)

        return dataCube, maskMat.max(0), bandMetadata

    def median(self, files=[], bands=[1], doReproject=True, maskName='mask',
                opener=Nansat, eResampleAlg=0, period=(None,None),
                vmin=-np.inf, vmax=np.inf):
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


        # add medians of all bands
        for band in bands:
            cube, mask, metadata = self._get_cube(files, band,
                                                    doReproject,
                                                    maskName,
                                                    opener,
                                                    eResampleAlg,
                                                    period, vmin, vmax)
            median = st.nanmedian(cube, axis=0)

            # add band and std with metadata
            self.add_band(array=median, parameters=metadata)

        self.add_band(array=mask, parameters={'name': 'mask'})

"""
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
            timeString += (dt.fromordinal(int(itimes[ars[i]])).
                           strftime('%Y-%m-%dZ%H:%M '))
        # add band with mask of coverage of each frame
        self.add_band(array=maxIndex, parameters={'name': 'date_index',
                                                  'values': timeString})

"""
