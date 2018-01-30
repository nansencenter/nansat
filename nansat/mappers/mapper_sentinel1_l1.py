#------------------------------------------------------------------------------
# Name:     mapper_sentinel1_l1.py
# Purpose:
#
# Author:       Morten Wergeland Hansen, Anton Korosov
#
# Created:  12.09.2014
# Copyright:    (c) NERSC
# License: GPL V3
#------------------------------------------------------------------------------
import warnings

import os
import glob
import zipfile
import numpy as np
from dateutil.parser import parse
import xml.etree.ElementTree as ET

try:
    import scipy
except:
    IMPORT_SCIPY = False
else:
    IMPORT_SCIPY = True

import json
import pythesint as pti

from nansat.vrt import VRT
from nansat.tools import gdal, initial_bearing
from nansat.exceptions import WrongMapperError, NansatReadError
from nansat.nsr import NSR
from nansat.node import Node


class Mapper(VRT):
    '''
        Create VRT with mapping of Sentinel-1 (A and B) stripmap mode (S1A_SM)
    '''

    def __init__(self, filename, gdalDataset, gdalMetadata,
                 manifestonly=False, **kwargs):

        if not os.path.split(filename.rstrip('/'))[1][:3] in ['S1A', 'S1B']:
            raise WrongMapperError('%s: Not Sentinel 1A or 1B' %filename)

        if zipfile.is_zipfile(filename):
            zz = zipfile.PyZipFile(filename)
            # Assuming the file names are consistent, the polarization
            # dependent data should be sorted equally such that we can use the
            # same indices consistently for all the following lists
            # THIS IS NOT THE CASE...
            mdsFiles = ['/vsizip/%s/%s' % (filename, fn)
                        for fn in zz.namelist() if 'measurement/s1' in fn]
            calFiles = ['/vsizip/%s/%s' % (filename, fn)
                        for fn in zz.namelist()
                        if 'annotation/calibration/calibration-s1' in fn]
            noiseFiles = ['/vsizip/%s/%s' % (filename, fn)
                          for fn in zz.namelist()
                          if 'annotation/calibration/noise-s1' in fn]
            annotationFiles = ['/vsizip/%s/%s' % (filename, fn)
                               for fn in zz.namelist()
                               if 'annotation/s1' in fn]
            manifestFile = ['/vsizip/%s/%s' % (filename, fn)
                            for fn in zz.namelist()
                            if 'manifest.safe' in fn]
            zz.close()
        else:
            mdsFiles = glob.glob('%s/measurement/s1*' % filename)
            calFiles = glob.glob('%s/annotation/calibration/calibration-s1*'
                                 % filename)
            noiseFiles = glob.glob('%s/annotation/calibration/noise-s1*'
                                   % filename)
            annotationFiles = glob.glob('%s/annotation/s1*'
                                        % filename)
            manifestFile = glob.glob('%s/manifest.safe' % filename)

        if (not mdsFiles or not calFiles or not noiseFiles or
                not annotationFiles or not manifestFile):
            raise WrongMapperError(filename)

        mdsDict = {}
        for ff in mdsFiles:
            mdsDict[
                os.path.splitext(os.path.basename(ff))[0].split('-')[3]] = ff

        self.calXMLDict = {}
        for ff in calFiles:
            self.calXMLDict[
                os.path.splitext(
                os.path.basename(ff))[0].split('-')[4]] = self.read_vsi(ff)

        self.noiseXMLDict = {}
        for ff in noiseFiles:
            self.noiseXMLDict[
                os.path.splitext(
                os.path.basename(ff))[0].split('-')[4]] = self.read_vsi(ff)

        self.annotationXMLDict = {}
        for ff in annotationFiles:
            self.annotationXMLDict[
                os.path.splitext(
                os.path.basename(ff))[0].split('-')[3]] = self.read_vsi(ff)

        self.manifestXML = self.read_vsi(manifestFile[0])

        missionName = {'S1A': 'SENTINEL-1A', 'S1B': 'SENTINEL-1B'}[
            os.path.split(filename.rstrip('/'))[1][:3]]

        # very fast constructor without any bands
        if manifestonly:
            self.init_from_manifest_only(self.manifestXML,
                                         self.annotationXMLDict[
                                         self.annotationXMLDict.keys()[0]],
                                         missionName)
            return

        gdalDatasets = {}
        for key in mdsDict.keys():
            # Open data files
            gdalDatasets[key] = gdal.Open(mdsDict[key])

        if not gdalDatasets:
            raise WrongMapperError('%s: No Sentinel-1 datasets found'
                    %mdsDict[key])

        # Check metadata to confirm it is Sentinel-1 L1
        metadata = gdalDatasets[list(mdsDict.keys())[0]].GetMetadata()

        if not 'TIFFTAG_IMAGEDESCRIPTION' in metadata.keys():
            raise WrongMapperError(filename)
        if (not 'Sentinel-1' in metadata['TIFFTAG_IMAGEDESCRIPTION']
                and not 'L1' in metadata['TIFFTAG_IMAGEDESCRIPTION']):
            raise WrongMapperError(filename)

        warnings.warn('Sentinel-1 level-1 mapper is not yet adapted to '
                      'complex data. In addition, the band names should be '
                      'updated for multi-swath data - '
                      'and there might be other issues.')

        if not IMPORT_SCIPY:
            raise NansatReadError('Sentinel-1 data cannot be read because scipy is not installed')

        # create empty VRT dataset with geolocation only
        for key in gdalDatasets:
            self._init_from_gdal_dataset(gdalDatasets[key], metadata=metadata)
            break

        # Read annotation, noise and calibration xml-files
        pol = {}
        it = 0
        for key in self.annotationXMLDict:
            xml = Node.create(self.annotationXMLDict[key])
            pol[key] = (xml.node('product').
                        node('adsHeader')['polarisation'].upper())
            it += 1
            if it == 1:
                # Get incidence angle
                pi = xml.node('generalAnnotation').node('productInformation')

                self.dataset.SetMetadataItem('ORBIT_DIRECTION',
                                              str(pi['pass']))
                (X, Y, lon, lat, inc, ele, numberOfSamples,
                numberOfLines) = self.read_geolocation_lut(
                                                self.annotationXMLDict[key])

                nX, nY = (Y==0).sum(), (X==0).sum()
                lon = np.array(lon).reshape(nY, nX)
                lat = np.array(lat).reshape(nY, nX)
                inc = np.array(inc).reshape(nY, nX)
                ele = np.array(ele).reshape(nY, nX)

                incVRT = VRT.from_array(inc)
                eleVRT = VRT.from_array(ele)
                incVRT = incVRT.get_resized_vrt(self.dataset.RasterXSize, self.dataset.RasterYSize, 2)
                eleVRT = eleVRT.get_resized_vrt(self.dataset.RasterXSize, self.dataset.RasterYSize, 2)
                self.band_vrts['incVRT'] = incVRT
                self.band_vrts['eleVRT'] = eleVRT

        for key in self.calXMLDict:
            calibration_LUT_VRTs, longitude, latitude = (
                self.get_LUT_VRTs(self.calXMLDict[key],
                                  'calibrationVectorList',
                                  ['sigmaNought', 'betaNought',
                                   'gamma', 'dn']
                                  ))
            self.band_vrts['LUT_sigmaNought_VRT_'+pol[key]] = (
                calibration_LUT_VRTs['sigmaNought'].
                get_resized_vrt(self.dataset.RasterXSize, self.dataset.RasterYSize, 1))
            self.band_vrts['LUT_betaNought_VRT_'+pol[key]] = (
                calibration_LUT_VRTs['betaNought'].
                get_resized_vrt(self.dataset.RasterXSize, self.dataset.RasterYSize, 1))
            self.band_vrts['LUT_gamma_VRT'] = calibration_LUT_VRTs['gamma']
            self.band_vrts['LUT_dn_VRT'] = calibration_LUT_VRTs['dn']

        for key in self.noiseXMLDict:
            noise_LUT_VRT = self.get_LUT_VRTs(self.noiseXMLDict[key],
                                              'noiseVectorList',
                                              ['noiseLut'])[0]
            self.band_vrts['LUT_noise_VRT_'+pol[key]] = (
                noise_LUT_VRT['noiseLut'].get_resized_vrt(self.dataset.RasterXSize, self.dataset.RasterYSize, 1))

        metaDict = []
        bandNumberDict = {}
        bnmax = 0
        for key in gdalDatasets.keys():
            dsPath, dsName = os.path.split(mdsDict[key])
            name = 'DN_%s' % pol[key]
            # A dictionary of band numbers is needed for the pixel function
            # bands further down. This is not the best solution. It would be
            # better to have a function in VRT that returns the number given a
            # band name. This function exists in Nansat but could perhaps be
            # moved to VRT? The existing nansat function could just call the
            # VRT one...
            bandNumberDict[name] = bnmax + 1
            bnmax = bandNumberDict[name]
            band = gdalDatasets[key].GetRasterBand(1)
            dtype = band.DataType
            metaDict.append({
                'src': {
                    'SourceFilename': mdsDict[key],
                    'SourceBand': 1,
                    'DataType': dtype,
                },
                'dst': {
                    'name': name,
                },
            })
        # add bands with metadata and corresponding values to the empty VRT
        self.create_bands(metaDict)

        '''
        Calibration should be performed as

        s0 = DN^2/sigmaNought^2,

        where sigmaNought is from e.g.
        annotation/calibration/calibration-s1a-iw-grd-hh-20140811t151231-20140811t151301-001894-001cc7-001.xml,
        and DN is the Digital Numbers in the tiff files.

        Also the noise should be subtracted.

        See
        https://sentinel.esa.int/web/sentinel/sentinel-1-sar-wiki/-/wiki/Sentinel%20One/Application+of+Radiometric+Calibration+LUT

        The noise correction/subtraction is implemented in an independent package "sentinel1denoised"
        See
        https://github.com/nansencenter/sentinel1denoised
        '''

        # Get look direction
        sat_heading = initial_bearing(longitude[:-1, :],
                                      latitude[:-1, :],
                                      longitude[1:, :],
                                      latitude[1:, :])
        look_direction = scipy.ndimage.interpolation.zoom(
            np.mod(sat_heading + 90, 360),
            (np.shape(longitude)[0] / (np.shape(longitude)[0]-1.), 1))

        # Decompose, to avoid interpolation errors around 0 <-> 360
        look_direction_u = np.sin(np.deg2rad(look_direction))
        look_direction_v = np.cos(np.deg2rad(look_direction))
        look_u_VRT = VRT.from_array(look_direction_u)
        look_v_VRT = VRT.from_array(look_direction_v)
        lookVRT = VRT.from_lonlat(longitude, latitude)
        lookVRT.create_band([{'SourceFilename': look_u_VRT.filename,
                               'SourceBand': 1},
                              {'SourceFilename': look_v_VRT.filename,
                               'SourceBand': 1}],
                             {'PixelFunctionType': 'UVToDirectionTo'}
                             )

        # Blow up to full size
        lookVRT = lookVRT.get_resized_vrt(self.dataset.RasterXSize, self.dataset.RasterYSize, 1)

        # Store VRTs so that they are accessible later
        self.band_vrts['look_u_VRT'] = look_u_VRT
        self.band_vrts['look_v_VRT'] = look_v_VRT
        self.band_vrts['lookVRT'] = lookVRT

        metaDict = []
        # Add bands to full size VRT
        for key in pol:
            name = 'LUT_sigmaNought_%s' % pol[key]
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append(
                {'src': {'SourceFilename':
                         (self.band_vrts['LUT_sigmaNought_VRT_' +
                          pol[key]].filename),
                         'SourceBand': 1
                         },
                 'dst': {'name': name
                         }
                 })
            name = 'LUT_noise_%s' % pol[key]
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append({
                'src': {
                    'SourceFilename': self.band_vrts['LUT_noise_VRT_' +
                                                   pol[key]].filename,
                    'SourceBand': 1
                },
                'dst': {
                    'name': name
                }
            })

        name = 'look_direction'
        bandNumberDict[name] = bnmax+1
        bnmax = bandNumberDict[name]
        metaDict.append({
            'src': {
                'SourceFilename': self.band_vrts['lookVRT'].filename,
                'SourceBand': 1
            },
            'dst': {
                'wkv': 'sensor_azimuth_angle',
                'name': name
            }
        })

        for key in gdalDatasets.keys():
            dsPath, dsName = os.path.split(mdsDict[key])
            name = 'sigma0_%s' % pol[key]
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append(
                {'src': [{'SourceFilename': self.filename,
                          'SourceBand': bandNumberDict['DN_%s' % pol[key]],
                          },
                         {'SourceFilename':
                          (self.band_vrts['LUT_sigmaNought_VRT_%s'
                           % pol[key]].filename),
                          'SourceBand': 1
                          }
                         ],
                 'dst': {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                         'PixelFunctionType': 'Sentinel1Calibration',
                         'polarization': pol[key],
                         'suffix': pol[key],
                         },
                 })
            name = 'beta0_%s' % pol[key]
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append(
                {'src': [{'SourceFilename': self.filename,
                          'SourceBand': bandNumberDict['DN_%s' % pol[key]]
                          },
                         {'SourceFilename':
                          (self.band_vrts['LUT_betaNought_VRT_%s'
                           % pol[key]].filename),
                          'SourceBand': 1
                          }
                         ],
                 'dst': {'wkv': 'surface_backwards_brightness_coefficient_of_radar_wave',
                         'PixelFunctionType': 'Sentinel1Calibration',
                         'polarization': pol[key],
                         'suffix': pol[key],
                         },
                 })

        self.create_bands(metaDict)

        # Add incidence angle as band
        name = 'incidence_angle'
        bandNumberDict[name] = bnmax+1
        bnmax = bandNumberDict[name]
        src = {'SourceFilename': self.band_vrts['incVRT'].filename,
               'SourceBand': 1}
        dst = {'wkv': 'angle_of_incidence',
               'name': name}
        self.create_band(src, dst)
        self.dataset.FlushCache()

        # Add elevation angle as band
        name = 'elevation_angle'
        bandNumberDict[name] = bnmax+1
        bnmax = bandNumberDict[name]
        src = {'SourceFilename': self.band_vrts['eleVRT'].filename,
               'SourceBand': 1}
        dst = {'wkv': 'angle_of_elevation',
               'name': name}
        self.create_band(src, dst)
        self.dataset.FlushCache()

        # Add sigma0_VV
        pp = [pol[key] for key in pol]
        if 'VV' not in pp and 'HH' in pp:
            name = 'sigma0_VV'
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            src = [{'SourceFilename': self.filename,
                    'SourceBand': bandNumberDict['DN_HH'],
                    },
                   {'SourceFilename': (self.band_vrts['LUT_sigmaNought_VRT_HH'].
                                       filename),
                    'SourceBand': 1,
                    },
                   {'SourceFilename': self.band_vrts['incVRT'].filename,
                    'SourceBand': 1}
                   ]
            dst = {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                   'PixelFunctionType': 'Sentinel1Sigma0HHToSigma0VV',
                   'polarization': 'VV',
                   'suffix': 'VV'}
            self.create_band(src, dst)
            self.dataset.FlushCache()

        # set time as acquisition start time
        n = Node.create(self.manifestXML)
        meta = n.node('metadataSection')
        for nn in meta.children:
            if nn.getAttribute('ID') == u'acquisitionPeriod':
                # set valid time
                self.dataset.SetMetadataItem(
                    'time_coverage_start',
                    parse((nn.node('metadataWrap').
                           node('xmlData').
                           node('safe:acquisitionPeriod')['safe:startTime'])
                          ).isoformat())
                self.dataset.SetMetadataItem(
                    'time_coverage_end',
                    parse((nn.node('metadataWrap').
                           node('xmlData').
                           node('safe:acquisitionPeriod')['safe:stopTime'])
                          ).isoformat())

        # Get dictionary describing the instrument and platform according to
        # the GCMD keywords
        mm = pti.get_gcmd_instrument('sar')
        ee = pti.get_gcmd_platform(missionName)

        # TODO: Validate that the found instrument and platform are indeed what we
        # want....

        self.dataset.SetMetadataItem('instrument', json.dumps(mm))
        self.dataset.SetMetadataItem('platform', json.dumps(ee))

    def get_LUT_VRTs(self, XML, vectorListName, LUT_list):
        n = Node.create(XML)
        vecList = n.node(vectorListName)
        X = []
        Y = []
        LUTs = {}
        for LUT in LUT_list:
            LUTs[LUT] = []
        xLengths = []
        for vec in vecList.children:
            xVec = list(map(int, vec['pixel'].split()))
            xLengths.append(len(xVec))
            X.append(xVec)
            Y.append(int(vec['line']))
            for LUT in LUT_list:
                LUTs[LUT].append(list(map(float, vec[LUT].split())))

        # truncate X and LUT to minimum length for all rows
        minLength = np.min(xLengths)
        X = [x[:minLength] for x in X]
        for LUT in LUT_list:
            LUTs[LUT] = [lut[:minLength] for lut in LUTs[LUT]]

        X = np.array(X)
        for LUT in LUT_list:
            LUTs[LUT] = np.array(LUTs[LUT])
        Ym = np.array([Y, ]*np.shape(X)[1]).transpose()

        lon, lat = self.transform_points(X.flatten(), Ym.flatten())
        longitude = lon.reshape(X.shape)
        latitude = lat.reshape(X.shape)

        LUT_VRTs = {}
        for LUT in LUT_list:
            LUT_VRTs[LUT] = VRT.from_array(LUTs[LUT])

        return LUT_VRTs, longitude, latitude

    def read_geolocation_lut(self, annotXML):
        ''' Read lon, lat, pixel, line, ia, ea from XML string <annotXML>'''
        xml = Node.create(annotXML)
        geolocationGridPointList = xml.node('geolocationGrid').node('geolocationGridPointList').children
        X = []
        Y = []
        lon = []
        lat = []
        inc = []
        ele = []
        for gridPoint in geolocationGridPointList:
            X.append(gridPoint['pixel'])
            Y.append(gridPoint['line'])
            lon.append(gridPoint['longitude'])
            lat.append(gridPoint['latitude'])
            inc.append(gridPoint['incidenceAngle'])
            ele.append(gridPoint['elevationAngle'])

        X = np.fromiter(X, float)
        Y = np.fromiter(Y, float)
        lon = np.fromiter(lon, float)
        lat = np.fromiter(lat, float)
        inc = np.fromiter(inc, float)
        ele = np.fromiter(ele, float)

        numberOfSamples = int(xml.node('imageAnnotation').node('imageInformation').node('numberOfSamples').value)
        numberOfLines = int(xml.node('imageAnnotation').node('imageInformation').node('numberOfLines').value)

        return X, Y, lon, lat, inc, ele, numberOfSamples, numberOfLines

    def init_from_manifest_only(self, manifestXML, annotXML, missionName):
        ''' Create fake VRT and add metadata only from the manifest.safe '''
        X, Y, lon, lat, inc, ele, numberOfSamples, numberOfLines = self.read_geolocation_lut(annotXML)

        super(Mapper, self).__init__(numberOfSamples, numberOfLines)
        doc = ET.fromstring(manifestXML)

        gcps = []
        for i in range(len(X)):
            gcps.append(gdal.GCP(lon[i], lat[i], 0, X[i], Y[i]))

        self.dataset.SetGCPs(gcps, NSR().wkt)
        self.dataset.SetMetadataItem('time_coverage_start',
                                     doc.findall(".//*[{http://www.esa.int/safe/sentinel-1.0}startTime]")[0][0].text)
        self.dataset.SetMetadataItem('time_coverage_end',
                                     doc.findall(".//*[{http://www.esa.int/safe/sentinel-1.0}stopTime]")[0][0].text)
        self.dataset.SetMetadataItem('platform', json.dumps(pti.get_gcmd_platform(missionName)))
        self.dataset.SetMetadataItem('instrument', json.dumps(pti.get_gcmd_instrument('SAR')))
        self.dataset.SetMetadataItem('Entry Title', missionName + ' SAR')
        self.dataset.SetMetadataItem('Data Center', 'ESA/EO')
        self.dataset.SetMetadataItem('ISO Topic Category', 'Oceans')
        self.dataset.SetMetadataItem('Summary', missionName + ' SAR data')


