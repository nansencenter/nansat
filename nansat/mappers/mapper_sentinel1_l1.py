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
from nansat.utils import gdal, initial_bearing
from nansat.exceptions import WrongMapperError, NansatReadError
from nansat.nsr import NSR
from nansat.node import Node


class Mapper(VRT):
    """ Create VRT with mapping of Sentinel-1 (A and B) stripmap mode (S1A_SM)

    Parameters
    ----------
    filename : str
        name of input Sentinel-1 L1 file
    gdalDataset : None
    gdalMetadata : None
    fast : bool
        Flag that triggers faster reading of metadata from Sentinel-1 file.
        If True, no bands are added to the dataset and georeference is not corrected.
        If False, all bands are added and GCPs are corrected if necessary
        (see Mapper.correct_geolocation_data for details).

    Note
    ----
    Creates self.dataset and populates it with S1 bands (when fast=False).
    """
    def __init__(self, filename, gdalDataset, gdalMetadata, fast=False, fixgcp=True, **kwargs):
        if not os.path.split(filename.rstrip('/'))[1][:3] in ['S1A', 'S1B']:
            raise WrongMapperError('%s: Not Sentinel 1A or 1B' %filename)

        if not IMPORT_SCIPY:
            raise NansatReadError('Sentinel-1 data cannot be read because scipy is not installed')

        if zipfile.is_zipfile(filename):
            zz = zipfile.PyZipFile(filename)
            # Assuming the file names are consistent, the polarization
            # dependent data should be sorted equally such that we can use the
            # same indices consistently for all the following lists
            # THIS IS NOT THE CASE...
            mds_files = ['/vsizip/%s/%s' % (filename, fn)
                        for fn in zz.namelist() if 'measurement/s1' in fn]
            calibration_files = ['/vsizip/%s/%s' % (filename, fn)
                        for fn in zz.namelist()
                        if 'annotation/calibration/calibration-s1' in fn]
            noise_files = ['/vsizip/%s/%s' % (filename, fn)
                          for fn in zz.namelist()
                          if 'annotation/calibration/noise-s1' in fn]
            annotation_files = ['/vsizip/%s/%s' % (filename, fn)
                               for fn in zz.namelist()
                               if 'annotation/s1' in fn]
            manifest_files = ['/vsizip/%s/%s' % (filename, fn)
                            for fn in zz.namelist()
                            if 'manifest.safe' in fn]
            zz.close()
        else:
            mds_files = glob.glob('%s/measurement/s1*' % filename)
            calibration_files = glob.glob('%s/annotation/calibration/calibration-s1*'
                                 % filename)
            noise_files = glob.glob('%s/annotation/calibration/noise-s1*'
                                   % filename)
            annotation_files = glob.glob('%s/annotation/s1*'
                                        % filename)
            manifest_files = glob.glob('%s/manifest.safe' % filename)

        if (not mds_files or not calibration_files or not noise_files or
            not annotation_files or not manifest_files):
            raise WrongMapperError(filename)

        # convert list of MDS files into dictionary. Keys - polarizations in upper case.
        mds_files = {os.path.basename(ff).split('-')[3].upper():ff for ff in mds_files}
        polarizations = list(mds_files.keys())

        # read annotation files
        self.annotation_data = self.read_annotation(annotation_files)
        if not fast and fixgcp:
            self.correct_geolocation_data()

        # read manifest file
        manifest_data = self.read_manifest_data(manifest_files[0])

        # very fast constructor without any bands only with some metadata and geolocation
        self._init_empty(manifest_data, self.annotation_data, filename)

        # skip adding bands in the fast mode and RETURN
        if fast:
            return

        # Open data files with GDAL
        gdalDatasets = {}
        for pol in polarizations:
            gdalDatasets[pol] = gdal.Open(mds_files[pol])

            if not gdalDatasets[pol]:
                raise WrongMapperError('%s: No Sentinel-1 datasets found' % mds_files[pol])

        # Check metadata to confirm it is Sentinel-1 L1
        metadata = gdalDatasets[polarizations[0]].GetMetadata()

        # create full size VRTs with incidenceAngle and elevationAngle
        annotation_vrts = self.vrts_from_arrays(self.annotation_data,
                                                ['incidenceAngle', 'elevationAngle'])
        self.band_vrts.update(annotation_vrts)

        # create full size VRTS with calibration LUT
        calibration_names = ['sigmaNought', 'betaNought']
        calibration_list_tag = 'calibrationVectorList'
        for calibration_file in calibration_files:
            pol = '_' + os.path.basename(calibration_file).split('-')[4].upper()
            xml = self.read_vsi(calibration_file)
            calibration_data = self.read_calibration(xml, calibration_list_tag, calibration_names, pol)
            calibration_vrts = self.vrts_from_arrays(calibration_data, calibration_names, pol, True, 1)
            self.band_vrts.update(calibration_vrts)

        # create full size VRTS with noise LUT
        for noise_file in noise_files:
            pol = '_' + os.path.basename(noise_file).split('-')[4].upper()
            xml = self.read_vsi(noise_file)
            if '<noiseVectorList' in xml:
                noise_list_tag = 'noiseVectorList'
                noise_name = 'noiseLut'
            elif '<noiseRangeVectorList' in xml:
                noise_list_tag = 'noiseRangeVectorList'
                noise_name = 'noiseRangeLut'
            noise_data = self.read_calibration(xml, noise_list_tag, [noise_name], pol)
            noise_vrts = self.vrts_from_arrays(noise_data, [noise_name], pol, True, 1)
            self.band_vrts.update(noise_vrts)

        #### Create metaDict: dict with metadata for all bands
        metaDict = []
        bandNumberDict = {}
        bnmax = 0
        for pol in polarizations:
            dsPath, dsName = os.path.split(mds_files[pol])
            name = 'DN_%s' % pol
            # A dictionary of band numbers is needed for the pixel function
            # bands further down. This is not the best solution. It would be
            # better to have a function in VRT that returns the number given a
            # band name. This function exists in Nansat but could perhaps be
            # moved to VRT? The existing nansat function could just call the
            # VRT one...
            bandNumberDict[name] = bnmax + 1
            bnmax = bandNumberDict[name]
            band = gdalDatasets[pol].GetRasterBand(1)
            dtype = band.DataType
            metaDict.append({
                'src': {
                    'SourceFilename': mds_files[pol],
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
        longitude, latitude = self.transform_points(calibration_data['pixel'].flatten(),
                                                    calibration_data['line'].flatten())
        longitude.shape = calibration_data['pixel'].shape
        latitude.shape = calibration_data['pixel'].shape
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
        for pol in polarizations:
            name = 'sigmaNought_%s' % pol
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append(
                {'src': {'SourceFilename':
                         (self.band_vrts[name].filename),
                         'SourceBand': 1
                         },
                 'dst': {'name': name
                         }
                 })
            name = 'noise_%s' % pol
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append({
                'src': {
                    'SourceFilename': self.band_vrts['%s_%s' % (noise_name, pol)].filename,
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

        for pol in polarizations:
            dsPath, dsName = os.path.split(mds_files[pol])
            name = 'sigma0_%s' % pol
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append(
                {'src': [{'SourceFilename': self.filename,
                          'SourceBand': bandNumberDict['DN_%s' % pol],
                          },
                         {'SourceFilename': self.band_vrts['sigmaNought_%s' % pol].filename,
                          'SourceBand': 1
                          }
                         ],
                 'dst': {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                         'PixelFunctionType': 'Sentinel1Calibration',
                         'polarization': pol,
                         'suffix': pol,
                         },
                 })
            name = 'beta0_%s' % pol
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append(
                {'src': [{'SourceFilename': self.filename,
                          'SourceBand': bandNumberDict['DN_%s' % pol]
                          },
                         {'SourceFilename': self.band_vrts['betaNought_%s' % pol].filename,
                          'SourceBand': 1
                          }
                         ],
                 'dst': {'wkv': 'surface_backwards_brightness_coefficient_of_radar_wave',
                         'PixelFunctionType': 'Sentinel1Calibration',
                         'polarization': pol,
                         'suffix': pol,
                         },
                 })

        self.create_bands(metaDict)

        # Add incidence angle as band
        name = 'incidence_angle'
        bandNumberDict[name] = bnmax+1
        bnmax = bandNumberDict[name]
        src = {'SourceFilename': self.band_vrts['incidenceAngle'].filename,
               'SourceBand': 1}
        dst = {'wkv': 'angle_of_incidence',
               'name': name}
        self.create_band(src, dst)
        self.dataset.FlushCache()

        # Add elevation angle as band
        name = 'elevation_angle'
        bandNumberDict[name] = bnmax+1
        bnmax = bandNumberDict[name]
        src = {'SourceFilename': self.band_vrts['elevationAngle'].filename,
               'SourceBand': 1}
        dst = {'wkv': 'angle_of_elevation',
               'name': name}
        self.create_band(src, dst)
        self.dataset.FlushCache()

        # Add sigma0_VV
        if 'VV' not in polarizations and 'HH' in polarizations:
            name = 'sigma0_VV'
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            src = [{'SourceFilename': self.filename,
                    'SourceBand': bandNumberDict['DN_HH'],
                    },
                   {'SourceFilename': (self.band_vrts['sigmaNought_HH'].
                                       filename),
                    'SourceBand': 1,
                    },
                   {'SourceFilename': self.band_vrts['incidenceAngle'].filename,
                    'SourceBand': 1}
                   ]
            dst = {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                   'PixelFunctionType': 'Sentinel1Sigma0HHToSigma0VV',
                   'polarization': 'VV',
                   'suffix': 'VV'}
            self.create_band(src, dst)
            self.dataset.FlushCache()


    def read_calibration(self, xml, vectorListName, variable_names, pol):
        """ Read calibration data from calibration or noise XML files
        Parameters
        ----------
        xml : str
            String with XML from calibration or noise files
        vectorListName : str
            tag of the element that contains lists with LUT values
        variable_names : list of str
            names of LUT variable to read
        pol : str
            HH, HV, etc

        Returns
        -------
        data : dict
            Calibration or noise data. Keys:
            The same as variable_names + 'pixel', 'line'
        """
        data = {}
        n = Node.create(xml)
        vecList = n.node(vectorListName)
        data['pixel'] = []
        data['line'] = []
        for var_name in variable_names:
            data[var_name+pol] = []
        xLengths = []
        for vec in vecList.children:
            xVec = list(map(int, vec['pixel'].split()))
            xLengths.append(len(xVec))
            data['pixel'].append(xVec)
            data['line'].append(int(vec['line']))
            for var_name in variable_names:
                data[var_name+pol].append(np.fromiter(vec[var_name].split(), float))

        # truncate data['pixel'] and var_name to minimum length for all rows
        minLength = np.min(xLengths)
        data['pixel'] = [x[:minLength] for x in data['pixel']]
        for var_name in variable_names:
            data[var_name+pol] = [d[:minLength] for d in data[var_name+pol]]

        data['pixel'] = np.array(data['pixel'])
        for var_name in variable_names:
            data[var_name+pol] = np.array(data[var_name+pol])
        data['line'] = np.array([data['line'], ]*np.shape(data['pixel'])[1]).transpose()

        return data

    def read_annotation(self, annotation_files):
        """ Read lon, lat, etc from annotation XML

        Parameters
        ----------
        annotation_files : list
            strings with names of annotation files

        Returns
        -------
        data : dict
            geolocation data from the XML as 2D np.arrays. Keys:
                pixel, line, longitude, latitude, height, incidenceAngle, elevationAngle: 2D arrays
                shape : tuple (shape of geolocation data arrays)
                x_size, y_size : int
                pol : list

        """
        variable_names = ['pixel', 'line', 'longitude', 'latitude', 'height', 'incidenceAngle', 'elevationAngle']
        data = {var_name:[] for var_name in variable_names}

        xml = self.read_vsi(annotation_files[0])
        xml = Node.create(xml)
        geolocation_points = xml.node('geolocationGrid').node('geolocationGridPointList').children
        # collect data from XML into dictionary with lists
        for point in geolocation_points:
            for var_name in variable_names:
                data[var_name].append(point[var_name])

        # convert lists to 1D arrays
        for var_name in variable_names:
            data[var_name] = np.fromiter(data[var_name], float)

        # get shape of geolocation matrix (number of occurence of minimal element)
        data['shape'] = (data['pixel']==0).sum(), (data['line']==0).sum()

        # convert 1D arrays to 2D
        for var_name in variable_names:
            data[var_name].shape = data['shape']

        # get raster dimentions
        image_info = xml.node('imageAnnotation').node('imageInformation')
        data['x_size'] = int(image_info.node('numberOfSamples').value)
        data['y_size'] = int(image_info.node('numberOfLines').value)

        # get list of polarizations
        data['pol'] = []
        for ff in annotation_files:
            p = os.path.basename(ff).split('-')[3]
            data['pol'].append(p.upper())

        return data

    def _init_empty(self, manifest_data, annotation_data, filename):
        """ Fast initialization from minimum of information

        Parameters
        ----------
        manifest_data : dict
            data from the manifest file (time_coverage_start, etc)
        annotation_data : dict
            data from annotation file (longitude, latitude, x_size, etc)
        'entry_id' is the unique id of each dataset which is set identical
            to the filename for this mapper : str

        Note
        ----
            Calls VRT.__init__, Adds GCPs, metadata
        """
        # init empty dataset
        super(Mapper, self).__init__(annotation_data['x_size'], annotation_data['y_size'])
        # add GCPs from (corrected) geolocation data
        gcps = Mapper.create_gcps(annotation_data['longitude'],
                                  annotation_data['latitude'],
                                  annotation_data['height'],
                                  annotation_data['pixel'],
                                  annotation_data['line'])
        self.dataset.SetGCPs(gcps, NSR().wkt)
        # set metadata
        self.dataset.SetMetadataItem('time_coverage_start', manifest_data['time_coverage_start'])
        self.dataset.SetMetadataItem('time_coverage_end', manifest_data['time_coverage_end'])
        platform_name = manifest_data['platform_family_name'] + manifest_data['platform_number']
        self.dataset.SetMetadataItem('platform', json.dumps(pti.get_gcmd_platform(platform_name)))
        self.dataset.SetMetadataItem('instrument', json.dumps(pti.get_gcmd_instrument('SAR')))
        self.dataset.SetMetadataItem('entry_title', platform_name + ' SAR')
        self.dataset.SetMetadataItem('data_center', json.dumps(pti.get_gcmd_provider('ESA/EO')))
        self.dataset.SetMetadataItem('iso_topic_category', json.dumps(pti.get_iso19115_topic_category('Oceans')))
        self.dataset.SetMetadataItem('summary', platform_name + ' SAR data')
        self.dataset.SetMetadataItem('entry_id',  os.path.splitext(os.path.basename(filename))[0].upper())
        self.dataset.FlushCache()

    def correct_geolocation_data(self):
        """ Correct lon/lat values in geolocation data for points high above ground (incorrect)

        Each GCP in Sentinel-1 L1 image (both in the GeoTIF files and Annotation LUT) have five
        coordinates: X, Y, Z (height), Pixel and Line. On some scenes that cover Greenland (and
        probably other lands) some GCPs have height above zero even over ocean. This is incorrect,
        because the radar signal comes actually from the surface and not from a point above the
        ground as stipulated in such GCPs. Correction of GCPs in this function is equivalelnt to
        reverse DEM correction of SAR data.

        Notes
        -------
        Updates 'pixel' and 'height' in self.annotation_data

        """
        PIXEL_SIZE = 40. # meters
        pix_corr = (self.annotation_data['height'] /
                    np.tan(np.radians(self.annotation_data['incidenceAngle'])) / PIXEL_SIZE)
        self.annotation_data['pixel'] += pix_corr
        self.annotation_data['height'][:] = 0

    @staticmethod
    def create_gcps(x, y, z, p, l):
        """ Create GCPs from geolocation data

        Parameters
        ----------
        x, y, z, p, l
        N-D arrays with value of X, Y, Z, Pixel and Line coordinates.
        X and Y are typically lon, lat, Z - height.

        Returns
        -------
        gcps : list with GDAL GCPs

        """
        gcps = []
        for xi, yi, zi, pi, li in zip(x.flat, y.flat, z.flat, p.flat, l.flat):
            gcps.append(gdal.GCP(xi, yi, zi, pi, li))
        return gcps

    def read_manifest_data(self, input_file):
        """ Read information (time_coverage_start, etc) manifest XML

        Parameters
        ----------
        input_file : str
            name of manifest file

        Returns
        -------
        data : dict
            manifest data. Keys:
                time_coverage_start
                time_coverage_end
                platform_familyName
                platform_number

        """

        data = {}
        xml = self.read_vsi(input_file)
        # set time as acquisition start time
        n = Node.create(xml)
        meta = n.node('metadataSection')
        for nn in meta.children:
            if str(nn.getAttribute('ID')) == 'acquisitionPeriod':
                # get valid time
                data['time_coverage_start'] = parse((nn.node('metadataWrap').
                                                     node('xmlData').
                                                     node('safe:acquisitionPeriod')['safe:startTime']
                                                     )).isoformat()
                data['time_coverage_end'] = parse((nn.node('metadataWrap').
                                                   node('xmlData').
                                                   node('safe:acquisitionPeriod')['safe:stopTime']
                                                   )).isoformat()
            if str(nn.getAttribute('ID')) == 'platform':
                data['platform_family_name'] = str(nn.node('metadataWrap').
                                                     node('xmlData').
                                                     node('safe:platform')['safe:familyName'])
                data['platform_number'] = str(nn.node('metadataWrap').
                                                     node('xmlData').
                                                     node('safe:platform')['safe:number'])
        return data

    def vrts_from_arrays(self, data, variable_names, pol='', resize=True, resample_alg=2):
        """ Convert input dict with arrays into dict with VRTs

        Parameters
        ----------
        data : dict
            2D arrays with data from LUT
        variable_names : list of str
            variable names that should be converted to VRTs
        pol : str
            HH, HV, etc
        resize : bool
            Shall VRT be zoomed to full size?
        resample_alg : int
            Index of resampling algorithm. See VRT.get_resized_vrt()

        Returns
        -------
        vrts : dict with (resized) VRTs

        """
        vrts = {}
        for var_name in variable_names:
            vrts[var_name+pol] = VRT.from_array(data[var_name+pol])
            if resize:
                vrts[var_name+pol] = vrts[var_name+pol].get_resized_vrt(self.dataset.RasterXSize,
                                                                        self.dataset.RasterYSize,
                                                                        resample_alg)
        return vrts
