#-------------------------------------------------------------------------------
# Name:		mapper_s1a_l1.py
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	12.09.2014
# Last modified:19.09.2014 17:08
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------

import os, zipfile
import numpy as np
import scipy
from dateutil.parser import parse

from nansat.vrt import VRT
from nansat.tools import gdal, WrongMapperError, initial_bearing
from nansat.nsr import NSR
from nansat.node import Node

class Mapper(VRT):
    '''
        Create VRT with mapping of Sentinel-1A stripmap mode (S1A_SM)
    '''

    def __init__(self, fileName, gdalDataset, gdalMetadata, product_type='RVL',
            GCP_COUNT=10, **kwargs):
        '''
        Parameters
        ----------
        product_type: string
            Sentinel-1 level-1 product, i.e. 
        GCP_COUNT : int
            number of GCPs along each dimention
        '''

        gdalDatasets = []
        mds_filenames = []
        polarizations = []
        if zipfile.is_zipfile(fileName):
            zz = zipfile.PyZipFile(fileName)
            mds = [fn for fn in zz.namelist() if 'measurement/s1a' in fn]
            if not mds:
                raise WrongMapperError
            for i,mm in enumerate(mds):
                # Open zip file using VSI
                mds_filenames.append('/vsizip/%s/%s' % (fileName, mm))
                gdalDatasets.append(gdal.Open(mds_filenames[i]))

                # read noise and calibration xml-files in folder annotation/calibration
                dsPath, dsName = os.path.split(mm)
                pol = dsName[11:13]
                polarizations.append(pol.upper())
                calFile = [fn for fn in zz.namelist() if \
                        'annotation/calibration/calibration-s1a' in fn and \
                        pol in fn]
                noiseFile = [fn for fn in zz.namelist() if \
                        'annotation/calibration/noise-s1a' in fn and \
                        pol in fn]
                annotationFile = [fn for fn in zz.namelist() if \
                        'annotation/s1a' in fn and \
                        pol in fn]
                calXML = self.read_xml('/vsizip/%s/%s' %(fileName,calFile[0]))
                noiseXML = self.read_xml('/vsizip/%s/%s'
                        %(fileName,noiseFile[0]))
                annotationXML = self.read_xml('/vsizip/%s/%s'
                        %(fileName,annotationFile[0]))
            zz.close()

        if not gdalDatasets:
            raise WrongMapperError('No Sentinel-1 datasets found')

        # Check metadata to confirm it is Sentinel-1 L1
        metadata = gdalDatasets[0].GetMetadata()
        if not 'TIFFTAG_IMAGEDESCRIPTION' in metadata.keys():
            raise WrongMapperError
        if not 'Sentinel-1' in metadata['TIFFTAG_IMAGEDESCRIPTION'] \
                and not 'L1' in metadata['TIFFTAG_IMAGEDESCRIPTION']:
            raise WrongMapperError

        ''' Other metadata:

            - incidence angles, pass direction + more in annotation folder
            - the pass direction can be found in manifest.safe as s1:pass

        '''

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalDatasets[0])

        # set time
        self._set_time(parse(metadata['TIFFTAG_DATETIME']))

        metaDict = []
        bandNumberDict = {}
        for i,dataset in enumerate(gdalDatasets):
            dsPath, dsName = os.path.split(mds[i])
            suffix = dsName[11:13].upper() # polarization
            name = 'DN_%s' %suffix
            # A dictionary of band numbers is needed for the pixel function
            # bands further down. This is not the best solution. It would be
            # better to have a function in VRT that returns the number given a
            # band name. This function exists in Nansat but could perhaps be
            # moved to VRT? The existing nansat function could just call the
            # VRT one...
            bandNumberDict[name] = i + 1
            bnmax = bandNumberDict[name]
            band = dataset.GetRasterBand(1)
            dtype = band.DataType
            metaDict.append({
                'src': {
                    'SourceFilename': mds_filenames[i],
                    'SourceBand': 1,
                    'DataType': dtype,
                },
                'dst': {
                    'name': name,
                    'SourceTransferType': gdal.GetDataTypeName(dtype),
                    'dataType': 6,
                },
            })
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        '''
        Calibration should be performed as

        s0 = DN^2/sigmaNought^2,
        
        where sigmaNought is from e.g.
        annotation/calibration/calibration-s1a-iw-grd-hh-20140811t151231-20140811t151301-001894-001cc7-001.xml,
        and DN is the Digital Numbers in the tiff files.

        Also the noise should be subtracted.

        See
        https://sentinel.esa.int/web/sentinel/sentinel-1-sar-wiki/-/wiki/Sentinel%20One/Application+of+Radiometric+Calibration+LUT
        '''
        # Get calibration data
        calibration_LUT_VRTs, longitude, latitude  = self.get_LUT_VRTs(calXML,
                'calibrationVectorList', 
                ['sigmaNought', 'betaNought', 'gamma', 'dn'])
        # Get noise data
        noise_LUT_VRTs = self.get_LUT_VRTs(noiseXML,
                'noiseVectorList',
                ['noiseLut'])[0]

        # Get look direction
        sat_heading = initial_bearing(longitude[:-1,:], latitude[:-1,:],
                longitude[1:,:], latitude[1:,:])
        look_direction = scipy.ndimage.interpolation.zoom( np.mod(sat_heading +
            90, 360), (np.shape(longitude)[0]/(np.shape(longitude)[0]-1.), 1) )
        # Decompose, to avoid interpolation errors around 0 <-> 360
        look_direction_u = np.sin(np.deg2rad(look_direction))
        look_direction_v = np.cos(np.deg2rad(look_direction))
        look_u_VRT = VRT(array=look_direction_u, lat=latitude,
                lon=longitude)
        look_v_VRT = VRT(array=look_direction_v, lat=latitude,
                lon=longitude)
        lookVRT = VRT(lat=latitude, lon=longitude)
        lookVRT._create_band(
                [{'SourceFilename': look_u_VRT.fileName, 'SourceBand': 1},
                 {'SourceFilename': look_v_VRT.fileName, 'SourceBand': 1}],
                 {'PixelFunctionType': 'UVToDirectionTo'})

        # Blow up to full size
        lookVRT = lookVRT.get_resized_vrt(self.dataset.RasterXSize,
                self.dataset.RasterYSize, eResampleAlg=1)
        LUT_sigmaNought_VRT = \
                calibration_LUT_VRTs['sigmaNought'].get_resized_vrt(
                self.dataset.RasterXSize,
                self.dataset.RasterYSize, eResampleAlg=1)
        LUT_betaNought_VRT = \
                calibration_LUT_VRTs['betaNought'].get_resized_vrt(
                self.dataset.RasterXSize,
                self.dataset.RasterYSize, eResampleAlg=1)
        LUT_noise_VRT = noise_LUT_VRTs['noiseLut'].get_resized_vrt(
                self.dataset.RasterXSize,
                self.dataset.RasterYSize, eResampleAlg=1)

        # Store VRTs so that they are accessible later
        self.subVRTs = {
                'look_u_VRT': look_u_VRT,
                'look_v_VRT': look_v_VRT,
                'lookVRT': lookVRT,
                'LUT_sigmaNought_VRT': LUT_sigmaNought_VRT,
                'LUT_betaNought_VRT': LUT_betaNought_VRT,
                'LUT_gamma_VRT': calibration_LUT_VRTs['gamma'],
                'LUT_dn_VRT': calibration_LUT_VRTs['dn'],
                'LUT_noise_VRT': LUT_noise_VRT,
            }

        metaDict = []
        # Add bands to full sized VRT
        name = 'LUT_sigmaNought'
        bandNumberDict[name] = bnmax+1
        bnmax = bandNumberDict[name]
        metaDict.append({
            'src': {
                'SourceFilename': self.subVRTs['LUT_sigmaNought_VRT'].fileName,
                'SourceBand': 1
            },
            'dst': {
                'name': name
            }
        })
        name = 'LUT_noise'
        bandNumberDict[name] = bnmax+1
        bnmax = bandNumberDict[name]
        metaDict.append({
            'src': {
                'SourceFilename': self.subVRTs['LUT_noise_VRT'].fileName,
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
                'SourceFilename': self.subVRTs['lookVRT'].fileName,
                'SourceBand': 1
            },
            'dst': {
                'wkv': 'sensor_azimuth_angle',
                'name': name
            }
        })
        for i,dataset in enumerate(gdalDatasets):
            dsPath, dsName = os.path.split(mds[i])
            suffix = dsName[11:13].upper() # polarization
            name = 'sigma0_'+suffix
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append({
                'src': [{
                    'SourceFilename': self.fileName,
                    'SourceBand': bandNumberDict['DN_%s' %suffix],
                },{
                    'SourceFilename': self.subVRTs['LUT_noise_VRT'].fileName,
                    'SourceBand': 1
                },{
                    'SourceFilename': self.subVRTs['LUT_sigmaNought_VRT'].fileName,
                    'SourceBand': 1
                }],
                'dst': {
                    'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                    'PixelFunctionType': 'Sentinel1Calibration',
                    'polarization': suffix,
                    'suffix': suffix,
                },
            })
            name = 'beta0_'+suffix
            bandNumberDict[name] = bnmax+1
            bnmax = bandNumberDict[name]
            metaDict.append({
                'src': [{
                    'SourceFilename': self.fileName,
                    'SourceBand': bandNumberDict['DN_%s' %suffix],
                },{
                    'SourceFilename': self.subVRTs['LUT_noise_VRT'].fileName,
                    'SourceBand': 1
                },{
                    'SourceFilename': self.subVRTs['LUT_betaNought_VRT'].fileName,
                    'SourceBand': 1
                }],
                'dst': {
                    'wkv': 'radar_brightness_coefficient',
                    'PixelFunctionType': 'Sentinel1Calibration',
                    'polarization': suffix,
                    'suffix': suffix,
                },
            })

        self._create_bands(metaDict)

        # Add incidence angle
        name = 'incidence_angle'
        bandNumberDict[name] = bnmax+1
        bnmax = bandNumberDict[name]
        src = [{
                'SourceFilename': self.fileName,
                'SourceBand': bandNumberDict['DN_%s' %polarizations[0]],
            },{
                'SourceFilename': self.subVRTs['LUT_noise_VRT'].fileName,
                'SourceBand': 1
            },{
                'SourceFilename': self.subVRTs['LUT_sigmaNought_VRT'].fileName,
                'SourceBand': 1
            },{
                 'SourceFilename': self.subVRTs['LUT_betaNought_VRT'].fileName,
                 'SourceBand': 1
        }]
        dst = {'wkv': 'angle_of_incidence',
               'PixelFunctionType': 'Sentinel1IncidenceAngle',
               '_FillValue': -10000,
               'name': name}

        #import pdb
        #pdb.set_trace()

        self._create_band(src, dst)
        self.dataset.FlushCache()

        ##check landsat mapper!

        ## Add sigma0_VV
        #if 'VV' not in pol and 'HH' in pol:
        #    name = 'sigma0_VV'
        #    bandNumberDict[name] = bnmax+1
        #    bnmax = bandNumberDict[name]
        #    src = [{'SourceFilename': self.fileName,
        #            'SourceBand': bandNumberDict['sigma0_HH'],
        #            'DataType': 6},
        #           {'SourceFilename': self.fileName,
        #            'SourceBand': bandNumberDict['beta0_HH'],
        #            'DataType': 6}]
        #    dst = {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
        #           'PixelFunctionType': 'Sigma0HHBetaToSigma0VV',
        #           'polarization': 'VV',
        #           'suffix': 'VV'}
        #    self._create_band(src, dst)
        #    self.dataset.FlushCache()

        self.n = Node.create(annotationXML)
        n = Node.create(annotationXML)
        ga = n.node('generalAnnotation')
        pi = ga.node('productInformation')
        self.dataset.SetMetadataItem('ORBIT_DIRECTION', str(pi['pass']))
        # Incidence angles are found in
        #<geolocationGrid>
        #    <geolocationGridPointList count="252">
        #          <geolocationGridPoint>
        geolocationGridPointList = n.node('geolocationGrid').children[0]
        for gridPoint in geolocationGridPointList.children:
            # Get the data and do like above for the LUTs


    def get_LUT_VRTs(self, XML, vectorListName, LUT_list):
        n = Node.create(XML)
        vecList = n.node(vectorListName)
        X = []
        Y = []
        LUTs = {}
        for LUT in LUT_list:
            LUTs[LUT] = []
        for vec in vecList.children:
            X.append(map(int, vec['pixel'].split()))
            Y.append(int(vec['line']))
            for LUT in LUT_list:
                LUTs[LUT].append(map(float, vec[LUT].split()))

        X = np.array(X)
        for LUT in LUT_list:
            LUTs[LUT] = np.array(LUTs[LUT])
        Ym = np.array([Y,]*np.shape(X)[1]).transpose()

        lon, lat = self.transform_points(X.flatten(), Ym.flatten())
        longitude = lon.reshape(X.shape)
        latitude = lat.reshape(X.shape)

        LUT_VRTs = {}
        for LUT in LUT_list:
            LUT_VRTs[LUT] = VRT(array=LUTs[LUT], lat=latitude, lon=longitude)

        return LUT_VRTs, longitude, latitude


