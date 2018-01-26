# Name:         mapper_CSKS.py
# Purpose:      Mapper for Cosmo-Skymed SAR data
# Authors:      Morten W. Hansen
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html
from dateutil.parser import parse
from struct import unpack

import numpy as np
import os

from nansat.tools import gdal, ogr, osr, parse_time
from nansat.exceptions import WrongMapperError
from nansat.vrt import VRT


class Mapper(VRT):
    ''' VRT with mapping of WKV for Cosmo-Skymed '''

    def __init__(self, filename, gdalDataset, gdalMetadata, **kwargs):
        ''' Create CSKS VRT '''

        if filename.split('/')[-1][0:4] != "CSKS":
            raise WrongMapperError

        # Get coordinates
        metadata = gdalMetadata['Estimated_Bottom_Left_Geodetic_Coordinates']
        bottom_left_lon = float(metadata.split(' ')[1])
        bottom_left_lat = float(metadata.split(' ')[0])
        metadata = gdalMetadata['Estimated_Bottom_Right_Geodetic_Coordinates']
        bottom_right_lon = float(metadata.split(' ')[1])
        bottom_right_lat = float(metadata.split(' ')[0])
        metadata = gdalMetadata['Estimated_Top_Left_Geodetic_Coordinates']
        top_left_lon = float(metadata.split(' ')[1])
        top_left_lat = float(metadata.split(' ')[0])
        metadata = gdalMetadata['Estimated_Top_Right_Geodetic_Coordinates']
        top_right_lon = float(metadata.split(' ')[1])
        top_right_lat = float(metadata.split(' ')[0])
        metadata = gdalMetadata['Scene_Centre_Geodetic_Coordinates']
        center_lon = float(metadata.split(' ')[1])
        center_lat = float(metadata.split(' ')[0])

        # Get sub-datasets
        subDatasets = gdalDataset.GetSubDatasets()

        # Get file names from dataset or subdataset
        if subDatasets.__len__() == 1:
            filenames = [filename]
        else:
            filenames = [f[0] for f in subDatasets]

        for i, elem in enumerate(filenames):
            if filenames[i][-3:] == 'QLK':
                filenames.pop(i)
        #print filenames

        subDataset = gdal.Open(filenames[0])

        # generate list of GCPs
        gcps = []
        # create GCP with X,Y,Z(?),pixel,line from lat/lon matrices
        gcp = gdal.GCP(float(bottom_left_lon),
                       float(bottom_left_lat), 0, 0, 0)
        gcps.append(gcp)
        #self.logger.debug('%d %d %d %f %f', 0, gcp.GCPPixel, gcp.GCPLine,
        #                  gcp.GCPX, gcp.GCPY)
        gcp = gdal.GCP(float(bottom_right_lon),
                       float(bottom_right_lat),
                       0,
                       subDataset.RasterXSize,
                       0)
        gcps.append(gcp)
        #self.logger.debug('%d %d %d %f %f', 1, gcp.GCPPixel, gcp.GCPLine,
        #                  gcp.GCPX, gcp.GCPY)
        gcp = gdal.GCP(float(top_left_lon),
                       float(top_left_lat),
                       0,
                       0,
                       subDataset.RasterYSize)
        gcps.append(gcp)
        #self.logger.debug('%d %d %d %f %f', 2, gcp.GCPPixel, gcp.GCPLine,
        #                  gcp.GCPX, gcp.GCPY)
        gcp = gdal.GCP(float(top_right_lon),
                       float(top_right_lat),
                       0,
                       subDataset.RasterXSize,
                       subDataset.RasterYSize)
        gcps.append(gcp)
        #self.logger.debug('%d %d %d %f %f', 3, gcp.GCPPixel, gcp.GCPLine,
        #                  gcp.GCPX, gcp.GCPY)
        gcp = gdal.GCP(float(center_lon),
                       float(center_lat),
                       0,
                       int(np.round(subDataset.RasterXSize/2.)),
                       int(round(subDataset.RasterYSize/2.)))
        gcps.append(gcp)
        #self.logger.debug('%d %d %d %f %f', 4, gcp.GCPPixel, gcp.GCPLine,
        #                  gcp.GCPX, gcp.GCPY)

        # append GCPs and lat/lon projection to the vsiDataset
        latlongSRS = osr.SpatialReference()
        latlongSRS.ImportFromProj4("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
        latlongSRSWKT = latlongSRS.ExportToWkt()

        # create empty VRT dataset with geolocation only
        # x_size, y_size, geo_transform, projection, gcps=None, gcp_projection='', **kwargs
        self._init_from_dataset_params(subDataset.RasterXSize, subDataset.RasterYSize,
                                        (0,1,0,subDataset.RasterYSize,0,-1),
                                        latlongSRSWKT, gcps, latlongSRSWKT)

        #print self.filename
        # Read all bands later
        #band='S01'
        #res='SBI'

        # Use only full size "original" datasets
        for i, elem in enumerate(filenames):
            if filenames[i][-3:] == 'SBI':
                # Add real and imaginary raw counts as bands
                src = {'SourceFilename': filenames[i],
                       'SourceBand': 1,
                       'DataType': gdal.GDT_Int16}
                dst = {'dataType': gdal.GDT_Float32,
                       'name': 'RawCounts_%s_real' %
                       gdalMetadata[filenames[i][-7:-4]+'_Polarisation']}
                self.create_band(src, dst)

                src = {'SourceFilename': filenames[i],
                       'SourceBand': 2,
                       'DataType': gdal.GDT_Int16}
                dst = {'dataType': gdal.GDT_Float32,
                       'name': 'RawCounts_%s_imaginary' %
                       gdalMetadata[filenames[i][-7:-4] + '_Polarisation']}
                self.create_band(src, dst)

                self.dataset.FlushCache()

        for i, elem in enumerate(filenames):
            if filenames[i][-3:] == 'SBI':
                # Calculate sigma0 scaling factor
                Rref = float(gdalMetadata['Reference_Slant_Range'])
                Rexp = float(gdalMetadata['Reference_Slant_Range_Exponent'])
                alphaRef = float(gdalMetadata['Reference_Incidence_Angle'])
                F = float(gdalMetadata['Rescaling_Factor'])
                K = float(gdalMetadata[filenames[i][-7:-4] +
                          '_Calibration_Constant'])
                Ftot = Rref**(2.*Rexp)
                Ftot *= np.sin(alphaRef*np.pi / 180.0)
                Ftot /= F**2.
                Ftot /= K

                #print Ftot

                src = [{'SourceFilename': self.filename,
                        'DataType': gdal.GDT_Float32,
                        'SourceBand': 2*i+1,
                        'ScaleRatio': np.sqrt(Ftot)},
                       {'SourceFilename': self.filename,
                        'DataType': gdal.GDT_Float32,
                        'SourceBand': 2*i+2,
                        'ScaleRatio': np.sqrt(Ftot)}]
                dst = {'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                       'PixelFunctionType': 'RawcountsToSigma0_CosmoSkymed_SBI',
                       'polarisation': gdalMetadata[filenames[i][-7:-4] +
                                                    '_Polarisation'],
                       'name': 'sigma0_%s' % gdalMetadata[filenames[i][-7:-4] +
                                                          '_Polarisation'],
                       'SatelliteID': gdalMetadata['Satellite_ID'],
                       'dataType': gdal.GDT_Float32}
                       #'pass': gdalMetadata['']
                       #         - I can't find this in the metadata...

                self.create_band(src, dst)


                self.dataset.FlushCache()

        self.dataset.SetMetadataItem('time_coverage_start',
                    parse_time(gdalMetadata['Scene_Sensing_Start_UTC']).isoformat())
        self.dataset.SetMetadataItem('time_coverage_end',
                    parse_time(gdalMetadata['Scene_Sensing_Stop_UTC']).isoformat())
