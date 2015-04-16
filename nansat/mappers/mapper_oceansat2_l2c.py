#-------------------------------------------------------------------------------
# Name:     mapper_oceansat2_l2c.py
# Purpose:
#
# Author:       Morten Wergeland Hansen
# Modified: Morten Wergeland Hansen
#
# Created:  13.02.2015
# Last modified:13.02.2015 11:38
# Copyright:    (c) NERSC
# License:
#-------------------------------------------------------------------------------
import datetime
from osgeo import gdal, osr
from nansat.tools import WrongMapperError
from nansat.nsr import NSR
from nansat.vrt import VRT

class Mapper(VRT):

    def __init__(self, fileName, gdalDataset, gdalMetadata, **kwargs):
        if 'Title' not in gdalMetadata.keys():
            raise WrongMapperError
        if not 'Oceansat OCM2 Level-2C' in gdalMetadata['Title']:
            raise WrongMapperError

        subDatasets = gdalDataset.GetSubDatasets()
        for subDataset in subDatasets:
            if 'clo' in subDataset[1]:
                gdalDataset = gdal.Open(subDataset[0])
                gdalMetadata = gdalDataset.GetMetadata()
                break
        lons = [float(gdalMetadata['Start Center Longitude']),
                float(gdalMetadata['Upper Left Longitude']),
                float(gdalMetadata['Upper Right Longitude']),
                float(gdalMetadata['Lower Left Longitude']),
                float(gdalMetadata['Lower Right Longitude']),
                ]
        lats = [float(gdalMetadata['Start Center Latitude']),
                float(gdalMetadata['Upper Left Latitude']),
                float(gdalMetadata['Upper Right Latitude']),
                float(gdalMetadata['Lower Left Latitude']),
                float(gdalMetadata['Lower Right Latitude']),
                ]

        srcSRS = NSR(4326)
        dstSRS = NSR('+proj=lcc +lat_1=5 +lon_0=80')
        transformer = osr.CoordinateTransformation(srcSRS, dstSRS)

        xs, ys = [], []
        for lon,lat in zip(lons, lats):
            x, y, z = transformer.TransformPoint(lon, lat, 0)
            xs.append(x)
            ys.append(y)

        pWidth = (xs[2] - xs[3]) / gdalDataset.RasterXSize
        pHeight = (ys[3] - ys[2]) / gdalDataset.RasterYSize

        srcGeoTransform = (
                xs[1], pWidth, 0,
                ys[1], 0, pHeight,
            )

        # initiate VRT for the LCC projection
        VRT.__init__(self,
                     srcGeoTransform=srcGeoTransform,
                     srcProjection=dstSRS.wkt,
                     srcRasterXSize=gdalDataset.RasterXSize,
                     srcRasterYSize=gdalDataset.RasterYSize)

        metaDict = [{
            'src': {'SourceFilename': subDataset[0], 'SourceBand': 1},
            'dst': {'wkv': 'mass_concentration_of_chlorophyll_a_in_sea_water'}
            }]

        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        year1 = int(gdalMetadata['Start Time'][:4])
        doy1 = int(gdalMetadata['Start Time'][4:7])
        date1 = datetime.datetime(year1, 1,1) + datetime.timedelta(doy1)

        year2 = int(gdalMetadata['End Time'][:4])
        doy2 = int(gdalMetadata['End Time'][4:7])
        date2 = datetime.datetime(year1, 1,1) + datetime.timedelta(doy1 + 0.5)

        self.dataset.SetMetadataItem('start_date', date1.isoformat())
        self.dataset.SetMetadataItem('stop_date', date2.isoformat())
        self.dataset.SetMetadataItem('sensor', 'OCM')
        self.dataset.SetMetadataItem('satellite', 'Oceansat2')
        self.dataset.SetMetadataItem('mapper', 'oceansat2_l2c')
