# Name:    nansat.py
# Purpose: Container of VRT classes
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2018
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
import os
import tempfile
from string import Template, ascii_uppercase, digits
from random import choice
import warnings
import pythesint as pti

import numpy as np

from nansat.node import Node
from nansat.nsr import NSR
from nansat.geolocation import Geolocation
from nansat.tools import add_logger, gdal, osr, OptionError, numpy_to_gdal_type, gdal_type_to_offset

# TODO: Think which variables we should rename

# TODO: Think which conventional names we should use (lon, lat - OK), vrt - ?

class VRT(object):
    """Wrapper around GDAL VRT-file

    The GDAL VRT-file is an XML-file. It contains all metadata, geo-reference
    information and information ABOUT each band including band metadata,
    reference to the bands in the source file.
    VRT-class performs all operation on VRT-files: create, copy, modify,
    read, write, add band, add GeoTransform, set Projection, etc. It uses
    either GDAL methods for these operations (e.g. Create, AddBand,
    SetMetadata, AutoCreateWarpedVRT, etc.) or reads/writes the XML-file
    directly (e.g. remove_geotransform, get_warped_vrt, etc).

    The core of the VRT object is GDAL dataset <self.dataset> generated
    by the GDAL VRT-Driver. The respective VRT-file is located in /vismem
    and has a random name.

    GDAL data model doesn't have place for geolocaion arrays therefore
    VRT-object has instance of Geolocation (self.geolocation)
    an object to keep information about Geolocation metadata:
    reference to file with source data, pixel and line step and offset, etc.

    Domain has an instance of VRT-class <self.vrt>. It keeps only geo-
    reference information.

    All Mappers inherit from VRT. When Nansat opens a file it loops through
    list of mappers, selects the one appropriate for the input file,
    and creates an instance of Mapper.

    Nansat has one instances of Mapper-class (>=VRT-class): self.vrt.
    It holds VRT-file in original projection (derived from the
    input file). After most of the operations with Nansat object
    (e.g. reproject, crop, resize, add_band) self.vrt is replaced with a new
    VRT object which has reference to the previous VRT object inside (self.vrt.vrt).

    """
    COMPLEX_SOURCE_XML = Template('''
            <$SourceType>
                <SourceFilename relativeToVRT="0">$Dataset</SourceFilename>
                <SourceBand>$SourceBand</SourceBand>
                <NODATA>$NODATA</NODATA>
                <ScaleOffset>$ScaleOffset</ScaleOffset>
                <ScaleRatio>$ScaleRatio</ScaleRatio>
                <LUT>$LUT</LUT>
                <SrcRect xOff="$xOff" yOff="$yOff" xSize="$xSize" ySize="$ySize"/>
                <DstRect xOff="0" yOff="0" xSize="$xSize" ySize="$ySize"/>
            </$SourceType> ''')

    RAW_RASTER_BAND_SOURCE_XML = Template('''
            <VRTDataset rasterXSize="$XSize" rasterYSize="$YSize">
              <VRTRasterBand dataType="$DataType"
                band="$BandNum" subClass="VRTRawRasterBand">
                <SourceFilename relativeToVRT="0">$SrcFileName</SourceFilename>
                <ImageOffset>0</ImageOffset>
                <PixelOffset>$PixelOffset</PixelOffset>
                <LineOffset>$LineOffset</LineOffset>
              </VRTRasterBand>
            </VRTDataset> ''')

    REPROJECT_TRANSFORMER = Template('''
        <ReprojectTransformer>
          <ReprojectionTransformer>
            <SourceSRS>$SourceSRS</SourceSRS>
            <TargetSRS>$TargetSRS</TargetSRS>
          </ReprojectionTransformer>
        </ReprojectTransformer> ''')

    filename = ''
    vrt = None

    @classmethod
    def from_gdal_dataset(cls, gdal_dataset, **kwargs):
        """Create VRT from GDAL Dataset with the same size/georeference but wihout bands.

        Parameters
        ----------
            gdal_dataset : gdal.Dataset
                input GDAL dataset
            **kwargs : dict
                arguments for VRT()
        Returns
        -------
            vrt : VRT

        """
        vrt = cls.__new__(cls)
        vrt._init_from_gdal_dataset(gdal_dataset, **kwargs)
        return vrt

    @classmethod
    def from_dataset_params(cls, x_size, y_size, geo_transform, projection,
                            gcps, gcp_projection, **kwargs):
        """Create VRT from GDAL Dataset parameters

        Create VRT with dataset wihout bands but with size/georeference corresponding to
        input parameters.

        Parameters
        ----------
            x_size : int
                X-size of dataset
            y_size : int
                Y-size of dataset
            geotransform : tuple with 6 floats
                informaton on affine transformtaion
            projection : str
                WKT representation of spatial reference system
            gcps : tuple or list with GDAL GCP objects
                GDAL Ground Control Points
            gcp_projection : str
                WKT representation of GCPs spatial reference system
            **kwargs : dict
                arguments for VRT()
        Returns
        -------
            vrt : VRT

        """
        vrt = cls.__new__(cls)
        vrt._init_from_dataset_params(x_size, y_size, geo_transform, projection,
                         gcps, gcp_projection, **kwargs)
        return vrt

    @classmethod
    def from_array(cls, array, **kwargs):
        """Create VRT from numpy array with dataset wih one band but without georeference.

        Parameters
        ----------
            array : numpy.ndarray
                array with data
            **kwargs : dict
                arguments for VRT()
        Returns
        -------
            vrt : VRT

        """
        vrt = cls.__new__(cls)
        vrt._init_from_array(array, **kwargs)
        return vrt

    @classmethod
    def from_lonlat(cls, lon, lat, **kwargs):
        """Create VRT from longitude, latitude arrays

        Create VRT with dataset without bands but with GEOLOCATION metadata and Geolocation
        object. Geolocation contains 2 2D arrays with lon/lat values given at regular pixel/line
        steps.

        Parameters
        ----------
            lon : numpy.ndarray
                array with longitudes
            lat : numpy.ndarray
                array with latitudes
            **kwargs : dict
                arguments for VRT()
        Returns
        -------
            vrt : VRT

        """
        vrt = cls.__new__(cls)
        vrt._init_from_lonlat(lon, lat)
        return vrt

    @classmethod
    def copy_dataset(cls, gdal_dataset, **kwargs):
        """Create VRT with bands and georefernce as a full copy of input GDAL Dataset
        Parameters
        ----------
            gdal_dataset : GDAL.Dataset
                input dataset
            **kwargs : dict
                arguments for VRT()
        Returns
        -------
            vrt : VRT

        """
        vrt = cls.__new__(cls)
        vrt._copy_from_dataset(gdal_dataset, **kwargs)
        return vrt

    def __init__(self, x_size=1, y_size=1, metadata=None, nomem=False, **kwargs):
        """Init VRT object with all attributes
        Parameters
        ----------
            x_size : int
                width of self.dataset
            y_size : int
                arguments for VRT()
            metadata : dict
                dictionray with metadata keys (str) and values (str)
            nomem : bool
                don't create VRT in VSI memory?
        Modifies
        -------
            adds self.logger, self.driver, self.filename, self.band_vrts, self.tps, self.vrt
            adds self.dataset - GDAL Dataset without bands and with size=(x_zie, y_size)
            adds metadata to self.dataset
            writes VRT file content to self.filename

        """
        # essential attributes
        self.logger = add_logger('Nansat')
        self.driver = gdal.GetDriverByName('VRT')
        self.filename = VRT._make_filename(nomem=nomem)
        self.band_vrts = dict()
        self.tps = False
        self.vrt = None

        # create dataset
        self.dataset = self.driver.Create(self.filename, x_size, y_size, bands=0)
        if isinstance(metadata, dict):
            self.dataset.SetMetadata(metadata)
        self.dataset.FlushCache()

    def _init_from_gdal_dataset(self, gdal_dataset, **kwargs):
        """Init VRT from GDAL Dataset with the same size/georeference but wihout bands.

        Parameters
        ----------
            gdal_dataset : gdal.Dataset
                input GDAL dataset
            **kwargs : dict
                arguments for VRT()
        Modifies
        -------
            self - adds all VRT attributes
            self.dataset - sets size and georeference

        """

        # set dataset geo-metadata
        VRT.__init__(self, gdal_dataset.RasterXSize, gdal_dataset.RasterYSize, **kwargs)
        self.dataset.SetGCPs(gdal_dataset.GetGCPs(), gdal_dataset.GetGCPProjection())
        self.dataset.SetProjection(gdal_dataset.GetProjection())
        self.dataset.SetGeoTransform(gdal_dataset.GetGeoTransform())
        metadata = gdal_dataset.GetMetadata()
        for key in metadata:
            self.dataset.SetMetadataItem(key, metadata[key])
        self._add_geolocation(Geolocation.from_dataset(gdal_dataset))
        self.dataset.SetMetadataItem('filename', self.filename)

        # write XML file contents
        self.dataset.FlushCache()

    def _init_from_dataset_params(self, x_size, y_size, geo_transform, projection,
                                    gcps=None, gcp_projection='', **kwargs):
        """Init VRT from GDAL Dataset parameters

        Init VRT with dataset wihout bands but with size/georeference corresponding to
        input parameters.

        Parameters
        ----------
            x_size : int
                X-size of dataset
            y_size : int
                Y-size of dataset
            geotransform : tuple with 6 floats
                informaton on affine transformtaion
            projection : str
                WKT representation of spatial reference system
            gcps : tuple or list with GDAL GCP objects
                GDAL Ground Control Points
            gcp_projection : str
                WKT representation of GCPs spatial reference system
            **kwargs : dict
                arguments for VRT()
        Modifes
        -------
            self - adds all VRT attributes
            self.dataset - sets size and georeference

        """
        VRT.__init__(self, x_size, y_size, **kwargs)
        # set dataset (geo-)metadata
        self.dataset.SetProjection(projection)
        self.dataset.SetGeoTransform(geo_transform)
        if isinstance(gcps, (list, tuple)):
            self.dataset.SetGCPs(gcps, gcp_projection)
        self.dataset.SetMetadataItem('filename', self.filename)

        # write file contents
        self.dataset.FlushCache()

    def _init_from_array(self, array, **kwargs):
        """Init VRT from numpy array with dataset wih one band but without georeference.

        Write contents of the array into flat binary file (VSI)
        Write VRT file with RawRastesrBand, which points to the binary file
        Open the VRT file as self.dataset with GDAL

        Parameters
        ----------
            array : numpy.ndarray
                array with data
            **kwargs : dict
                arguments for VRT()

        Modifies
        ---------
            binary file is written (VSI)
            VRT file is written (VSI)
            self - adds all VRT attributes
            self.dataset is updated

        """
        VRT.__init__(self, **kwargs)
        # create flat binary file (in VSI) from numpy array
        array_type = array.dtype.name
        array_shape = array.shape
        binary_file = self.filename.replace('.vrt', '.raw')
        ofile = gdal.VSIFOpenL(binary_file, 'wb')
        gdal.VSIFWriteL(array.tostring(), len(array.tostring()), 1, ofile)
        gdal.VSIFCloseL(ofile)
        array = None

        # convert Numpy datatype to gdal datatype and pixel offset
        gdal_data_type = numpy_to_gdal_type[array_type]
        pixel_offset = gdal_type_to_offset[gdal_data_type]

        # create XML contents of VRT-file
        line_offset = str(int(pixel_offset) * array_shape[1])
        contents = self.RAW_RASTER_BAND_SOURCE_XML.substitute(
            XSize=array_shape[1],
            YSize=array_shape[0],
            DataType=gdal_data_type,
            BandNum=1,
            SrcFileName=binary_file,
            PixelOffset=pixel_offset,
            LineOffset=line_offset)

        # write XML contents to VRT-file
        self.write_xml(contents)
        self.dataset.SetMetadataItem('filename', self.filename)
        self.dataset.FlushCache()

    def _init_from_lonlat(self, lon, lat, **kwargs):
        """Init VRT from longitude, latitude arrays

        Init VRT with dataset without bands but with GEOLOCATION metadata and Geolocation
        object. Geolocation contains 2 2D arrays with lon/lat values given at regular pixel/line
        steps.

        Parameters
        ----------
            lon : numpy.ndarray
                array with longitudes
            lat : numpy.ndarray
                array with latitudes
            **kwargs : dict
                arguments for VRT() and VRT._lonlat2gcps
        Modifes
        -------
            self - adds all VRT attributes
            self.dataset - sets size and georeference
            self.geolocation - add Geolocation object with all attributes

        """
        VRT.__init__(self, lon.shape[1], lon.shape[0], **kwargs)
        self.dataset.SetGCPs(VRT._lonlat2gcps(lon, lat, **kwargs), NSR().wkt)
        self._add_geolocation(Geolocation(VRT.from_array(lon), VRT.from_array(lat)))
        self.dataset.SetMetadataItem('filename', self.filename)
        self.dataset.FlushCache()

    def _copy_from_dataset(self, gdal_dataset, **kwargs):
        """Init VRT with bands and georefernce as a full copy of input GDAL Dataset
        Parameters
        ----------
            gdal_dataset : GDAL.Dataset
                input dataset
            **kwargs : dict
                arguments for VRT()
        Modifes
        -------
            self - adds all VRT attributes
            self.dataset - sets size and georeference

        """
        # set dataset geo-metadata
        VRT.__init__(self, gdal_dataset.RasterXSize, gdal_dataset.RasterYSize, **kwargs)
        self.dataset = self.driver.CreateCopy(self.filename, gdal_dataset)
        self.dataset.SetMetadataItem('filename', self.filename)

        # write XMl file contents
        self.dataset.FlushCache()

    def __del__(self):
        """Destructor deletes VRT and RAW files"""
        self.dataset = None
        gdal.Unlink(self.filename)
        gdal.Unlink(self.filename.replace('vrt', 'raw'))

    def __repr__(self):
        str_out = os.path.split(self.filename)[1]
        if self.vrt is not None:
            str_out += '=>%s' % self.vrt.__repr__()
        return str_out

    def _create_bands(self, metadata_dict):
        """ Generic function called from the mappers to create bands
        in the VRT dataset from an input dictionary of metadata

        Parameters
        ----------
        metadata_dict : list of dict with params of input bands and generated bands.
            Each dict has:
                'src' : dictionary with parameters of the sources:
                'dst' : dictionary with parameters of the generated bands

        Modifies
        ---------
        Adds bands to the self.dataset based on info in metaDict

        See Also
        ---------
        VRT._create_band()

        """
        for band_dict in metadata_dict:
            src = band_dict['src']
            dst = band_dict.get('dst', None)
            self._create_band(src, dst)
            self.logger.debug('Creating band - OK!')
        self.dataset.FlushCache()

    def _create_band(self, src, dst=None):
        """ Add band to self.dataset:

        Get parameters of the source band(s) from input
        Generate source XML for the VRT, add options of creating
        Call GDALDataset.AddBand
        Set source and options
        Add metadata

        Parameters
        ----------
        src : dict with parameters of sources:
            SourceFilename
            SourceBand
            ScaleRatio
            ScaleOffset
            NODATA
            LUT
            SourceType
            DataType
            ImageOffset (RawVRT)
            PixelOffset (RawVRT)
            LineOffset (RawVRT)
            ByteOrder (RawVRT)
            xSize
            ySize
        dst : dict with parameters of the created band
            name
            dataType
            wkv
            suffix
            AnyOtherMetadata
            PixelFunctionType: - band will be a pixel function defined by the
                                 corresponding name/value.
                                 In this case src may be list of
                                 dicts with parameters for each source.
                               - in case the dst band has a different datatype
                                 than the source band it is important to add a
                                 SourceTransferType parameter in dst
            SourceTransferType

        Returns
        --------
        name : string, name of the added band

        Examples
        --------
        vrt._create_band({'SourceFilename': filename, 'SourceBand': 1})
        vrt._create_band({'SourceFilename': filename, 'SourceBand': 2,
                          'ScaleRatio': 0.0001},
                         {'name': 'LAT', 'wkv': 'latitude'})
        vrt._create_band({'SourceFilename': filename, 'SourceBand': 2},
                         {'suffix': '670',
                          'wkv': 'brightness_temperature'})
        vrt._create_band([{'SourceFilename': filename, 'SourceBand': 1},
                          {'SourceFilename': filename, 'SourceBand': 1}],
                         {'PixelFunctionType': 'NameOfPixelFunction'})

        """
        self.logger.debug('INPUTS: %s, %s " ' % (str(src), str(dst)))
        # Make sure src is list, ready for loop
        if type(src) == dict:
            srcs = [src]
        elif type(src) in [list, tuple]:
            srcs = src
        else:
            raise AttributeError('Wrong src type (%s)! Should be dict or list/tuple of dict'%type(src))

        # Check if dst is given, or create empty dict
        if dst is None:
            dst = {}

        srcs = list(map(VRT._make_source_bands_xml, srcs))

        options = VRT._set_add_band_options(srcs, dst)


# TODO:
#   check individual dst keys in idividual functions

        # set destination dataType (if not given in input parameters)
        if 'dataType' not in dst:
            if (len(srcs) > 1 or float(srcs[0]['ScaleRatio']) != 1.0 or
                    len(srcs[0]['LUT']) > 0 or 'DataType' not in srcs[0]):
                # if pixel function
                # if scaling is applied
                # if LUT
                # if source band not available: float32
                dst['dataType'] = gdal.GDT_Float32
            else:
                # otherwise take the DataType from source
                dst['dataType'] = srcs[0]['DataType']

        # get metadata from WKV using PyThesInt
        wkv_exists = False
        if 'wkv' in dst:
            try:
                wkv = pti.get_wkv_variable(dst['wkv'])
            except IndexError:
                pass
            else:
                wkv_exists = True

        # join wkv[short_name] and dst[suffix] if both given
        if ('name' not in dst and wkv_exists):
            if 'suffix' in dst:
                dstSuffix =  '_' + dst['suffix']
            else:
                dstSuffix = ''
            dst['name'] = wkv['short_name'] + dstSuffix

        # create list of available bands (to prevent duplicate names)
        band_names = []
        for i in range(self.dataset.RasterCount):
            band_names.append(self.dataset.GetRasterBand(i + 1).
                             GetMetadataItem('name'))

        # if name is not given add 'band_00N'
        if 'name' not in dst:
            for n in range(999):
                band_name = 'band_%03d' % n
                if band_name not in band_names:
                    dst['name'] = band_name
                    break
        # if name already exist add '_00N'
        elif dst['name'] in band_names:
            for n in range(999):
                band_name = dst['name'] + '_%03d' % n
                if band_name not in band_names:
                    dst['name'] = band_name
                    break

        self.logger.debug('dst[name]:%s' % dst['name'])

        # Add Band
        self.dataset.AddBand(int(dst['dataType']), options=options)
        dst_raster_band = self.dataset.GetRasterBand(self.dataset.RasterCount)

        # Append sources to destination dataset
        if len(srcs) == 1 and srcs[0]['SourceBand'] > 0:
            # only one source
            dst_raster_band.SetMetadataItem('source_0', str(srcs[0]['XML']), 'new_vrt_sources')
        elif len(srcs) > 1:
            # several sources for PixelFunction
            metadataSRC = {}
            for i, src in enumerate(srcs):
                metadataSRC['source_%d' % i] = src['XML']
            dst_raster_band.SetMetadata(metadataSRC, 'vrt_sources')

        # set metadata from WKV
        if wkv_exists:
            dst_raster_band = VRT._put_metadata(dst_raster_band, wkv)

        # set metadata from provided parameters
        # remove and add params
        dst['SourceFilename'] = srcs[0]['SourceFilename']
        dst['SourceBand'] = str(srcs[0]['SourceBand'])
        dst_raster_band = VRT._put_metadata(dst_raster_band, dst)

        # return name of the created band
        return dst['name']

    def _create_complex_bands(self, filenames):
        """Create bands with complex data type bands with real and imag components
        Parameters
        ---------
            filenames : list of used filenames

        Modifies
        --------
            self.dataset - adds complex bands; removes real and imag bands

        """
        rm_bands = []
        # loop to find real data band
        for i in range(self.dataset.RasterCount):
            band = self.dataset.GetRasterBand(i + 1)
            band_name = band.GetMetadataItem('name')
            if band_name.endswith('_real'):
                real_band_no = i
                real_band_type = band.GetMetadataItem('DataType')
                complex_band_name = band_name.replace('_real', '')
                # loop to find imag data band
                for j in range(self.dataset.RasterCount):
                    band = self.dataset.GetRasterBand(j + 1)
                    band_name = band.GetMetadataItem('name')
                    # find an imaginary data band corresponding to the real
                    # data band and create complex data band from the bands
                    if band_name == complex_band_name + '_imag':
                        imag_band_no = j
                        imag_band_type = band.GetMetadataItem('DataType')
                        dst = band.GetMetadata()
                        dst['name'] = complex_band_name
                        dst['PixelFunctionType'] = 'ComplexData'
                        dst['dataType'] = 10
                        src = [{'SourceFilename': filenames[real_band_no],
                                'SourceBand':  1,
                                'DataType': real_band_type},
                               {'SourceFilename': filenames[imag_band_no],
                                'SourceBand': 1,
                                'DataType': imag_band_type}]
                        self._create_band(src, dst)
                        self.dataset.FlushCache()
                        rm_bands.append(real_band_no + 1)
                        rm_bands.append(imag_band_no + 1)
                        break

        # Delete real and imaginary bands
        if len(rm_bands) != 0:
            self.delete_bands(rm_bands)

    def _add_swath_mask_band(self):
        """ Create a new band where all values = 1

        Modifies
        ---------
        Single band 'swathmask' with ones is added to the self.dataset

        """
        self._create_band(
            src=[{
                'SourceFilename': self.filename,
                'SourceBand':  1,
                'DataType': gdal.GDT_Byte}],
            dst={
                'dataType': gdal.GDT_Byte,
                'wkv': 'swath_binary_mask',
                'PixelFunctionType': 'OnesPixelFunc',
            })

    def _remove_strings_in_metadata_keys(self, gdal_metadata):
        if not gdal_metadata:
            raise WrongMapperError

        for key in gdal_metadata.keys():
            newkey = key.replace('NC_GLOBAL#', '')
            gdal_metadata[newkey] = gdal_metadata.pop(key)

        return gdal_metadata

    def _add_geolocation(self, geolocation):
        """ Add GEOLOCATION to the VRT

        Parameters
        -----------
            geolocation: Geolocation
                with grids of X/Y coordinates

        Modifes
        --------
        Add geolocation to self
        Sets GEOLOCATION metadata

        """
        self.geolocation = geolocation
        self.dataset.SetMetadata(geolocation.data, 'GEOLOCATION')
        self.dataset.FlushCache()

    def _remove_geolocation(self):
        """ Remove GEOLOCATION from the VRT

        Modifes
        --------
        Set self.geolocationArray to None
        Sets GEOLOCATION metadata to ''

        """
        self.geolocation.data = dict()

        # add GEOLOCATION metadata (empty if geolocation is empty)
        self.dataset.SetMetadata('', 'GEOLOCATION')
        self.dataset.FlushCache()

    def _remove_geotransform(self):
        """Remove GeoTransfomr from VRT Object

        Modifies
        ---------
        The tag <GeoTransform> is revoved from the VRT-file

        """
        # read XML content from VRT
        tmp_vrt_xml = self.xml
        # find and remove GeoTransform
        node0 = Node.create(tmp_vrt_xml)
        node0.delNode('GeoTransform')
        # Write the modified elemements back into temporary VRT
        self.write_xml(node0.rawxml())

    def _create_fake_gcps(self, dst_gcps, dst_srs, skip_gcps):
        """Create GCPs with reference self.pixel/line ==> dst.pixel/line

        GCPs from a destination image (dst_gcps) are converted to a gcp of source
        image (src_gcps) this way:

        srcGCPPixel = srcPixel
        srcGCPLine = srcLine
        srcGCPX = dstGCPPixel = f(srcSRS, dstGCPX, dstGCPY)
        srcGCPY = dstGCPLine = f(srcSRS, dstGCPX, dstGCPY)

        Parameters
        -----------
        gcps : list
            GDAL GCPs
        skip_gcps : int
            See nansat.reproject() for explanation

        Returns
        --------
        gcps : dict
            {'gcps': list with GDAL GCPs, 'srs': fake stereo WKT}

        """
        # create transformer. converts lat/lon to pixel/line of SRC image
        src_transformer = gdal.Transformer(self.dataset, None,
                                          ['SRC_SRS=' + self.get_projection(),
                                           'DST_SRS=' + dst_srs.wkt])

        # create 'fake' GCPs
        fake_gcps = []
        for g in dst_gcps[::skip_gcps]:
            # transform DST lat/lon to SRC pixel/line
            succ, point = src_transformer.TransformPoint(1, g.GCPX, g.GCPY)
            srcPixel = point[0]
            srcLine = point[1]

            # swap coordinates in GCPs:
            # pix1/line1 -> lat/lon  =>=>  pix2/line2 -> pix1/line1
            fake_gcps.append(gdal.GCP(g.GCPPixel, g.GCPLine,
                                     0, srcPixel, srcLine))

        return {'gcps': fake_gcps, 'srs': NSR('+proj=stere').wkt}

    def copy(self):
        """Create and return a full copy of a VRT instance"""
        if self.dataset.RasterCount == 0:
            vrt = VRT.from_gdal_dataset(self.dataset)
        else:
            vrt = VRT.copy_dataset(self.dataset)

        vrt.band_vrts = dict(self.band_vrts)
        vrt.tps = bool(self.tps)

        # recursive copy of vrt.vrt
        if self.vrt is not None:
            vrt.vrt = self.vrt.copy()
            # make reference from the new vrt to the copy of vrt.vrt
            new_vrt_xml = vrt.xml.replace(os.path.split(self.vrt.filename)[1],
                                          os.path.split(vrt.vrt.filename)[1])
            vrt.write_xml(new_vrt_xml)
        return vrt

    @property
    def xml(self):
        """Read XML content of the VRT-file using VSI

        Returns
        --------
        string : XMl Content which is read from the VSI file
        """
        self.dataset.FlushCache()
        return VRT.read_vsi(self.filename)

    def write_xml(self, vsi_file_content=None):
        """Write XML content into a VRT dataset

        Parameters
        -----------
        vsiFileContent: string, optional
            XML Content of the VSI file to write

        Modifies
        ---------
        self.dataset
            If XML content was written, self.dataset is re-opened

        """
        vsi_file = gdal.VSIFOpenL(self.filename, 'w')
        gdal.VSIFWriteL(vsi_file_content, len(vsi_file_content), 1, vsi_file)
        gdal.VSIFCloseL(vsi_file)
        # re-open self.dataset with new content
        self.dataset = gdal.Open(self.filename)

    def export(self, filename):
        """Export VRT file as XML into given <filename>"""
        self.driver.CreateCopy(filename, self.dataset)

# TODO:
#   split superfunctional get_warped_vrt into more specific methods
#       get_warped_vrt_geotransform
#       get_warped_vrt_gcp
#       etc...
#   check actual usage of get_warped_vrt and limit input params to the actually used ones only
#   replace keywork arguments with required arguments where possible

    def get_warped_vrt(self, dst_srs=None, resample_alg=0,
                       x_size=0, y_size=0, block_size=None,
                       geo_transform=None, working_data_type=None,
                       use_geolocation=True,
                       use_gcps=True, skip_gcps=1,
                       use_geotransform=True,
                       dst_gcps=[], dst_geolocation=None):

        """Create VRT object with WarpedVRT

        Modifies the input VRT according to the input options
        Creates simple WarpedVRT with AutoCreateWarpedVRT
        Modifies the WarpedVRT according to the input options
        The original VRT is stored as WarpedVRT.vrt

        The function tries to use geolocation array by default;
        if not present (or canceled) tries to use GCPs;
        if not present (or canceled) tries to use GeoTransform
        (either from input dataset or calculates a new one with dx=1,dy=-1).
        Three switches (use_geolocation, use_gcps, use_geotransform)
        allow to select which method to apply for warping. E.g.:
        # #1: srcVRT has Geolocation, geolocation array is used
        warpedVRT = srcVRT.get_warped_vrt(dst_srs, x_size, y_size,
                                             geo_transform)
        # #2: srcVRT has Geolocation, geolocation is not used,
        # either GCPs (if present) or GeoTransform is used
        warpedVRT = srcVRT.get_warped_vrt(dst_srs, x_size, y_size,
                                             geo_transform,
                                             use_geolocation=False)
        # #3: srcVRT has Geolocation or GCPs, geolocation is
        # not used, and GCPs are not used either.
        # Only input GeoTranform is used
        warpedVRT = srcVRT.get_warped_vrt(dst_srs, x_size, y_size,
                                             geo_transform,
                                             use_geolocation=False,
                                             use_gcps=False)

        # #4: srcVRT has whatever georeference, geolocation is not used,
        # GCPs are not used, GeoTransform is not used either.
        # Artificial GeoTranform is calculated: (0, 1, 0, srcVRT.x_size, -1)
        # Warping becomes pure affine resize
        warpedVRT = srcVRT.get_warped_vrt(dst_srs, x_size, y_size,
                                             geo_transform,
                                             use_geolocation=False,
                                             use_gcps=False.,
                                             use_geotransform=false)

        If destination image has GCPs (provided in <dst_gcps>): fake GCPs for
        referencing line/piex of SRC image and X/Y of DST image are created
        and added to the SRC image. After warping dst_gcps are added to
        the WarpedVRT

        If destination image has geolocation (provided in
        <dst_geolocation>):this geolocation is added to the WarpedVRT


        Parameters
        -----------
        dst_srs : string
            WKT of the destination projection
        resample_alg : int (GDALResampleAlg)
            0 : NearestNeighbour,
            1 : Bilinear,
            2 : Cubic,
            3 : CubicSpline,
            4 : Lancoz
        x_size, y_size : int
            width and height of the destination rasetr
        geo_transform : tuple with 6 floats
            destination GDALGeoTransfrom
        dst_gcps : list with GDAL GCPs
            GCPs of the destination image
        dst_geolocation : Geolocation object
            Geolocation of the destination object
        use_geolocation : bool (True)
            Use geolocation in input dataset (if present) for warping
        use_gcps : Boolean (True)
            Use GCPs in input dataset (if present) for warping
        skip_gcps : int
            See nansat.reproject() for explanation
        use_geotransform : Boolean (True)
            Use GeoTransform in input dataset for warping or make artificial
            GeoTransform : (0, 1, 0, srcVRT.x_size, -1)

        Returns
        --------
        warpedVRT : VRT object with WarpedVRT

        """

        # VRT to be warped
        srcVRT = self.copy()

        # srs to be used in AutoCreateWarpedVRT
        acwvSRS = dst_srs

        # if destination GCPs are given: create and add fake GCPs to src
        if len(dst_gcps) > 0 and use_gcps:
            fake_gcps = srcVRT._create_fake_gcps(dst_gcps, NSR(dst_srs), skip_gcps)
            srcVRT.dataset.SetGCPs(fake_gcps['gcps'], fake_gcps['srs'])
            # don't use geolocation
            use_geolocation = False
            acwvSRS = None

        # prepare VRT.dataset for warping.
        # Select if GEOLOCATION,
        # or GCPs, or GeoTransform from the original
        # dataset are used
        if len(self.geolocation.data) > 0 and use_geolocation:
            # use GEOLOCATION by default
            # (remove GCP and GeoTransform)
            srcVRT.dataset.SetGCPs([], '')
            srcVRT._remove_geotransform()
        elif len(srcVRT.dataset.GetGCPs()) > 0 and use_gcps:
            # fallback to GCPs
            # (remove Geolocation and GeoTransform)
            srcVRT.dataset.SetMetadata('', 'GEOLOCATION')
            srcVRT._remove_geotransform()
        elif use_geotransform:
            # fallback to GeoTransform in input VRT
            # (remove Geolocation and GCP)
            srcVRT.dataset.SetMetadata('', 'GEOLOCATION')
            srcVRT.dataset.SetGCPs([], '')
        else:
            # fallback to simplest GeoTransform
            # (remove Geolocation and GCP and replace GeoTransform)
            srcVRT.dataset.SetMetadata('', 'GEOLOCATION')
            srcVRT.dataset.SetGCPs([], '')
            srcVRT.dataset.SetGeoTransform((0, 1, 0, srcVRT.dataset.RasterYSize, 0, -1))
        # create Warped VRT GDAL Dataset
        self.logger.debug('Run AutoCreateWarpedVRT...')
        warped_dataset = gdal.AutoCreateWarpedVRT(srcVRT.dataset, None,
                                             acwvSRS, resample_alg)
        # TODO: implement the below option for proper handling of
        # stereo projections
        # warpedVRT = gdal.AutoCreateWarpedVRT(srcVRT.dataset, '',
        #                                      dst_srs, resample_alg)

        # check if Warped VRT was created
        if warped_dataset is None:
            raise AttributeError('Cannot create warpedVRT!')

        # create VRT object from Warped VRT GDAL Dataset
        self.logger.debug('create VRT object from Warped VRT GDAL Dataset')
        warpedVRT = VRT.copy_dataset(warped_dataset)

        # set x/y size, geo_transform, block_size
        self.logger.debug('set x/y size, geo_transform, block_size')

        # Modify rasterXsize, rasterYsize and geotranforms in the warped VRT
        node0 = Node.create(warpedVRT.xml)

        if x_size > 0:
            node0.replaceAttribute('rasterXSize', str(x_size))
        if y_size > 0:
            node0.replaceAttribute('rasterYSize', str(y_size))

        if geo_transform is not None:
            invGeotransform = gdal.InvGeoTransform(geo_transform)
            # convert proper string style and set to the GeoTransform element
            node0.node('GeoTransform').value = str(geo_transform).strip('()')
            node0.node('DstGeoTransform').value = str(geo_transform).strip('()')
            node0.node('DstInvGeoTransform').value = (
                str(invGeotransform[1]).strip('()'))

            if node0.node('SrcGeoLocTransformer'):
                node0.node('BlockXSize').value = str(x_size)
                node0.node('BlockYSize').value = str(y_size)

            if block_size is not None:
                node0.node('BlockXSize').value = str(block_size)
                node0.node('BlockYSize').value = str(block_size)

            if working_data_type is not None:
                node0.node('WorkingDataType').value = working_data_type

        """
        # TODO: test thoroughly and implement later
        if srcSRS is not None and dst_srs is not None:
            rt = self.REPROJECT_TRANSFORMER.substitute(SourceSRS=None,
                                                      TargetSRS=None)
            print 'rt', rt
            rtNode = Node.create(rt)
            print 'rtNode.xml()', rtNode.xml()
            giptNode = node0.node('GenImgProjTransformer')
            print 'giptNode', giptNode
            giptNode += rtNode
            print 'node0.xml()', node0.xml()
        """
        # overwrite XML of the warped VRT file with uprated size
        # and geotranform
        warpedVRT.write_xml(node0.rawxml())

        # apply thin-spline-transformation option
        if use_gcps and self.tps:
            tmpVRTXML = warpedVRT.xml
            tmpVRTXML = tmpVRTXML.replace('GCPTransformer', 'TPSTransformer')
            warpedVRT.write_xml(tmpVRTXML)

        """
        # TODO: implement the below option for proper handling stereo
        # projections over the pole get source projection from GCPs or
        # from dataset (TODO: or from Geolocation)
        if len(srcVRT.dataset.GetGCPs()) == 0:
            srcSRS = srcVRT.dataset.GetProjection()
        else:
            srcSRS = srcVRT.dataset.GetGCPProjection()
        # modify the VRT XML file
        """

        # if given, add dst GCPs
        self.logger.debug('if given, add dst GCPs')
        if len(dst_gcps) > 0:
            warpedVRT.dataset.SetGCPs(dst_gcps, dst_srs)
            warpedVRT._remove_geotransform()
            warpedVRT.dataset.SetProjection('')

        # if given, add dst Geolocation
        self.logger.debug('# if given, add dst Geolocation')
        if dst_geolocation is not None:
            warpedVRT._remove_geotransform()
            warpedVRT._add_geolocation(dst_geolocation)
            warpedVRT.dataset.SetProjection('')

        # Copy self to warpedVRT
        warpedVRT.vrt = self.copy()

        # replace the reference from srcVRT to self
        self.logger.debug('replace the reference from srcVRT to self')
        rawFileName = str(os.path.basename(warpedVRT.vrt.filename))
        warpedXML = str(warpedVRT.xml)
        node0 = Node.create(warpedXML)
        node1 = node0.node('GDALWarpOptions')
        node1.node('SourceDataset').value = '/vsimem/' + rawFileName
        warpedVRT.write_xml(node0.rawxml())

        return warpedVRT

    def copyproj(self, filename):
        """ Copy geoloctation data from given VRT to a figure file

        Useful for adding geolocation information to figure
        files produced e.g. by Figure class, which contain no geolocation.
        Analogue to utility gdalcopyproj.py.

        Parameters
        -----------
        filename : string
            Name of file to which the geolocation data shall be written

        """
        figDataset = gdal.Open(filename, gdal.GA_Update)
        figDataset.SetGeoTransform(self.dataset.GetGeoTransform())
        figDataset.SetProjection(self.dataset.GetProjection())
        gcps = self.dataset.GetGCPs()
        if len(gcps) != 0:
            figDataset.SetGCPs(gcps, self.dataset.GetGCPProjection())
        figDataset = None  # Close and write output file

    def delete_band(self, band_num):
        """ Delete a band from the given VRT

        Parameters
        ----------
        band_num : int
            band number

        """
        node0 = Node.create(self.xml)
        node0.delNode('VRTRasterBand', options={'band': band_num})
        node0.delNode('BandMapping', options={'src': band_num})
        self.write_xml(node0.rawxml())

    def delete_bands(self, band_nums):
        """ Delete bands

        Parameters
        ----------
        bandNums : list
            elements are int

        """
        band_nums.sort(reverse=True)
        for i in band_nums:
            self.delete_band(i)

    def get_shifted_vrt(self, shift_degree):
        """ Roll data in bands westwards or eastwards

        Create shiftVRT which references self. Modify georeference
        of shiftVRT to account for the roll. Add as many bands as in self
        but for each band create two complex sources: for western
        and eastern parts. Keep self in shiftVRT.vrt

        Parameters
        ----------
        shift_degree : float
            rolling angle, how far east/west to roll

        Returns
        -------
        shiftVRT : VRT object with rolled bands

        """
        # Copy self into self.vrt
        shiftVRT = VRT.from_gdal_dataset(self.dataset)
        shiftVRT.vrt = self.copy()

        if shift_degree < 0:
            shift_degree += 360.0

        geo_transform = shiftVRT.vrt.dataset.GetGeoTransform()
        shiftPixel = int(shift_degree / float(geo_transform[1]))
        geo_transform = list(geo_transform)
        geo_transform[0] = round(geo_transform[0] + shift_degree, 3)
        newEastBorder = geo_transform[0] + (geo_transform[1] *
                                           shiftVRT.dataset.RasterXSize)
        if newEastBorder > 360.0:
            geo_transform[0] -= 360.0
        shiftVRT.dataset.SetGeoTransform(tuple(geo_transform))

        # Add bands to self
        for iBand in range(shiftVRT.vrt.dataset.RasterCount):
            src = {'SourceFilename': shiftVRT.vrt.filename,
                   'SourceBand': iBand + 1}
            dst = shiftVRT.vrt.dataset.GetRasterBand(iBand+1).GetMetadata()
            shiftVRT._create_band(src, dst)

        # read xml and create the node
        XML = shiftVRT.xml
        node0 = Node.create(XML)

        # divide into two bands and switch the bands
        for i in range(len(node0.nodeList('VRTRasterBand'))):
            # create i-th 'VRTRasterBand' node
            node1 = node0.node('VRTRasterBand', i)
            # modify the 1st band
            shiftStr = str(shiftPixel)
            sizeStr = str(shiftVRT.vrt.dataset.RasterXSize - shiftPixel)
            node1.node('ComplexSource').node('DstRect').replaceAttribute('xOff', shiftStr)
            node1.node('ComplexSource').node('DstRect').replaceAttribute('xSize', sizeStr)
            node1.node('ComplexSource').node('SrcRect').replaceAttribute('xSize', sizeStr)

            # add the 2nd band
            xmlSource = node1.rawxml()
            cloneNode = Node.create(xmlSource).node('ComplexSource')
            cloneNode.node('SrcRect').replaceAttribute('xOff', sizeStr)
            cloneNode.node('DstRect').replaceAttribute('xOff', str(0))
            cloneNode.node('SrcRect').replaceAttribute('xSize', shiftStr)
            cloneNode.node('DstRect').replaceAttribute('xSize', shiftStr)

            # get VRTRasterBand with inserted ComplexSource
            node1 = node1.insert(cloneNode.rawxml())
            node0.replaceNode('VRTRasterBand', i, node1)

        # write down XML contents
        shiftVRT.write_xml(node0.rawxml())

        return shiftVRT

    def get_sub_vrt(self, steps=1):
        """Return sub-VRT from given depth

        Iteratively copy self.vrt into self until
        self.vrt is None or steps == 0

        Parameters
        -----------
        steps : int
            How many sub VRTs to restore

        Returns
        -------
        self : if no deeper VRTs found
        self.vrt : if deeper VRTs are found

        Modifies
        --------
        self
        self.vrt

        """

        # check if self is the last valid (deepest) VRT
        if self.vrt is None:
            return self

        # check if required depth of restoration is met
        if steps == 0:
            return self

        # decrease the depth of restoration
        steps -= 1

        # return restored sub-VRT
        return self.vrt.get_sub_vrt(steps)

    def get_super_vrt(self):
        """Create vrt with subVRT

        copy of self in vrt.vrt and change references from vrt to vrt.vrt

        """
        # create new vrt that refers to a copy of self
        superVRT = VRT.from_gdal_dataset(self.dataset)
        superVRT.vrt = self
        superVRT.tps = self.tps

        # add bands to the new vrt
        for iBand in range(superVRT.vrt.dataset.RasterCount):
            src = {'SourceFilename': superVRT.vrt.filename,
                   'SourceBand': iBand + 1}
            dst = superVRT.vrt.dataset.GetRasterBand(iBand + 1).GetMetadata()
            # remove PixelFunctionType from metadata to prevent its application
            if 'PixelFunctionType' in dst:
                dst.pop('PixelFunctionType')
            superVRT._create_band(src, dst)
        superVRT.dataset.FlushCache()

        return superVRT

    def get_subsampled_vrt(self, new_raster_x_size, new_raster_y_size, resample_alg):
        """Create VRT and replace step in the source"""

        subsamVRT = self.get_super_vrt()

        # Get XML content from VRT-file
        node0 = Node.create(subsamVRT.xml)

        # replace rasterXSize in <VRTDataset>
        node0.replaceAttribute('rasterXSize', str(new_raster_x_size))
        node0.replaceAttribute('rasterYSize', str(new_raster_y_size))

        # replace xSize in <DstRect> of each source
        for iNode1 in node0.nodeList('VRTRasterBand'):
            for sourceName in ['ComplexSource', 'SimpleSource']:
                for iNode2 in iNode1.nodeList(sourceName):
                    iNodeDstRect = iNode2.node('DstRect')
                    iNodeDstRect.replaceAttribute('xSize',
                                                  str(new_raster_x_size))
                    iNodeDstRect.replaceAttribute('ySize',
                                                  str(new_raster_y_size))
            # if method=-1, overwrite 'ComplexSource' to 'AveragedSource'
            if resample_alg == -1:
                iNode1.replaceTag('ComplexSource', 'AveragedSource')
                iNode1.replaceTag('SimpleSource', 'AveragedSource')
                # if the values are complex number, give a warning
                if iNode1.getAttribute('dataType').startswith('C'):
                    warnings.warn(
                        'Band %s : The imaginary parts of complex numbers '
                        'are lost when resampling by averaging '
                        '(resample_alg=-1)' % iNode1.getAttribute('band'))

        # Write the modified elemements into VRT
        subsamVRT.write_xml(node0.rawxml())

        return subsamVRT

    def transform_points(self, col_vector, row_vector, dst2src=0,
                         dst_srs=NSR(), dst_ds=None, options=None):
        """Transform given lists of X,Y coordinates into lat/lon

        Parameters
        -----------
        col_vector, row_vector : lists
            X and Y coordinates with any coordinate system
        dst2src : 0 or 1
            1 for inverse transformation, 0 for forward transformation.
        dst_srs : NSR
            destination SRS.
        dst_ds : GDAL Dataset
            destination dataset. The default is None.
            It means transform ownPixLin <--> ownXY.
        option : string
            if 'METHOD=GEOLOC_ARRAY', specify here.

        Returns
        --------
        lon_vector, lat_vector : numpy arrays
            X and Y coordinates in degree of lat/lon

        """
        # get source SRS (either Projection or GCPProjection)
        srcWKT = self.get_projection()

        # prepare options
        if options is None:
            options = ['SRC_SRS=' + srcWKT, 'DST_SRS=' + dst_srs.wkt]
            # add TPS method if we have GCPs and self.tps is True
            if self.tps and len(self.dataset.GetGCPs()) > 0:
                options.append('METHOD=GCP_TPS')

        # create transformer
        transformer = gdal.Transformer(self.dataset, dst_ds, options)

        # convert lists with X,Y coordinates to 2D numpy array
        xy = np.array([col_vector, row_vector]).transpose()

        # transfrom coordinates
        lonlat = transformer.TransformPoints(dst2src, xy)[0]

        # convert return to lon,lat vectors
        lonlat = np.array(lonlat)
        if lonlat.shape[0] > 0:
            lon_vector = lonlat[:, 0]
            lat_vector = lonlat[:, 1]
        else:
            lon_vector, lat_vector = [], []

        return lon_vector, lat_vector

    def get_projection(self):
        """Get projection from self.dataset

        Get projection from GetProjection() or GetGCPProjection().
        If both are empty, raise error

        Return
        -------
        projection : projection or GCPprojection

        Raises
        -------
        ProjectionError : occurs when the projection is empty.

        TODO: see issue #190 in nansat...

        """
        # get projection or GCPProjection
        projection = self.dataset.GetProjection()
        if projection == '':
            projection = self.dataset.GetGCPProjection()

        return projection

    def get_resized_vrt(self, x_size, y_size, resample_alg=1):

        """ Resize VRT

        Create Warped VRT with modified RasterXSize, RasterYSize, GeoTransform.
        The returned VRT object has a copy of its original source VRT in its
        own vrt object (e.g. warpedVRT.vrt = originalVRT.copy()).

        Parameters
        -----------
        x_size, y_size : int
            new size of the VRT object
        resample_alg : GDALResampleAlg
            see also gdal.AutoCreateWarpedVRT

        Returns
        --------
        warped_vrt : Resized VRT object

        """
        # modify GeoTransform: set resolution from new X/Y size
        geo_transform = (0.5,
                        float(self.dataset.RasterXSize - 1.0) / float(x_size),
                        0,
                        self.dataset.RasterYSize - 0.5,
                        0,
                        - float(self.dataset.RasterYSize - 1.0) / float(y_size))

        # update size and GeoTranform in XML of the warped VRT object
        warped_vrt = self.get_warped_vrt(x_size=x_size, y_size=y_size,
                                        geo_transform=geo_transform,
                                        use_geolocation=False,
                                        use_gcps=False,
                                        use_geotransform=False,
                                        resample_alg=resample_alg)

        return warped_vrt

    def reproject_GCPs(self, dst_srs):
        """Reproject all GCPs to a new spatial reference system

        Necessary before warping an image if the given GCPs
        are in a coordinate system which has a singularity
        in (or near) the destination area (e.g. poles for lonlat GCPs)

        Parameters
        ----------
        dst_srs : proj4, WKT, NSR, EPSG
            Destiination SRS given as any NSR input parameter

        Modifies
        --------
            Reprojects all GCPs to new SRS and updates GCPProjection
        """
        # Make tranformer from GCP SRS to destination SRS
        dst_srs = NSR(dst_srs)
        srcSRS = NSR(self.dataset.GetGCPProjection())
        transformer = osr.CoordinateTransformation(srcSRS, dst_srs)

        # Reproject all GCPs
        srcGCPs = self.dataset.GetGCPs()
        dst_gcps = []
        for srcGCP in srcGCPs:
            (x, y, z) = transformer.TransformPoint(srcGCP.GCPX,
                                                   srcGCP.GCPY,
                                                   srcGCP.GCPZ)
            dstGCP = gdal.GCP(x, y, z, srcGCP.GCPPixel,
                              srcGCP.GCPLine, srcGCP.Info, srcGCP.Id)
            dst_gcps.append(dstGCP)

        # Update dataset
        self.dataset.SetGCPs(dst_gcps, dst_srs.wkt)

    @staticmethod
    def read_vsi(filename):
        """Read text from input <filename:str> using VSI and return <content:str>."""
        # open
        vsiFile = gdal.VSIFOpenL(filename, 'r')
        # get file size
        gdal.VSIFSeekL(vsiFile, 0, 2)
        vsiFileSize = gdal.VSIFTellL(vsiFile)
        # fseek to start again
        gdal.VSIFSeekL(vsiFile, 0, 0)
        # read
        vsiFileContent = gdal.VSIFReadL(vsiFileSize, 1, vsiFile)
        gdal.VSIFCloseL(vsiFile)
        return vsiFileContent

    @staticmethod
    def _make_source_bands_xml(src_in):
        """Check parameters of band source, set defaults and generate XML for VRT

        Parameters
        -------
            src_in : dict
                dict with band source parameters (SourceFilename, SourceBand, etc)
        Returns
        -------
            src : dict
                updated dict with XML entry
        """
        # check if SourceFilename is given
        if 'SourceFilename' not in src_in:
            raise AttributeError('SourceFilename not given!')
        # set default values
        src = {'SourceBand': 1,
               'LUT': '',
               'NODATA': '',
               'SourceType': 'ComplexSource',
               'ScaleRatio': 1.0,
               'ScaleOffset': 0.0}
        src.update(src_in)

        # find DataType of source (if not given in src)
        if src['SourceBand'] > 0 and 'DataType' not in src:
            raster_band = gdal.Open(src['SourceFilename']).GetRasterBand(src['SourceBand'])
            src['DataType'] = raster_band.DataType

        if 'xSize' not in src or 'ySize' not in src:
            ds = gdal.Open(src['SourceFilename'])
            src['xSize'] = ds.RasterXSize
            src['ySize'] = ds.RasterYSize

        # TODO:
        #   write XML from dictionary using a standard method (not a filling a template)

        # create XML for each source
        src['XML'] = VRT.COMPLEX_SOURCE_XML.substitute(
            Dataset=src['SourceFilename'],
            SourceBand=src['SourceBand'],
            SourceType=src['SourceType'],
            NODATA=src['NODATA'],
            ScaleOffset=src['ScaleOffset'],
            ScaleRatio=src['ScaleRatio'],
            LUT=src['LUT'],
            xSize=src['xSize'],
            ySize=src['ySize'],
            xOff=src.get('xOff', 0),
            yOff=src.get('yOff', 0),)

        return src

    @staticmethod
    def _set_add_band_options(srcs, dst):
        """Generate options for gdal.AddBand based on input band src and dst parameters"""
        if 'PixelFunctionType' in dst and len(dst['PixelFunctionType']) > 0:
            # in case of PixelFunction
            options = ['subClass=VRTDerivedRasterBand',
                       'PixelFunctionType=%s' % dst['PixelFunctionType']]
            if 'SourceTransferType' in dst:
                options.append('SourceTransferType=%s' % dst['SourceTransferType'])
        elif len(srcs) == 1 and srcs[0]['SourceBand'] == 0:
            # in case of VRTRawRasterBand
            options = ['subclass=VRTRawRasterBand',
                       'SourceFilename=%s' % srcs[0]['SourceFilename'],
                       'ImageOffset=%d' % srcs[0]['ImageOffset'],
                       'PixelOffset=%d' % srcs[0]['PixelOffset'],
                       'LineOffset=%d' % srcs[0]['LineOffset'],
                       'ByteOrder=%s' % srcs[0]['ByteOrder']]
        else:
            # in common case
            options = []

        return options

    @staticmethod
    def _lonlat2gcps(lon, lat, n_gcps=100, **kwargs):
        """ Create list of GCPs from given grids of latitude and longitude

        take <numOfGCPs> regular pixels from inpt <lat> and <lon> grids
        Create GCPs from these pixels
        Create latlong GCPs projection

        Parameters
        -----------
        lat : Numpy grid
            array of latitudes
        lon : Numpy grid
            array of longitudes (should be the same size as lat)
        numOfGCPs : int, optional, default = 100
            number of GCPs to create

        Returns
        --------
        gcsp : List with GDAL GCPs

        """
        # estimate step of GCPs
        gcpSize = np.sqrt(n_gcps)
        step0 = max(1, int(float(lat.shape[0]) / gcpSize))
        step1 = max(1, int(float(lat.shape[1]) / gcpSize))

        # generate list of GCPs
        gcps = []
        k = 0
        for i0 in range(0, lat.shape[0], step0):
            for i1 in range(0, lat.shape[1], step1):
                # create GCP with X,Y,pixel,line from lat/lon matrices
                gcp = gdal.GCP(float(lon[i0, i1]), float(lat[i0, i1]), 0, i1, i0)
                gcps.append(gcp)
                k += 1

        return gcps

    @staticmethod
    def _put_metadata(raster_band, metadata_dict):
        """ Put all metadata into a raster band

        Take metadata from metadataDict and put to the GDAL Raster Band

        Parameters
        ----------
        rasterBand : GDALRasterBand
            destination band without metadata

        metadataDict : dictionary
            keys are names of metadata, values are values

        Returns
        --------
        rasterBand : GDALRasterBand
            destination band with metadata

        """
        for key in metadata_dict:
            try:
                meta_value = str(metadata_dict[key])
                meta_key = str(key)
            except UnicodeEncodeError:
                self.logger.error('Cannot add %s to metadata' % key)
            else:
                raster_band.SetMetadataItem(meta_key, meta_value)

        return raster_band

    @staticmethod
    def _make_filename(extention='vrt', nomem=False):
        """Create random VSI file name

        Parameters
        ----------
        extention : string
            extension of the file

        Returns
        -------
        random file name

        """
        if nomem:
            fd, filename = tempfile.mkstemp(suffix='.'+extention)
            os.close(fd)
        else:
            allChars = ascii_uppercase + digits
            randomChars = ''.join(choice(allChars) for x in range(10))
            filename = '/vsimem/%s.%s' % (randomChars, extention)
        return filename
