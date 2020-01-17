# Name:    vrt.py
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
from __future__ import absolute_import, unicode_literals, division
import os
import tempfile
from string import Template, ascii_uppercase, digits
from random import choice
import warnings
import pythesint as pti

import osr
import gdal
import numpy as np

from nansat.node import Node
from nansat.nsr import NSR
from nansat.geolocation import Geolocation
from nansat.utils import add_logger, numpy_to_gdal_type, gdal_type_to_offset, remove_keys

from nansat.exceptions import NansatProjectionError

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

    Notes
    --------
    adds self.logger, self.driver, self.filename, self.band_vrts, self.tps, self.vrt
    adds self.dataset - GDAL Dataset without bands and with size=(x_zie, y_size)
    adds metadata to self.dataset
    writes VRT file content to self.filename

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

    # instance attributes
    filename = ''
    vrt = None
    dataset = None
    logger = None
    driver = None
    band_vrts = None
    tps = None
    geolocation = None

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
    def from_lonlat(cls, lon, lat, add_gcps=True, **kwargs):
        """Create VRT from longitude, latitude arrays

        Create VRT with dataset without bands but with GEOLOCATION metadata and Geolocation
        object. Geolocation contains 2 2D arrays with lon/lat values given at regular pixel/line
        steps. GCPs can bee created from lon/lat arrays and added to the dataset

        Parameters
        ----------
        lon : numpy.ndarray
            array with longitudes
        lat : numpy.ndarray
            array with latitudes
        add_gcps : bool
            Create GCPs from lon/lat arrays and add to dataset
        **kwargs : dict
            arguments for VRT()

        Returns
        -------
        vrt : VRT

        """
        vrt = cls.__new__(cls)
        vrt._init_from_lonlat(lon, lat, add_gcps, **kwargs)
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
        """Init VRT object with all attributes"""
        if metadata is None:
            metadata = dict()

        # essential attributes
        self.logger = add_logger('Nansat')
        self.driver = gdal.GetDriverByName(str('VRT'))
        self.filename = str(VRT._make_filename(nomem=nomem))
        self.band_vrts = dict()
        self.tps = False
        self.vrt = None

        # create dataset
        self.dataset = self.driver.Create(self.filename, x_size, y_size, bands=0)
        self.dataset.SetMetadata(metadata)
        self.dataset.FlushCache()

    def _init_from_gdal_dataset(self, gdal_dataset, geolocation=None, **kwargs):
        """Init VRT from GDAL Dataset with the same size/georeference but wihout bands/metadata.

        Parameters
        ----------
        gdal_dataset : gdal.Dataset
            input GDAL dataset
        geolocation : geolocation.Geolocation
            optional Geolocation arrays
        **kwargs : dict
            arguments for VRT()

        Notes
        -------
        self - adds all VRT attributes
        self.dataset - sets size and georeference

        """
        # WORKAROUND for not providing metadata explicitly. Should be avoided.
        # get metadata from input and gdal_dataset
        # metadata = kwargs.pop('metadata', dict())
        # metadata.update(gdal_dataset.GetMetadata())
        # set dataset parameters and metadata
        VRT.__init__(self, gdal_dataset.RasterXSize, gdal_dataset.RasterYSize, **kwargs)
        self.dataset.SetGCPs(gdal_dataset.GetGCPs(), gdal_dataset.GetGCPProjection())
        self.dataset.SetProjection(gdal_dataset.GetProjection())
        self.dataset.SetGeoTransform(gdal_dataset.GetGeoTransform())
        if geolocation is None:
            geolocation = Geolocation.from_dataset(gdal_dataset)
        self._add_geolocation(geolocation)
        self.dataset.SetMetadataItem(str('filename'), self.filename)

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

        Notes
        -------
        self - adds all VRT attributes
        self.dataset - sets size and georeference

        """
        # x_size, y_size, geo_transform, projection, gcps=None, gcp_projection='', **kwargs
        VRT.__init__(self, x_size, y_size, **kwargs)
        # set dataset (geo-)metadata
        self.dataset.SetProjection(str(projection))
        self.dataset.SetGeoTransform(geo_transform)
        if isinstance(gcps, (list, tuple)):
            self.dataset.SetGCPs(gcps, str(gcp_projection))
        self.dataset.SetMetadataItem(str('filename'), self.filename)

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

        Notes
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
        ofile = gdal.VSIFOpenL(str(binary_file), str('wb'))
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
        self.dataset.SetMetadataItem(str('filename'), self.filename)
        self.dataset.FlushCache()

    def _init_from_lonlat(self, lon, lat, add_gcps=True, **kwargs):
        """Init VRT from longitude, latitude arrays

        Init VRT with dataset without bands but with GEOLOCATION metadata and Geolocation
        object. Geolocation contains 2 2D arrays with lon/lat values given at regular pixel/line
        steps. GCPs can be created from lon/lat values and added to the dataset.

        Parameters
        ----------
        lon : numpy.ndarray
            array with longitudes
        lat : numpy.ndarray
            array with latitudes
        add_gcps : bool
            Add GCPs to dataset
        **kwargs : dict
            arguments for VRT() and VRT._lonlat2gcps

        Notes
        -------
        self - adds all VRT attributes
        self.dataset - sets size and georeference
        self.geolocation - add Geolocation object with all attributes

        """
        VRT.__init__(self, lon.shape[1], lon.shape[0], **kwargs)
        if add_gcps:
            self.dataset.SetGCPs(VRT._lonlat2gcps(lon, lat, **kwargs), NSR().wkt)
        self._add_geolocation(Geolocation(VRT.from_array(lon), VRT.from_array(lat)))
        self.dataset.SetMetadataItem(str('filename'), self.filename)
        self.dataset.FlushCache()

    def _copy_from_dataset(self, gdal_dataset, geolocation=None, **kwargs):
        """Init VRT with bands and georefernce as a full copy of input GDAL Dataset

        Parameters
        ----------
        gdal_dataset : GDAL.Dataset
            input dataset
        **kwargs : dict
            arguments for VRT()

        Notes
        --------
        self - adds all VRT attributes
        self.dataset - sets size and georeference

        """
        # set dataset geo-metadata
        VRT.__init__(self, gdal_dataset.RasterXSize, gdal_dataset.RasterYSize, **kwargs)
        self.dataset = self.driver.CreateCopy(self.filename, gdal_dataset)
        self.dataset.SetMetadataItem(str('filename'), self.filename)
        if geolocation is None:
            geolocation = Geolocation.from_dataset(gdal_dataset)
        self._add_geolocation(geolocation)
        # write XMl file contents
        self.dataset.FlushCache()

    def __del__(self):
        """Destructor deletes VRT and RAW files"""
        self.dataset = None

        if gdal.VSIStatL(self.filename) is not None:
            gdal.Unlink(self.filename)

        if gdal.VSIStatL(self.filename.replace('vrt', 'raw')) is not None:
            gdal.Unlink(self.filename.replace('vrt', 'raw'))

    def __repr__(self):
        str_out = os.path.split(self.filename)[1]
        if self.vrt is not None:
            str_out += '=>%s' % self.vrt.__repr__()
        return str_out

    def _create_complex_bands(self, filenames):
        """Create bands with complex data type bands with real and imag components
        Parameters
        ---------
            filenames : list of used filenames

        Notes
        --------
        self.dataset - adds complex bands; removes real and imag bands

        """
        rm_bands = []
        # loop to find real data band
        for i in range(self.dataset.RasterCount):
            band = self.dataset.GetRasterBand(i + 1)
            band_name = band.GetMetadataItem(str('name'))
            if band_name.endswith('_real'):
                real_band_no = i
                real_band_type = band.GetMetadataItem(str('DataType'))
                complex_band_name = band_name.replace('_real', '')
                # loop to find imag data band
                for j in range(self.dataset.RasterCount):
                    band = self.dataset.GetRasterBand(j + 1)
                    band_name = band.GetMetadataItem(str('name'))
                    # find an imaginary data band corresponding to the real
                    # data band and create complex data band from the bands
                    if band_name == complex_band_name + '_imag':
                        imag_band_no = j
                        imag_band_type = band.GetMetadataItem(str('DataType'))
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
                        self.create_band(src, dst)
                        self.dataset.FlushCache()
                        rm_bands.append(real_band_no + 1)
                        rm_bands.append(imag_band_no + 1)
                        break

        # Delete real and imaginary bands
        if len(rm_bands) != 0:
            self.delete_bands(rm_bands)

    def _add_swath_mask_band(self):
        """Create a new band where all values = 1

        Notes
        ---------
        Single band 'swathmask' with ones is added to the self.dataset

        """
        self.create_band(
            src=[{
                'SourceFilename': self.filename,
                'SourceBand':  1,
                'DataType': gdal.GDT_Byte}],
            dst={
                'dataType': gdal.GDT_Byte,
                'wkv': 'swath_binary_mask',
                'PixelFunctionType': 'OnesPixelFunc',
            })

    def _add_geolocation(self, geolocation):
        """ Add GEOLOCATION to the VRT

        Parameters
        -----------
            geolocation: Geolocation
                with grids of X/Y coordinates

        Notes
        --------
        Add geolocation to self
        Sets GEOLOCATION metadata

        """
        if geolocation is None:
            return
        self.geolocation = geolocation
        self.dataset.SetMetadata(geolocation.data, str('GEOLOCATION'))
        self.dataset.FlushCache()

    def _remove_geolocation(self):
        """ Remove GEOLOCATION from the VRT

        Notes
        --------
        Set self.geolocationArray to None
        Sets GEOLOCATION metadata to ''

        """
        self.geolocation = None

        # add GEOLOCATION metadata (empty if geolocation is empty)
        self.dataset.SetMetadata(str(''), str('GEOLOCATION'))
        self.dataset.FlushCache()

    def _remove_geotransform(self):
        """Remove GeoTransfomr from VRT Object

        Notes
        ---------
        The tag <GeoTransform> is revoved from the VRT-file

        """
        # read XML content from VRT
        # find and remove GeoTransform
        node0 = Node.create(self.xml)
        node0.delNode('GeoTransform')
        # Write the modified elemements back into temporary VRT
        self.write_xml(node0.rawxml())

    def _set_fake_gcps(self, dst_srs, dst_gcps, skip_gcps):
        """Create GCPs with reference self.pixel/line ==> dst.pixel/line and set to self.dataset

        GCPs from a destination image (dst_gcps) are converted to a gcp of source
        image (src_gcps) this way:

        src_gcpPixel = srcPixel
        src_gcpLine = srcLine
        src_gcpX = dstGCPPixel = f(src_srs, dstGCPX, dstGCPY)
        src_gcpY = dstGCPLine = f(src_srs, dstGCPX, dstGCPY)

        Parameters
        -----------
        dst_srs : str
            destination SRS
        dst_gcps : list or tuple
            GDAL GCPs
        skip_gcps : int
            number of GCPs to skip

        Returns
        --------
        srs : str or None
            if GCPs are set, returns None
            if GCPs are absent, retuns dst_srs

        """
        if len(dst_gcps) == 0:
            return dst_srs
        # create transformer. converts lat/lon to pixel/line of SRC image
        src_transformer = gdal.Transformer(self.dataset, None,
                                          ['SRC_SRS=' + self.get_projection()[0],
                                           'DST_SRS=' + NSR(dst_srs).wkt])

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

        self.dataset.SetGCPs(fake_gcps, NSR('+proj=stere').wkt)
        return None

    def _set_geotransform_for_resize(self):
        """Prepare VRT.dataset for resizing."""
        self.dataset.SetMetadata(str(''), str('GEOLOCATION'))
        self.dataset.SetGCPs([], str(''))
        self.dataset.SetGeoTransform((0, 1, 0, self.dataset.RasterYSize, 0, -1))

    def _set_gcps_geolocation_geotransform(self):
        """Prepare VRT.dataset for warping."""
        # Select if GEOLOCATION, or GCPs, or GeoTransform from the original dataset are used
        if self.geolocation is not None and len(self.geolocation.data) > 0:
            # use GEOLOCATION by default (remove GCP and GeoTransform)
            self.dataset.SetGCPs([], str(''))
            self._remove_geotransform()
        elif len(self.dataset.GetGCPs()) > 0:
            # otherwise fallback to GCPs (remove Geolocation and GeoTransform)
            self.dataset.SetMetadata(str(''), str('GEOLOCATION'))
            self._remove_geotransform()
        else:
            # otherwise fallback to GeoTransform in input VRT (remove Geolocation and GCP)
            self.dataset.SetMetadata(str(''), str('GEOLOCATION'))
            self.dataset.SetGCPs([], str(''))
        self.dataset.FlushCache()

    def _update_warped_vrt_xml(self, x_size, y_size, geo_transform, block_size, working_data_type):
        """Update rasterXsize, rasterYsize, geotransform, block_size and working_data_type"""
        node0 = Node.create(str(self.xml))
        node0.replaceAttribute('rasterXSize', str(x_size))
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

        # overwrite XML file with updated size, geotranform, etc
        self.write_xml(node0.rawxml())

    def _create_band_name(self, dst):
        """Create band name based on destination band dictionary <dst>"""
        band_name = dst.get('name', None)

        # try to get metadata from WKV using PyThesInt if it exists
        wkv_dst = dst.get('wkv', None)
        try:
            wkv_pti = pti.get_wkv_variable(str(wkv_dst))
        except IndexError:
            # IndexError is raised when PyThesInt doesn't find the requested WKV.
            # In that case and empty dict without any metadata is created
            wkv = {}
        else:
            # If WKV was found by PyThesInt, a dict with metadata is created
            wkv = dict(wkv_pti)

        if band_name is None:
            band_name = wkv.get('short_name', 'band')
            if 'suffix' in dst:
                 band_name += '_' + dst['suffix']

        # create list of available bands (to prevent duplicate names)
        band_names = [self.dataset.GetRasterBand(i + 1).GetMetadataItem(str('name'))
                        for i in range(self.dataset.RasterCount)]

        # check if name already exist and add '_NNNN'
        dst_band_name = band_name
        for n in range(9999):
            if dst_band_name not in band_names:
                break
            dst_band_name = '%s_%04d' % (band_name, n)

        return dst_band_name, wkv

    def _find_complex_band(self):
        """Find complex data bands"""
        # find complex bands
        for i in range(1, self.dataset.RasterCount+1):
            if self.dataset.GetRasterBand(i).DataType in [8, 9, 10, 11]:
                return i
        return None

    def leave_few_bands(self, bands=None):
        """Leave only given bands in VRT"""
        if bands is None:
            return

        rm_bands = []
        for i in range(1, self.dataset.RasterCount+1):
            band_name = self.dataset.GetRasterBand(i).GetMetadata().get('name', '')
            if i not in bands and band_name not in bands:
                rm_bands.append(i)
        # delete bands from VRT
        self.delete_bands(rm_bands)

    def split_complex_bands(self):
        """Recursevly find complex bands and relace by real and imag components"""
        rm_metadata= ['dataType', 'PixelFunctionType']

        i = self._find_complex_band()
        if i is None:
            return

        band = self.dataset.GetRasterBand(i)
        band_array = band.ReadAsArray()
        band_metadata_orig = remove_keys(band.GetMetadata(), rm_metadata)
        band_name_orig = band_metadata_orig.get('name', 'complex_%003d'%i)
        # Copy metadata, modify 'name' and create VRTs for real and imag parts of each band
        band_description = []
        for part in ['real', 'imag']:
            band_metadata = band_metadata_orig.copy()
            band_name = band_name_orig + '_' + part
            band_metadata['name'] = band_name
            self.band_vrts[band_name] = VRT.from_array(eval('band_array.' + part))
            band_description.append({'src': {
                     'SourceFilename': self.band_vrts[band_name].filename,
                     'SourceBand':  1},
                     'dst': band_metadata})
        # create bands with parts
        self.create_bands(band_description)
        # delete the complex band
        self.delete_band(i)
        # find and split more complex bands
        self.split_complex_bands()

    def create_geolocation_bands(self):
        """Create bands from Geolocation"""
        if self.geolocation is not None and len(self.geolocation.data) > 0:
            self.create_band(
                {'SourceFilename': self.geolocation.data['X_DATASET'],
                 'SourceBand': int(self.geolocation.data['X_BAND'])},
                {'wkv': 'longitude',
                 'name': 'longitude'})
            self.create_band(
                {'SourceFilename': self.geolocation.data['Y_DATASET'],
                 'SourceBand': int(self.geolocation.data['Y_BAND'])},
                {'wkv': 'latitude',
                 'name': 'latitude'})
        self.dataset.FlushCache()

    def fix_band_metadata(self, rm_metadata):
        """Add NETCDF_VARNAME and remove <rm_metadata> in metadata for each band"""
        for iBand in range(self.dataset.RasterCount):
            band = self.dataset.GetRasterBand(iBand + 1)
            metadata = remove_keys(band.GetMetadata(), rm_metadata)
            if 'name' in metadata:
                metadata['NETCDF_VARNAME'] = metadata['name']
            band.SetMetadata(metadata)
        self.dataset.FlushCache()

    def fix_global_metadata(self, rm_metadata):
        """Remove unwanted global metadata and escape special characters"""
        metadata = remove_keys(self.dataset.GetMetadata(), rm_metadata)
        # Apply escaping to metadata strings to preserve special characters (in XML/HTML format)
        metadata_escaped = {}
        for key, val in list(metadata.items()):
            # Keys not escaped - this may be changed if needed...
            metadata_escaped[key] = gdal.EscapeString(val, gdal.CPLES_XML)
        self.dataset.SetMetadata(metadata_escaped)
        self.dataset.FlushCache()

    def hardcopy_bands(self):
        """Make 'hardcopy' of bands: evaluate array from band and put into original band"""
        bands = range(1, self.dataset.RasterCount+1)
        for i in bands:
            self.band_vrts[i] = VRT.from_array(self.dataset.GetRasterBand(i).ReadAsArray())

        node0 = Node.create(str(self.xml))
        for i, iNode1 in enumerate(node0.nodeList('VRTRasterBand')):
            iNode1.node('SourceFilename').value = self.band_vrts[i+1].filename
            iNode1.node('SourceBand').value = str(1)
        self.write_xml(node0.rawxml())

    def prepare_export_gtiff(self):
        """Prepare dataset for export using GTiff driver"""
        if len(self.dataset.GetGCPs()) > 0:
            self._remove_geotransform()
        return False

    def prepare_export_netcdf(self):
        """Prepare dataset for export using netCDF driver"""
        gcps = self.dataset.GetGCPs()
        srs = str(self.get_projection()[0])
        if len(gcps) > 0:
            self._remove_geotransform()
            self.dataset.SetMetadataItem(str('NANSAT_GCPProjection'),
                                         str(srs.replace(',', '|').replace('"', '&')))
            add_gcps = True
        else:
            self.dataset.SetMetadataItem(str('NANSAT_Projection'),
                                         str(srs.replace(',', '|').replace('"', '&')))

            # add GeoTransform metadata
            geo_transform_str = str(self.dataset.GetGeoTransform()).replace(',', '|')
            self.dataset.SetMetadataItem(str('NANSAT_GeoTransform'), str(geo_transform_str))
            add_gcps = False
        return add_gcps

    def copy(self):
        """Create and return a full copy of a VRT instance with new filenames

        If self.dataset has no bands, the copy is created also without bands.
        If self.dataset has bands, the copy is created from the dataset with all bands.
        If self has attribute 'vrt' (a sub-VRT object, result of get_super_vrt) it's contents is
        also copied into sub-VRT of the copy.
        Other attributes of self, such as tps flag and band_vrts are also copied.

        """
        if self.dataset.RasterCount == 0:
            new_vrt = VRT.from_gdal_dataset(self.dataset, geolocation=self.geolocation,
                                                          metadata=self.dataset.GetMetadata())
        else:
            new_vrt = VRT.copy_dataset(self.dataset, geolocation=self.geolocation,
                                                     metadata=self.dataset.GetMetadata())
            replace_filenames = [(self.filename, new_vrt.filename)]
            if self.vrt is not None:
                # recursive copy of sub-VRT object
                new_vrt.vrt = self.vrt.copy()
                replace_filenames.append((self.vrt.filename, new_vrt.vrt.filename))

            # change reference from original filenames to the new ones
            for replace_filename in replace_filenames:
                new_vrt_xml = new_vrt.xml.replace(os.path.basename(replace_filename[0]),
                                          os.path.basename(replace_filename[1]))
                new_vrt.write_xml(new_vrt_xml)

            # copy VRTs of bands
            new_vrt.band_vrts = dict(self.band_vrts)
        # copy the thin spline transformation option
        new_vrt.tps = bool(self.tps)

        return new_vrt

    @property
    def xml(self):
        """Read XML content of the VRT-file using VSI

        Returns
        --------
        string : XMl Content which is read from the VSI file
        """
        self.dataset.FlushCache()
        return VRT.read_vsi(self.filename)

    def create_bands(self, metadata_dict):
        """ Generic function called from the mappers to create bands
        in the VRT dataset from an input dictionary of metadata

        Parameters
        ----------
        metadata_dict : list of dict with params of input bands and generated bands.
            Each dict has:
            'src' : dictionary with parameters of the sources:
            'dst' : dictionary with parameters of the generated bands

        Notes
        ---------
        Adds bands to the self.dataset based on info in metaDict

        See Also
        ---------
        VRT.create_band()

        """
        for band_dict in metadata_dict:
            src = band_dict['src']
            dst = band_dict.get('dst', None)
            self.create_band(src, dst)
            self.logger.debug('Creating band - OK!')
        self.dataset.FlushCache()

    def create_band(self, src, dst=None):
        """ Add band to self.dataset:

        Get parameters of the source band(s) from input
        Generate source XML for the VRT, add options of creating
        Call GDALDataset.AddBand
        Set source and options
        Add metadata

        Parameters
        ----------
        src : dict with parameters of sources
            SourceFilename,
            SourceBand,
            ScaleRatio,
            ScaleOffset,
            NODATA,
            LUT,
            SourceType,
            DataType,
            ImageOffset (RawVRT),
            PixelOffset (RawVRT),
            LineOffset (RawVRT),
            ByteOrder (RawVRT),
            xSize,
            ySize
        dst : dict with parameters of the created band
            name,
            dataType,
            wkv,
            suffix,
            AnyOtherMetadata,
            PixelFunctionType (
            1) band will be a pixel function defined by the
            corresponding name/value.
            In this case src may be list of
            dicts with parameters for each source.
            2) in case the dst band has a different datatype
            than the source band it is important to add a
            SourceTransferType parameter in dst),
            SourceTransferType

        Returns
        --------
        name : string, name of the added band

        Examples
        --------
        >>> vrt.create_band({'SourceFilename': filename, 'SourceBand': 1})

        >>> vrt.create_band({'SourceFilename': filename, 'SourceBand': 2,
                             'ScaleRatio': 0.0001},
                            {'name': 'LAT', 'wkv': 'latitude'})

        >>> vrt.create_band({'SourceFilename': filename, 'SourceBand': 2},
                            {'suffix': '670',
                             'wkv': 'brightness_temperature'})

        >>> vrt.create_band([{'SourceFilename': filename, 'SourceBand': 1},
                             {'SourceFilename': filename, 'SourceBand': 1}],
                             {'PixelFunctionType': 'NameOfPixelFunction'})

        """
        self.logger.debug('INPUTS: %s, %s " ' % (str(src), str(dst)))
        # Make sure src is list, ready for loop
        if type(src) == dict:
            srcs = [src]
        elif type(src) in [list, tuple]:
            srcs = src

        # Check if dst is given, or create empty dict
        if dst is None:
            dst = {}

        srcs = list(map(VRT._make_source_bands_xml, srcs))
        options = VRT._set_add_band_options(srcs, dst)
        dst['dataType'] = VRT._get_dst_band_data_type(srcs, dst)
        dst['name'], wkv = self._create_band_name(dst)

        # Add Band
        self.dataset.AddBand(int(dst['dataType']), options=options)
        dst_raster_band = self.dataset.GetRasterBand(self.dataset.RasterCount)

        # Append sources to destination dataset
        if len(srcs) == 1 and srcs[0]['SourceBand'] > 0:
            # only one source
            dst_raster_band.SetMetadataItem(str('source_0'), str(srcs[0]['XML']),
                                            str('new_vrt_sources'))
        elif len(srcs) > 1:
            # several sources for PixelFunction
            metadataSRC = {}
            for i, src in enumerate(srcs):
                metadataSRC['source_%d' % i] = src['XML']
            dst_raster_band.SetMetadata(metadataSRC, str('vrt_sources'))

        # set metadata from WKV
        dst_raster_band = VRT._put_metadata(dst_raster_band, wkv)

        # set metadata from provided parameters
        # remove and add params
        dst['SourceFilename'] = srcs[0]['SourceFilename']
        dst['SourceBand'] = str(srcs[0]['SourceBand'])
        dst_raster_band = VRT._put_metadata(dst_raster_band, dst)

        # return name of the created band
        return dst['name']

    def write_xml(self, vsi_file_content=None):
        """Write XML content into a VRT dataset

        Parameters
        -----------
        vsi_fileContent: string, optional
            XML Content of the VSI file to write

        Notes
        -----
        self.dataset
            If XML content was written, self.dataset is re-opened

        """
        vsi_file = gdal.VSIFOpenL(self.filename, str('w'))
        gdal.VSIFWriteL(str(vsi_file_content), len(vsi_file_content), 1, vsi_file)
        gdal.VSIFCloseL(vsi_file)
        # re-open self.dataset with new content
        self.dataset = gdal.Open(self.filename)

    def export(self, filename):
        """Export VRT file as XML into given <filename>"""
        self.driver.CreateCopy(filename, self.dataset)

    def _get_sub_filenames(self, gdal_dataset):
        """ Get filenames of subdatasets

        Parameters
        ----------
        gdal_dataset : gdal.Dataset
            Main gdal dataset containing sub-datasets

        Returns
        -------
        filenames : list
           List of filenames of each sub-dataset within the netCDF file

        """
        return [f[0] for f in gdal_dataset.GetSubDatasets()]

    def get_warped_vrt(self, dst_srs, x_size, y_size, geo_transform,
                       resample_alg=0,
                       dst_gcps=[],
                       skip_gcps=1,
                       block_size=None,
                       working_data_type=None,
                       resize_only=False):

        """Create warped (reprojected) VRT object

        1. Create a simple warped dataset using GDAL.AutoCreateWarpedVRT.
        2. Create VRT from the warped dataset.
        3. Modify the warped VRT according to the input options (size, geotransform, GCPs, etc)
        4. Keep the original VRT in the attribute vrt

        For reprojecting the function tries to use geolocation by default,
        if geolocation is not present it tries to use GCPs,
        if GCPs are not present it uses GeoTransform.

        If destination image has GCPs (provided in <dst_gcps>): fake GCPs for
        referencing line/pixel of SRC image and X/Y of DST image are created
        and added to the SRC image. After warping the dst_gcps are added to
        the warped VRT.

        Parameters
        -----------
        dst_srs : string
            WKT of the destination projection
        x_size, y_size : int
            dimentions of the destination rasetr
        geo_transform : tuple with 6 floats
            destination GDALGeoTransfrom
        resample_alg : int (GDALResampleAlg)
            0 : NearestNeighbour,
            1 : Bilinear,
            2 : Cubic,
            3 : CubicSpline,
            4 : Lancoz
        dst_gcps : list with GDAL GCPs
            GCPs of the destination image
        skip_gcps : int
            Using TPS can be very slow if the number of GCPs are large.
            If this parameter is given, only every [skip_gcp] GCP is used,
            improving calculation time at the cost of accuracy.
        block_size : int
            BlockSize to use for resampling. Larger blocksize reduces speed but increases accuracy.
        working_data_type : str
            'Float32', 'Int16', etc.
        resize_only : bool
            Create warped_vrt which will be used for resizing only?

        Returns
        --------
        warped_vrt : VRT object with warped dataset and with vrt


        """
        # VRT to be warped
        src_vrt = self.copy()

        # if destination GCPs are given: create and add fake GCPs to src
        dst_wkt = src_vrt._set_fake_gcps(dst_srs, dst_gcps, skip_gcps)

        if resize_only:
            src_vrt._set_geotransform_for_resize()
        else:
            src_vrt._set_gcps_geolocation_geotransform()

        # create Warped VRT GDAL Dataset
        warped_dataset = gdal.AutoCreateWarpedVRT(src_vrt.dataset, None, dst_wkt, resample_alg)

        # check if Warped VRT was created
        if warped_dataset is None:
            raise NansatGDALError('Cannot create warpedVRT!')

        # create VRT object from Warped VRT GDAL Dataset
        warped_vrt = VRT.copy_dataset(warped_dataset)

        # set x/y size, geo_transform, block_size
        warped_vrt._update_warped_vrt_xml(x_size, y_size, geo_transform, block_size, working_data_type)

        # apply thin-spline-transformation option
        if self.tps:
            warped_vrt.write_xml(warped_vrt.xml.replace('GCPTransformer', 'TPSTransformer'))

        # if given, add dst GCPs
        if len(dst_gcps) > 0:
            warped_vrt.dataset.SetGCPs(dst_gcps, dst_srs)
            warped_vrt._remove_geotransform()
            warped_vrt.dataset.SetProjection(str(''))

        # Copy self to warpedVRT
        warped_vrt.vrt = self.copy()

        # replace the reference from src_vrt to warped_vrt.vrt
        node0 = Node.create(str(warped_vrt.xml))
        node1 = node0.node('GDALWarpOptions')
        node1.node('SourceDataset').value = '/vsimem/' + str(os.path.basename(warped_vrt.vrt.filename))
        warped_vrt.write_xml(node0.rawxml())

        return warped_vrt

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

        Create shift_vrt which references self. Modify georeference
        of shift_vrt to account for the roll. Add as many bands as in self
        but for each band create two complex sources: for western
        and eastern parts. Keep self in shift_vrt.vrt

        Parameters
        ----------
        shift_degree : float
            rolling angle, how far east/west to roll

        Returns
        -------
        shift_vrt : VRT object with rolled bands

        """
        # Copy self into self.vrt
        shift_vrt = self.get_super_vrt()

        # cancel shift if not longlat projection
        projection = self.get_projection()
        if projection[1] != 'dataset' or 'AUTHORITY["EPSG","4326"]]' not in projection[0]:
            return shift_vrt

        if shift_degree < 0:
            shift_degree += 360.0

        geo_transform = shift_vrt.vrt.dataset.GetGeoTransform()
        shift_pixel = int(shift_degree / float(geo_transform[1]))
        geo_transform = list(geo_transform)
        geo_transform[0] = round(geo_transform[0] + shift_degree, 3)
        new_east_border = geo_transform[0] + (geo_transform[1] * shift_vrt.dataset.RasterXSize)
        if new_east_border > 360.0:
            geo_transform[0] -= 360.0
        shift_vrt.dataset.SetGeoTransform(tuple(geo_transform))

        # read xml and create the node
        node0 = Node.create(shift_vrt.xml)

        # divide into two bands and switch the bands
        for i in range(len(node0.nodeList('VRTRasterBand'))):
            # create i-th 'VRTRasterBand' node
            node1 = node0.node('VRTRasterBand', i)
            # modify the 1st band
            shift_str = str(shift_pixel)
            size_str = str(shift_vrt.vrt.dataset.RasterXSize - shift_pixel)
            node1.node('ComplexSource').node('DstRect').replaceAttribute('xOff', shift_str)
            node1.node('ComplexSource').node('DstRect').replaceAttribute('xSize', size_str)
            node1.node('ComplexSource').node('SrcRect').replaceAttribute('xSize', size_str)

            # add the 2nd band
            clone_node = Node.create(node1.rawxml()).node('ComplexSource')
            clone_node.node('SrcRect').replaceAttribute('xOff', size_str)
            clone_node.node('DstRect').replaceAttribute('xOff', str(0))
            clone_node.node('SrcRect').replaceAttribute('xSize', shift_str)
            clone_node.node('DstRect').replaceAttribute('xSize', shift_str)

            # get VRTRasterBand with inserted ComplexSource
            node1 = node1.insert(clone_node.rawxml())
            node0.replaceNode('VRTRasterBand', i, node1)

        # write down XML contents
        shift_vrt.write_xml(node0.rawxml())

        return shift_vrt

    def get_sub_vrt(self, steps=1):
        """Returns sub-VRT from given depth

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

        Notes
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
        """Create a new VRT object with a reference to the current object (self)

        Create a new VRT (super_vrt) with exactly the same structure (number of bands, raster size,
        metadata) as the current object (self). Create a copy of the current object and add it as
        an attribute of the new object (super_vrt.vrt). Bands in the new object will refer to the
        same bands in the current object. Recursively copy all vrt attributes of the current
        object (self.vrt.vrt.vrt...) into the new object (super_vrt.vrt.vrt.vrt.vrt...).


        Returns
        -------
        super_vrt : VRT
            a new VRT object with copy of self in super_vrt.vrt

        """
        # create new vrt that refers to a copy of self
        super_vrt = VRT.from_gdal_dataset(self.dataset, geolocation=self.geolocation,
                                                        metadata=self.dataset.GetMetadata())
        super_vrt.vrt = self.copy()
        super_vrt.tps = self.tps

        # add bands to the new vrt
        for i in range(super_vrt.vrt.dataset.RasterCount):
            src = {'SourceFilename': super_vrt.vrt.filename, 'SourceBand': i + 1}
            dst = super_vrt.vrt.dataset.GetRasterBand(i + 1).GetMetadata()
            # remove PixelFunctionType from metadata to prevent its application
            if 'PixelFunctionType' in dst:
                dst.pop('PixelFunctionType')
            super_vrt.create_band(src, dst)
        super_vrt.dataset.FlushCache()

        return super_vrt

    def get_subsampled_vrt(self, new_raster_x_size, new_raster_y_size, resample_alg):
        """Create VRT and replace step in the source"""

        subsamp_vrt = self.get_super_vrt()

        # Get XML content from VRT-file
        node0 = Node.create(str(subsamp_vrt.xml))

        # replace rasterXSize in <VRTDataset>
        node0.replaceAttribute('rasterXSize', str(new_raster_x_size))
        node0.replaceAttribute('rasterYSize', str(new_raster_y_size))

        # replace xSize in <DstRect> of each source
        for iNode1 in node0.nodeList('VRTRasterBand'):
            for sourceName in ['ComplexSource', 'SimpleSource']:
                for iNode2 in iNode1.nodeList(sourceName):
                    iNodeDstRect = iNode2.node('DstRect')
                    iNodeDstRect.replaceAttribute('xSize', str(new_raster_x_size))
                    iNodeDstRect.replaceAttribute('ySize', str(new_raster_y_size))
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
        subsamp_vrt.write_xml(node0.rawxml())

        return subsamp_vrt

    def transform_points(self, col_vector, row_vector, dst2src=0,
                         dst_srs=None, dst_ds=None, options=None):
        """Transform input pixel/line coordinates into lon/lat (or opposite)

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
        if dst_srs is None:
            dst_srs = NSR()
        # get source SRS (either Projection or GCPProjection or Metadata(GEOLOCATION)[SRS])
        src_wkt = self.get_projection()[0]

        # prepare options
        if options is None:
            options = ['SRC_SRS=' + src_wkt, 'DST_SRS=' + dst_srs.wkt]
            # add TPS method if we have GCPs and self.tps is True
            if self.tps and len(self.dataset.GetGCPs()) > 0:
                options.append('METHOD=GCP_TPS')

        # create transformer
        transformer = gdal.Transformer(self.dataset, dst_ds, options)

        # convert lists with X,Y coordinates to 2D numpy array
        xy = np.array([col_vector, row_vector]).transpose()

        # transfrom coordinates (TransformPoints returns list of (X, Y, Z) tuples)
        lonlat = transformer.TransformPoints(dst2src, xy)[0]

        # convert to Nx3 numpy array (keep second dimention to allow empty inputs)
        lonlat = np.array(lonlat)
        lonlat.shape = int(lonlat.size/3), 3

        # convert to vectors with lon,lat values
        lon_vector, lat_vector, _ = lonlat.T

        return lon_vector, lat_vector

    def get_projection(self):
        """Get projection (spatial reference system) of the dataset

        Uses dataset.GetProjection() or dataset.GetGCPProjection()
        or dataset.GetMetadata('GEOLOCATION')['SRS']

        Returns
        -------
        projection : projection WKT
        source : str ['gcps' or 'dataset' or 'geolocation']

        Raises
        -------
        NansatProjectionError : occurs when the projection is empty.

        """
        projection = self.dataset.GetProjection()
        if projection != '':
            return projection, str('dataset')

        projection = self.dataset.GetGCPProjection()
        if projection != '':
            return projection, str('gcps')

        projection = self.dataset.GetMetadata(str('GEOLOCATION')).get(str('SRS'), '')
        if projection != '':
            return projection, str('geolocation')

        raise NansatProjectionError


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
        warped_vrt = self.get_warped_vrt(None, x_size, y_size, geo_transform, resample_alg,
                                         resize_only=True)

        return warped_vrt

    def reproject_gcps(self, dst_srs):
        """Reproject all GCPs to a new spatial reference system

        Necessary before warping an image if the given GCPs
        are in a coordinate system which has a singularity
        in (or near) the destination area (e.g. poles for lonlat GCPs)

        Parameters
        ----------
        dst_srs : proj4, WKT, NSR, EPSG
            Destiination SRS given as any NSR input parameter

        Notes
        -----
        Reprojects all GCPs to new SRS and updates GCPProjection

        """
        # transform coordinates of original GCPs
        dst_srs = NSR(dst_srs)
        src_srs = NSR(self.dataset.GetGCPProjection())
        # make three tuples with X, Y and Z values of GCPs
        src_points = list(zip(*[(gcp.GCPX, gcp.GCPY, gcp.GCPZ) for gcp in self.dataset.GetGCPs()]))
        dst_points = VRT.transform_coordinates(src_srs, src_points, dst_srs)
        # create new GCPs
        dst_gcps = []
        for p in zip(self.dataset.GetGCPs(), *dst_points):
            dst_gcp = gdal.GCP(p[1], p[2], p[3],
                               p[0].GCPPixel,  p[0].GCPLine, p[0].Info, p[0].Id)
            dst_gcps.append(dst_gcp)
        # Update dataset
        self.dataset.SetGCPs(dst_gcps, dst_srs.wkt)
        self.dataset.FlushCache()

    @staticmethod
    def transform_coordinates(src_srs, src_points, dst_srs):
        """ Transform coordinates of points from one spatial reference to another

        Parameters
        ----------
        src_srs : nansat.NSR
            Source spatial reference system
        src_points : tuple of two or three N-D arrays
            Coordinates of points in the source spatial reference system. A tuple with (X, Y) or
            (X, Y, Z) coordinates arrays. Each coordinate Each array can be a list, 1D, 2D, N-D
            array.
        dst_srs : nansat.NSR
            Destination spatial reference

        Returns
        -------
        dst_points : tuple of two or three N-D arrays
            Coordinates of points in the destination spatial reference system. A tuple with (X, Y) or
            (X, Y, Z) coordinates arrays. Each coordinate Each array can be 1D, 2D, N-D
            array. Shape of output arrays corrrrespond to shape of inputs.

        """
        transformer = osr.CoordinateTransformation(src_srs, dst_srs)
        src_shape = np.array(src_points[0]).shape
        src_points = list(zip(*[np.array(xyz).flatten() for xyz in src_points]))
        dst_points = transformer.TransformPoints(src_points)
        dst_x, dst_y, dst_z = np.array(list(zip(*dst_points)))
        return dst_x.reshape(src_shape), dst_y.reshape(src_shape), dst_z.reshape(src_shape)

    def set_offset_size(self, axis, offset, size):
        """Set offset and  size in VRT dataset and band attributes

        Parameters
        ----------
        axis : str
            name of axis ('x' or 'y')
        offset : int
            value of offset to put into VRT
        size : int
            value of size to put into VRT

        Notes
        -----
        Changes VRT file, sets new offset and size

        """
        node0 = Node.create(str(self.xml))

        # change size
        node0.node('VRTDataset').replaceAttribute('raster%sSize'%str(axis).upper(), str(size))

        # replace x/y-Off and x/y-Size
        #   in <SrcRect> and <DstRect> of each source
        for iNode1 in node0.nodeList('VRTRasterBand'):
            iNode2 = iNode1.node('ComplexSource')

            iNode3 = iNode2.node('SrcRect')
            iNode3.replaceAttribute('%sOff'%str(axis).lower(), str(offset))
            iNode3.replaceAttribute('%sSize'%str(axis).lower(), str(size))

            iNode3 = iNode2.node('DstRect')
            iNode3.replaceAttribute('%sSize'%str(axis).lower(), str(size))

        # write modified XML
        self.write_xml(node0.rawxml())

    def shift_cropped_gcps(self, x_offset, x_size, y_offset, y_size):
        """Modify GCPs to fit the size/offset of cropped image"""
        gcps = self.dataset.GetGCPs()
        if len(gcps) == 0:
            return

        dst_gcps = []
        i = 0
        # keep current GCPs
        for igcp in gcps:
            if (0 < igcp.GCPPixel - x_offset < x_size and 0 < igcp.GCPLine - y_offset < y_size):
                i += 1
                dst_gcps.append(gdal.GCP(igcp.GCPX, igcp.GCPY, 0,
                                         igcp.GCPPixel - x_offset,
                                         igcp.GCPLine - y_offset, str(''), str(i)))
        n_gcps = i
        if n_gcps < 100:
            # create new 100 GPCs (10 x 10 regular matrix)
            pix_array, lin_array = np.mgrid[0:x_size:10j, 0:y_size:10j]
            pix_array = pix_array.flatten() + x_offset
            lin_array = lin_array.flatten() + y_offset

            lon_array, lat_array = self.vrt.transform_points(pix_array, lin_array,
                                                dst_srs=NSR(self.vrt.dataset.GetGCPProjection()))

            for i in range(len(lon_array)):
                dst_gcps.append(gdal.GCP(lon_array[i], lat_array[i], 0,
                                         pix_array[i] - x_offset,
                                         lin_array[i] - y_offset,
                                         str(''), str(n_gcps+i+1)))

        # set new GCPss
        self.dataset.SetGCPs(dst_gcps, self.vrt.dataset.GetGCPProjection())
        # remove geotranform which was automatically added
        self._remove_geotransform()

    def shift_cropped_geo_transform(self, x_offset, x_size, y_offset, y_size):
        """Modify GeoTransform to fit the size/offset of the cropped image"""
        if len(self.vrt.dataset.GetGCPs()) != 0:
            return
        # shift upper left corner coordinates
        geoTransfrom = self.dataset.GetGeoTransform()
        geoTransfrom = list(map(float, geoTransfrom))
        geoTransfrom[0] += geoTransfrom[1] * x_offset
        geoTransfrom[3] += geoTransfrom[5] * y_offset
        self.dataset.SetGeoTransform(geoTransfrom)
        self.dataset.FlushCache()

    @staticmethod
    def read_vsi(filename):
        """Read text from input <filename:str> using VSI and return <content:str>."""
        # open
        vsi_file = gdal.VSIFOpenL(str(filename), str('r'))
        # get file size
        gdal.VSIFSeekL(vsi_file, 0, 2)
        vsi_file_size = gdal.VSIFTellL(vsi_file)
        # fseek to start again
        gdal.VSIFSeekL(vsi_file, 0, 0)
        # read
        vsi_file_content = gdal.VSIFReadL(vsi_file_size, 1, vsi_file)
        gdal.VSIFCloseL(vsi_file)
        return str(vsi_file_content.decode())

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
            raise KeyError('SourceFilename must be available in dict:src_in')
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

        take <n_gcps> regular pixels from inpt <lat> and <lon> grids
        Create GCPs from these pixels
        Create latlong GCPs projection

        Parameters
        -----------
        lat : Numpy grid
            array of latitudes
        lon : Numpy grid
            array of longitudes (should be the same size as lat)
        n_gcps : int, optional, default = 100
            number of GCPs to create

        Returns
        --------
        gcsp : List with GDAL GCPs

        """
        # estimate step of GCPs
        gcp_size = np.sqrt(n_gcps)
        step0 = max(1, int(float(lat.shape[0]) / gcp_size))
        step1 = max(1, int(float(lat.shape[1]) / gcp_size))

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
            all_chars = ascii_uppercase + digits
            random_chars = ''.join(choice(all_chars) for x in range(10))
            filename = '/vsimem/%s.%s' % (random_chars, extention)
        return filename

    @staticmethod
    def _get_dst_band_data_type(srcs, dst):
        """Get destination dataType based on source and destnation band dictionaries"""
        if 'dataType' in dst:                                               # data type is preset
            data_type = dst['dataType']
        elif (len(srcs) > 1 or                                                     # pixel function
              ('ScaleRatio' in srcs[0] and float(srcs[0]['ScaleRatio']) != 1.0) or # scaling applied
              ('LUT' in srcs[0] and len(srcs[0]['LUT'])) > 0 or                    # if LUT exists
              'DataType' not in srcs[0]):               # if source band not available
            data_type = gdal.GDT_Float32
        else:
            # otherwise take the DataType from source
            data_type = srcs[0]['DataType']
        return data_type

    @staticmethod
    def _remove_strings_in_metadata_keys(gdal_metadata, rm_strings):
        """Remove unwanted metadata"""

        new_meta = {}
        for key, value in gdal_metadata.items():
            for rms in rm_strings:
                if rms in key:
                    key = key.replace(rms, '')
            new_meta[key] = value

        return new_meta
