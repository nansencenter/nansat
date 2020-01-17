# Name:  nansat.py
# Purpose: Container of Nansat class
# Authors:      Anton Korosov, Knut-Frode Dagestad, Morten W. Hansen, Artem Moiseev,
#               Asuka Yamakawa, Alexander Myasoyedov,
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
from __future__ import absolute_import, print_function
import os
import glob
import sys
import tempfile
import datetime
import pkgutil
import warnings
from xml.sax import saxutils

import numpy as np
from numpy import nanmedian
from numpy.lib.recfunctions import append_fields
from netCDF4 import Dataset

try:
    from matplotlib import cm
except ImportError:
    MATPLOTLIB_IS_INSTALLED = False
else:
    MATPLOTLIB_IS_INSTALLED = True

from nansat.domain import Domain
from nansat.exporter import Exporter
from nansat.figure import Figure
from nansat.vrt import VRT
from nansat.utils import add_logger, gdal, parse_time
from nansat.node import Node
from nansat.pointbrowser import PointBrowser

from nansat.exceptions import NansatGDALError, WrongMapperError, NansatReadError

import collections
if hasattr(collections, 'OrderedDict'):
    from collections import OrderedDict
else:
    from ordereddict import OrderedDict

# container for all mappers
nansatMappers = None


class Nansat(Domain, Exporter):
    """Container for geospatial data. Performs all high-level operations.

    n = Nansat(filename) opens the file with satellite or model data for
    reading, adds scientific metadata to bands, and prepares the data for
    further handling.

    Parameters
    -----------
    filename : str
        uri of the input file or OpeNDAP datastream
    mapper : str
        name of the mapper from nansat/mappers dir. E.g.
        'sentinel1_l1', 'asar', 'hirlam', 'meris_l1', 'meris_l2', etc.
    log_level : int
        Level of logging. See: http://docs.python.org/howto/logging.html
    kwargs : additional arguments for mappers

    Examples
    --------
        >>> n1 = Nansat(filename)
        >>> n2 = Nansat(sentinel1_filename, mapper='sentinel1_l1')
        >>> array1 = n1[1]
        >>> array2 = n2['sigma0_HV']

    Notes
    -----
    The instance of Nansat class (the object <n>) contains information
    about geographical reference of the data (e.g raster size, pixel
    resolution, type of projection, etc) and about bands with values of
    geophysical variables (e.g. water leaving radiance, normalized radar
    cross section, chlrophyll concentraion, etc). The object <n> has methods
    for high-level operations with data. E.g.:
    * reading data from file (Nansat.__getitem__);
    * visualization (Nansat.write_figure);
    * changing geographical reference (Nansat.reproject);
    * exporting (Nansat.export)
    * and much more...

    Nansat inherits from Domain (container of geo-reference information)
    Nansat uses instance of VRT (wraper around GDAL VRT-files)
    Nansat uses instance of Figure (collection of methods for visualization)

    """

    FILL_VALUE = 9.96921e+36
    ALT_FILL_VALUE = -10000.

    # instance attributes
    logger = None
    filename = None
    name = None
    path = None
    vrt = None
    mapper = None

    @classmethod
    def from_domain(cls, domain, array=None, parameters=None, log_level=30):
        """Create Nansat object from input Domain [and array with data]

        Parameters
        ----------
        domain : Domain
            Defines spatial reference system and geographical extent.
        array : numpy NDarray
            Data for the first band. Shape must correspond to shape of <domain>
        parameters : dict
            Metadata for the first band. May contain 'name', 'wkv' and other keys.
        log_level : int
            Level of logging.

        """
        n = cls.__new__(cls)
        n._init_from_domain(domain, array, parameters, log_level)
        return n

    def __init__(self, filename='', mapper='', log_level=30, **kwargs):
        """Create Nansat object

        Notes
        -----
        self.mapper : str
            name of the used mapper
        self.filename : file name
            set file name given by the argument
        self.vrt : VRT
            Wrapper around VRT file and GDAL dataset with satellite raster data
        self.logger : logging.Logger
            logger for output debugging info
        self.name : str
            name of object (for writing KML)
        self.path : str
            path to input file

        """
        if filename == '':
            raise ValueError('Nansat is called without valid parameters! Use: Nansat(filename)')

        self._init_empty(filename, log_level)
        # Create VRT object with mapping of variables
        self.vrt = self._get_mapper(mapper, **kwargs)

    def __getitem__(self, band_id):
        """Returns the band as a NumPy array, by overloading []

        Parameters
        -----------
        band_id : int or str
            If int, array from band with number <band_id> is returned
            If string, array from band with metadata 'name' equal to
            <band_id> is returned

        Returns
        --------
        a : NumPy array


        """
        # get band
        band = self.get_GDALRasterBand(band_id)
        # get expression from metadata
        expression = band.GetMetadata().get('expression', '')
        # get data
        band_data = band.ReadAsArray()
        if band_data is None:
            raise NansatGDALError('Cannot read array from band %s' % str(band_data))

        # execute expression if any
        if expression != '':
            band_data = eval(expression)

        all_float_flag = band_data.dtype.char in np.typecodes['AllFloat']
        # Set invalid and missing data to np.nan (for floats only)
        if '_FillValue' in band.GetMetadata() and all_float_flag:
            band_data = self._fill_with_nan(band, band_data)

        # replace infs with np.NAN
        if np.size(np.where(np.isinf(band_data))) > 0:
            band_data[np.isinf(band_data)] = np.nan

        # erase out-of-swath pixels with np.Nan (if not integer)
        if self.has_band('swathmask') and all_float_flag:
            swathmask = self.get_GDALRasterBand('swathmask').ReadAsArray()
            band_data[swathmask == 0] = np.nan

        return band_data

    def __repr__(self):
        """Creates string with basic info about the Nansat object"""
        out_str = '{separator}{filename}{separator}Mapper: {mapper}{bands}{separator}{domain}'
        return out_str.format(separator=self.OUTPUT_SEPARATOR, filename=self.filename,
                              bands=self.list_bands(False), mapper=self.mapper,
                              domain=Domain.__repr__(self))

    def _init_empty(self, filename, log_level):
        """Init empty Nansat object

        Parameters
        ----------
        filename : str
            Name of input file.
        log_level : int
            Level of logging verbosity.

        Notes
        --------
            self.logger
                adds logging.Logger
            self.filename
                adds full path to input file
            self.name
                adds name of file
            self.path
                adds path to file

        """
        # create logger
        self.logger = add_logger('Nansat', log_level)
        # set input file name
        self.filename = filename
        # name, for compatibility with some Domain methods
        self.name = os.path.basename(filename)
        self.path = os.path.dirname(filename)

    def _init_from_domain(self, domain, array=None, parameters=None, log_level=30):
        """Init Nansat object from input Domain and optionally array with band values

        Parameters
        ----------
        domain : Domain
            Defines spatial reference system and geographical extent.
        array : numpy NDarray
            Data for the first band. Shape must correspond to shape of <domain>
        parameters : dict
            Metadata for the first band. May contain 'name', 'wkv' and other keys.
        log_level : int
            Level of logging.

        """
        self._init_empty('', log_level)
        self.vrt = VRT.from_gdal_dataset(domain.vrt.dataset, geolocation=domain.vrt.geolocation)
        self.mapper = ''
        if array is not None:
            self.add_band(array=array, parameters=parameters)

    def _fill_with_nan(self, band, band_data):
        """Fill input array with fill value taen from input band metadata"""
        fill_value = float(band.GetMetadata()['_FillValue'])
        band_data[band_data == fill_value] = np.nan
        # quick hack to avoid problem with wrong _FillValue - see issue
        # #123
        if fill_value == self.FILL_VALUE:
            band_data[band_data == self.ALT_FILL_VALUE] = np.nan

        return band_data


    def add_band(self, array, parameters=None, nomem=False):
        """Add band from numpy array with metadata.

        Create VRT object which contains VRT and RAW binary file and append it
        to self.vrt.band_vrts

        Parameters
        -----------
        array : ndarray
            new band data. Shape should be equal to shape
        parameters : dict
            band metadata: wkv, name, etc. (or for several bands)
        nomem : bool
            saves the vrt to a tempfile on disk?

        Notes
        -----
        Creates VRT object with VRT-file and RAW-file. Adds band to the self.vrt.

        Examples
        --------
            >>> n.add_band(array, {'name': 'new_data'}) # add new band and metadata, keep in memory
            >>> n.add_band(array, nomem=True) # add new band, keep on disk

        """
        self.add_bands([array], [parameters], nomem)

    def add_bands(self, arrays, parameters=None, nomem=False):
        """Add bands from numpy arrays with metadata.

        Create VRT object which contains VRT and RAW binary file and append it
        to self.vrt.band_vrts

        Parameters
        -----------
        arrays : list of ndarrays
            new band data. Shape should be equal to shape
        parameters : list of dict
            band metadata: wkv, name, etc. (or for several bands)
        nomem : bool
            saves the vrt to a tempfile on disk?

        Notes
        -----
        Creates VRT object with VRT-file and RAW-file. Adds band to the self.vrt.

        Examples
        --------
            >>> n.add_bands([array1, array2]) # add new bands, keep in memory

        """
        # replace empty parameters with list of empty dictionaries
        if parameters is None:
            parameters = [{}] * len(arrays)

        self.vrt = self.vrt.get_super_vrt()

        # create VRTs from arrays and generate band_metadata
        band_metadata = []
        for array, parameter in zip(arrays, parameters):
            vrt = VRT.from_array(array, nomem=nomem)
            band_metadata.append({
                'src': {'SourceFilename': vrt.filename, 'SourceBand': 1},
                'dst': parameter
                })
            self.vrt.band_vrts[vrt.filename] = vrt

        band_name = self.vrt.create_bands(band_metadata)

    def bands(self):
        """Make a dictionary with all metadata from all bands

        Returns
        --------
        b : dictionary
            key = N, value = dict with all band metadata

        """
        band_metadata = {}
        for band_num in range(self.vrt.dataset.RasterCount):
            band_metadata[band_num + 1] = self.get_metadata(band_id=band_num + 1)

        return band_metadata

    def has_band(self, band):
        """Check if self has band with name <band>
        Parameters
        ----------
            band : str
                name or standard_name of the band to check

        Returns
        -------
            True/False if band exists or not

        """
        for band_num in self.bands():
            band_meta = self.bands()[band_num]
            if band_meta['name'] == band:
                return True
            elif 'standard_name' in band_meta and band_meta['standard_name'] == band:
                return True

    def _get_resize_shape(self, factor, width, height, dst_pixel_size):
        """Estimate new shape either from factor or destination width/height or pixel size"""
        src_shape = np.array(self.shape(), np.float)
        # estimate factor if either width or height is given and factor is not given
        if width is not None:
            factor = width / src_shape[1]
        if height is not None:
            factor = height / src_shape[0]

        # estimate factor if pixelsize is given
        if dst_pixel_size is not None:
            src_pixel_size = np.array(self.get_pixelsize_meters(), np.float)[::-1]
            factor = (src_pixel_size / float(dst_pixel_size)).mean()

        factor = float(factor)
        dst_shape = np.floor(src_shape * factor)
        self.logger.info('New shape: ({0}, {1})'.format(dst_shape[0], dst_shape[1]))
        return factor, dst_shape

    def resize(self, factor=None, width=None, height=None, pixelsize=None, resample_alg=-1):
        """Proportional resize of the dataset.

        The dataset is resized as (x_size*factor, y_size*factor)
        If desired width, height or pixelsize is specified,
        the scaling factor is calculated accordingly.
        If GCPs are given in a dataset, they are also rewritten.

        Parameters
        -----------
        factor : float, optional, default=1
            Scaling factor for width and height

            - > 1 means increasing domain size
            - < 1 means decreasing domain size

        width : int, optional
            Desired new width in pixels
        height : int, optional
            Desired new height in pixels
        pixelsize : float, optional
            Desired new pixelsize in meters (approximate).
            A factor is calculated from ratio of the
            current pixelsize to the desired pixelsize.
        resample_alg : int (GDALResampleAlg), optional

            - -1 : Average (default),
            - 0 : NearestNeighbour
            - 1 : Bilinear,
            - 2 : Cubic,
            - 3 : CubicSpline,
            - 4 : Lancoz

        Notes
        -----
        self.vrt.dataset : VRT dataset of VRT object
            raster size are modified to downscaled size.
            If GCPs are given in the dataset, they are also overwritten.

        """
        factor, dst_shape = self._get_resize_shape(factor, width, height, pixelsize)
        if resample_alg <= 0:
            self.vrt = self.vrt.get_subsampled_vrt(dst_shape[1], dst_shape[0], resample_alg)
        else:
            # update size and GeoTranform in XML of the warped VRT object
            self.vrt = self.vrt.get_resized_vrt(dst_shape[1], dst_shape[0], resample_alg)

        # resize gcps
        gcps = self.vrt.vrt.dataset.GetGCPs()
        if len(gcps) > 0:
            gcpPro = self.vrt.vrt.dataset.GetGCPProjection()
            for gcp in gcps:
                gcp.GCPPixel *= factor
                gcp.GCPLine *= factor
            self.vrt.dataset.SetGCPs(gcps, gcpPro)
            self.vrt._remove_geotransform()
        else:
            # change resultion in geotransform to keep spatial extent
            geoTransform = list(self.vrt.vrt.dataset.GetGeoTransform())
            geoTransform[1] = float(geoTransform[1])/factor
            geoTransform[5] = float(geoTransform[5])/factor
            geoTransform = list(map(float, geoTransform))
            self.vrt.dataset.SetGeoTransform(geoTransform)

        # set global metadata
        subMetaData = self.vrt.vrt.dataset.GetMetadata()
        subMetaData.pop('filename')
        self.set_metadata(subMetaData)

        return factor

    def get_GDALRasterBand(self, band_id=1):
        """Get a GDALRasterBand of a given Nansat object

        If str is given find corresponding band number
        If int is given check if band with this number exists.
        Get a GDALRasterBand from vrt.

        Parameters
        -----------
        band_id : int or str
            - if int - a band number of the band to fetch
            - if str band_id = {'name': band_id}

        Returns
        --------
        GDAL RasterBand

        Example
        -------
            >>> b = n.get_GDALRasterBand(1)
            >>> b = n.get_GDALRasterBand('sigma0')

        """
        # get band number
        bandNumber = self.get_band_number(band_id)
        # the GDAL RasterBand of the corresponding band is returned
        return self.vrt.dataset.GetRasterBand(bandNumber)

    def list_bands(self, do_print=True):
        """Show band information of the given Nansat object

        Show serial number, longName, name and all parameters
        for each band in the metadata of the given Nansat object.

        Parameters
        -----------
        do_print : boolean
            print on screen?

        Returns
        --------
        outString : String
            formatted string with bands info

        """
        # get dictionary of bands metadata
        bands = self.bands()
        outString = ''

        for b in bands:
            # print band number, name
            outString += 'Band : %d %s\n' % (b, bands[b].get('name', ''))
            # print band metadata
            for i in bands[b]:
                outString += '  %s: %s\n' % (i, bands[b][i])
        if do_print:
            # print to screeen
            print(outString)
        else:
            return outString

    def reproject(self, dst_domain=None, resample_alg=0,
                  block_size=None, tps=None, skip_gcps=1, addmask=True,
                  **kwargs):
        """
        Change projection of the object based on the given Domain

        Create superVRT from self.vrt with AutoCreateWarpedVRT() using
        projection from the dst_domain.
        Modify XML content of the warped vrt using the Domain parameters.
        Generate warpedVRT and replace self.vrt with warpedVRT.
        If current object spans from 0 to 360 and dst_domain is west of 0,
        the object is shifted by 180 westwards.

        Parameters
        -----------
        dst_domain : domain
            destination Domain where projection and resolution are set
        resample_alg : int (GDALResampleAlg)

            - 0 : NearestNeighbour
            - 1 : Bilinear
            - 2 : Cubic,
            - 3 : CubicSpline
            - 4 : Lancoz

        block_size : int
            size of blocks for resampling. Large value decrease speed
            but increase accuracy at the edge
        tps : bool
            Apply Thin Spline Transformation if source or destination has GCPs
            Usage of TPS can also be triggered by setting self.vrt.tps=True
            before calling to reproject.
            This options has priority over self.vrt.tps
        skip_gcps : int
            Using TPS can be very slow if the number of GCPs are large.
            If this parameter is given, only every [skip_gcp] GCP is used,
            improving calculation time at the cost of accuracy.
            If not given explicitly, 'skip_gcps' is fetched from the
            metadata of self, or from dst_domain (as set by mapper or user).
            [defaults to 1 if not specified, i.e. using all GCPs]
        addmask : bool
            If True, add band 'swathmask'. 1 - valid data, 0 no-data.
            This band is used to replace no-data values with np.nan

        Notes
        -----
        self.vrt : VRT object with dataset replaced to warpedVRT dataset
        Integer data is returnd by integer. Round off to decimal place.
        If you do not want to round off, convert the data types to
        GDT_Float32, GDT_Float64, or GDT_CFloat32.

        See Also
        ---------
        `<http://www.gdal.org/gdalwarp.html>`_

        """
        # This is time consuming and therefore not done...:
        #if not self.overlaps(dst_domain):
        #    raise ValueError('Source and destination domains do not overlap')

        # if self spans from 0 to 360 AND dst_domain is west of 0:
        #     shift self westwards by 180 degrees
        # check span
        srcCorners = self.get_corners()
        if round(min(srcCorners[0])) == 0 and round(max(srcCorners[0])) == 360:
            # check intersection of src and dst
            dstCorners = dst_domain.get_corners()
            if min(dstCorners[0]) < 0:
                # shift
                self.vrt = self.vrt.get_shifted_vrt(-180)

        # get projection of destination dataset
        dstSRS = dst_domain.vrt.dataset.GetProjection()

        # get destination GCPs
        dstGCPs = dst_domain.vrt.dataset.GetGCPs()
        if len(dstGCPs) > 0:
            # get projection of destination GCPs
            dstSRS = dst_domain.vrt.dataset.GetGCPProjection()

        x_size = dst_domain.vrt.dataset.RasterXSize
        y_size = dst_domain.vrt.dataset.RasterYSize

        geoTransform = dst_domain.vrt.dataset.GetGeoTransform()

        # set trigger for using TPS
        if tps is True:
            self.vrt.tps = True
        elif tps is False:
            self.vrt.tps = False

        # Reduce number of GCPs for faster reprojection
        # when using TPS (if requested)
        src_skip_gcps = self.vrt.dataset.GetMetadataItem('skip_gcps')
        dst_skip_gcps = dst_domain.vrt.dataset.GetMetadataItem('skip_gcps')
        kwargs['skip_gcps'] = skip_gcps  # default (use all GCPs)
        if dst_skip_gcps is not None:  # ...or use setting from dst
            kwargs['skip_gcps'] = int(dst_skip_gcps)
        if src_skip_gcps is not None:  # ...or use setting from src
            kwargs['skip_gcps'] = int(src_skip_gcps)

        # add band that masks valid values with 1 and nodata with 0
        # after reproject
        # TODD: REFACTOR: replace with VRT._add_swath_mask_band
        if addmask:
            self.vrt = self.vrt.get_super_vrt()
            src = [{
                'SourceFilename': self.vrt.vrt.filename,
                'SourceBand':  1,
                'DataType': gdal.GDT_Byte
            }]
            dst = {
                'dataType': gdal.GDT_Byte,
                'wkv': 'swath_binary_mask',
                'PixelFunctionType': 'OnesPixelFunc',
            }
            self.vrt.create_band(src=src, dst=dst)
            self.vrt.dataset.FlushCache()

        # create Warped VRT
        self.vrt = self.vrt.get_warped_vrt(dstSRS, x_size, y_size, geoTransform,
                                           resample_alg=resample_alg,
                                           dst_gcps=dstGCPs,
                                           block_size=block_size, **kwargs)

        # set global metadata from subVRT
        subMetaData = self.vrt.vrt.dataset.GetMetadata()
        subMetaData.pop('filename')
        self.set_metadata(subMetaData)

    def undo(self, steps=1):
        """Undo reproject, resize, add_band or crop of Nansat object

        Restore the self.vrt from self.vrt.vrt

        Parameters
        -----------
        steps : int
            How many steps back to undo

        Notes
        ------
        Modifies self.vrt

        """
        self.vrt = self.vrt.get_sub_vrt(steps)

    def watermask(self, mod44path=None, dst_domain=None, **kwargs):
        """
        Create numpy array with watermask (water=1, land=0)

        250 meters resolution watermask from MODIS 44W Product:
        http://www.glcf.umd.edu/data/watermask/

        Watermask is stored as tiles in TIF(LZW) format and a VRT file
        All files are stored in one directory.
        A tarball with compressed TIF and VRT files should be additionally
        downloaded from the Nansat documentation page:
        http://nansat.readthedocs.io/en/latest/source/features.html#differentiating-between-land-and-water

        The method :
            Gets the directory either from input parameter or from environment
            variable MOD44WPATH
            Open Nansat object from the VRT file
            Reprojects the watermask onto the current object using reproject()
            or reproject_on_jcps()
            Returns the reprojected Nansat object

        Parameters
        -----------
        mod44path : string
            path with MOD44W Products and a VRT file
        dst_domain : Domain
            destination domain other than self
        tps : Bool
            Use Thin Spline Transformation in reprojection of watermask?
            See also Nansat.reproject()
        skip_gcps : int
            Factor to reduce the number of GCPs by and increase speed
            See also Nansat.reproject()

        Returns
        --------
        watermask : Nansat object with water mask in current projection

        See Also
        ---------
        `<http://www.glcf.umd.edu/data/watermask/>`_

        `<http://nansat.readthedocs.io/en/latest/source/features.html#differentiating-between-land-and-water>`_

        """
        mod44DataExist = True
        # check if path is given in input param or in environment
        if mod44path is None:
            mod44path = os.getenv('MOD44WPATH')
        if mod44path is None:
            mod44DataExist = False
        # check if VRT file exist
        elif not os.path.exists(mod44path + '/MOD44W.vrt'):
            mod44DataExist = False
        self.logger.debug('MODPATH: %s' % mod44path)

        if not mod44DataExist:
            raise IOError('250 meters resolution watermask from MODIS '
                    '44W Product does not exist - see Nansat '
                    'documentation to get it (the path is % s)' % mod44path)

        # MOD44W data does exist: open the VRT file in Nansat
        watermask = Nansat(mod44path + '/MOD44W.vrt', mapper='MOD44W',
                           log_level=self.logger.level)
        # reproject on self or given Domain
        if dst_domain is None:
            dst_domain = self
        lon, lat = dst_domain.get_border()
        watermask.crop_lonlat([lon.min(), lon.max()], [lat.min(), lat.max()])
        watermask.reproject(dst_domain, addmask=False, **kwargs)

        return watermask

    def write_figure(self, filename='', bands=1, clim=None, addDate=False,
                     array_modfunc=None, **kwargs):
        """Save a raster band to a figure in graphical format.

        Get numpy array from the band(s) and band information specified
        either by given band number or band id.
        -- If three bands are given, merge them and create PIL image.
        -- If one band is given, create indexed image
        Create Figure object and:
        Adjust the array brightness and contrast using the given min/max or
        histogram.
        Apply logarithmic scaling of color tone.
        Generate and append legend.
        Save the PIL output image in PNG or any other graphical format.
        If the filename extension is 'tif', the figure file is converted
        to GeoTiff

        Parameters
        -----------
        filename : str
            Output file name. if one of extensions 'png', 'PNG', 'tif',
            'TIF', 'bmp', 'BMP', 'jpg', 'JPG', 'jpeg', 'JPEG' is included,
            specified file is created. otherwise, 'png' file is created.
        bands : integer or string or list (elements are integer or string),
            default = 1
            the size of the list has to be 1 or 3.
            if the size is 3, RGB image is created based on the three bands.
            Then the first element is Red, the second is Green,
            and the third is Blue.
        clim : list with two elements or 'hist' to specify range of colormap
            None (default) : min/max values are fetched from WKV,
            fallback-'hist'
            [min, max] : min and max are numbers, or
            [[min, min, min], [max, max, max]]: three bands used
            'hist' : a histogram is used to calculate min and max values
        addDate : boolean
            False (default) : no date will be aded to the caption
            True : the first time of the object will be added to the caption
        array_modfunc : None
            None (default) : figure created using array in provided band
            function : figure created using array modified by provided function
        **kwargs : parameters for Figure().

        Notes
        ---------
        if filename is specified, creates image file

        Returns
        -------
        Figure : Figure object
            filename extension define format (default format is png)

        Example
        --------
            >>> n.write_figure('test.jpg') # write indexed image
            >>> n.write_figure('test_rgb_hist.jpg', clim='hist', bands=[1, 2, 3]) # RGB image
            >>> n.write_figure('r09_log3_leg.jpg', logarithm=True, legend=True,
                                gamma=3, titleString='Title', fontSize=30,
                                numOfTicks=15) # add legend
            >>> n.write_figure(filename='transparent.png', bands=[3],
                               mask_array=wmArray,
                               mask_lut={0: [0,0,0]},
                               clim=[0,0.15], cmapName='gray',
                               transparency=[0,0,0]) # write transparent image
        See also
        --------
        Figure()

        `<http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps>`_

        """
        # convert <bands> from integer, or string, or list of strings
        # into list of integers
        if isinstance(bands, list):
            for i, band in enumerate(bands):
                bands[i] = self.get_band_number(band)
        else:
            bands = [self.get_band_number(bands)]

        # == create 3D ARRAY ==
        array = None
        for band in bands:
            # get array from band and reshape to (1,height,width)
            iArray = self[band]
            if array_modfunc:
                iArray = array_modfunc(iArray)
            iArray = iArray.reshape(1, iArray.shape[0], iArray.shape[1])
            # create new 3D array or append band
            if array is None:
                array = iArray
            else:
                array = np.append(array, iArray, axis=0)

        # == CREATE FIGURE object and parse input parameters ==
        fig = Figure(array, **kwargs)
        array = None

        # == PREPARE cmin/cmax ==
        # check if cmin and cmax are given as the arguments
        if 'cmin' in kwargs.keys() and 'cmax' in kwargs.keys():
            clim = [kwargs['cmin'], kwargs['cmax']]

        # try to get clim from WKV if it is not given as the argument
        # if failed clim will be evaluated from histogram
        if clim is None:
            clim = [[], []]
            for i, iBand in enumerate(bands):
                try:
                    defValue = (self.vrt.dataset.GetRasterBand(iBand).
                                GetMetadataItem('minmax').split(' '))
                except:
                    clim = 'hist'
                    break
                clim[0].append(float(defValue[0]))
                clim[1].append(float(defValue[1]))

        # Estimate color min/max from histogram
        if clim == 'hist':
            clim = fig.clim_from_histogram(**kwargs)

        # modify clim to the proper shape [[min], [max]]
        # or [[min, min, min], [max, max, max]]
        if (len(clim) == 2 and
           ((isinstance(clim[0], float)) or (isinstance(clim[0], int))) and
           ((isinstance(clim[1], float)) or (isinstance(clim[1], int)))):
            clim = [[clim[0]], [clim[1]]]

        # if the len(clim) is not same as len(bands), the 1st element is used.
        for i in range(2):
            if len(clim[i]) != len(bands):
                clim[i] = [clim[i][0]] * len(bands)

        self.logger.info('clim: %s ' % clim)

        # == PREPARE caption ==
        if 'caption' in kwargs:
            caption = kwargs['caption']
        else:
            # get longName and units from vrt
            band = self.get_GDALRasterBand(bands[0])
            longName = band.GetMetadata().get('long_name', '')
            units = band.GetMetadata().get('units', '')

            # make caption from longname, units
            caption = longName + ' [' + units + ']'

        # add DATE to caption
        if addDate:
            caption += self.time_coverage_start.strftime(' %Y-%m-%d')

        self.logger.info('caption: %s ' % caption)

        # == PROCESS figure ==
        fig.process(cmin=clim[0], cmax=clim[1], caption=caption)

        # == finally SAVE to a image file ==
        fig.save(filename, **kwargs)
        # If tiff image, convert to GeoTiff
        if filename[-3:] == 'tif':
            self.vrt.copyproj(filename)
        return fig


    def write_geotiffimage(self, filename, band_id=1):
        """Writes an 8-bit GeoTiff image for a given band.

        The colormap is fetched from the metadata item 'colormap'. Fallback colormap is 'gray'.

        Color limits are fetched from the metadata item 'minmax'. If 'minmax' is not specified, min
        and max of the raster data is used.

        The method can be replaced by using nansat.write_figure(). However, write_figure uses PIL,
        which does not allow Tiff compression. This gives much larger files.

        Parameters
        -----------
        filename : str
        band_id : int or str

        """
        bandNo = self.get_band_number(band_id)
        band = self.get_GDALRasterBand(band_id)

        minmax = band.GetMetadataItem('minmax')
        # Get min and max from band histogram if not given (from wkv)
        if minmax is None:
            (rmin, rmax) = band.ComputeRasterMinMax()
            minmax = str(rmin) + ' ' + str(rmax)

        bMin = float(minmax.split(' ')[0])
        bMax = float(minmax.split(' ')[1])
        # Make colormap from WKV information
        try:
            colormap = band.GetMetadataItem('colormap')
            cmap = cm.get_cmap(colormap, 256)
            cmap = cmap(np.arange(256)) * 255
        except:
            if not MATPLOTLIB_IS_INSTALLED:
                self.logger.debug('Geotiff is only available in gray '
                                  'since matplotlib is not installed.')
            cmap = np.vstack([np.arange(256.),
                              np.arange(256.),
                              np.arange(256.),
                              np.ones(256)*255]).T
        colorTable = gdal.ColorTable()
        for i in range(cmap.shape[0]):
            colorEntry = (int(cmap[i, 0]), int(cmap[i, 1]),
                          int(cmap[i, 2]), int(cmap[i, 3]))
            colorTable.SetColorEntry(i, colorEntry)
        # Write Tiff image, with data scaled to values between 0 and 255
        outDataset = gdal.GetDriverByName('Gtiff').Create(filename,
                                                          band.XSize,
                                                          band.YSize, 1,
                                                          gdal.GDT_Byte,
                                                          ['COMPRESS=LZW'])
        data = self.__getitem__(bandNo)
        scaledData = ((data - bMin) / (bMax - bMin)) * 255
        outDataset.GetRasterBand(1).WriteArray(scaledData)
        outDataset.GetRasterBand(1).SetMetadata(band.GetMetadata())
        try:
            outDataset.GetRasterBand(1).SetColorTable(colorTable)
        except:
            # Happens after reprojection, a possible bug?
            print('Could not set color table')
            print(colorTable)
        outDataset = None
        self.vrt.copyproj(filename)

    @property
    def time_coverage_start(self):
        return parse_time(self.get_metadata('time_coverage_start'))

    @property
    def time_coverage_end(self):
        return parse_time(self.get_metadata('time_coverage_end'))

    def get_metadata(self, key=None, band_id=None, unescape=True):
        """Get metadata from self.vrt.dataset

        Parameters
        ----------
        key : str
            name of the metadata key. If not givem all metadata is returned
        band_id : int or str
            number or name of band to get metadata from.
            If not given, global metadata is returned
        unescape : bool
            Replace '&quot;', '&amp;', '&lt;' and '&gt;' with these symbols " & < > ?

        Returns
        --------
        metadata : str
            string with metadata if key is given and found
        metadata : dict
            dictionary with all metadata if key is not given

        Raises
        ------
        ValueError, if key is not found

        """
        # get all metadata from dataset or from band
        if band_id is None:
            metadata = self.vrt.dataset.GetMetadata()
        else:
            metadata = self.get_GDALRasterBand(band_id).GetMetadata()

        # remove escapes of special characters
        if unescape:
            for i in metadata:
                metadata[i] = saxutils.unescape(metadata[i], {'&quot;': '"'})

        # get all metadata or from a key
        if key is not None:
            try:
                metadata = metadata[key]
            except KeyError:
                raise ValueError('%s does not have metadata %s' % (self.filename, key))

        return metadata

    def set_metadata(self, key='', value='', band_id=None):
        """Set metadata to self.vrt.dataset

        Parameters
        -----------
        key : string or dictionary with strings
            name of the metadata, or dictionary with metadata names, values
        value : string
            value of metadata
        band_id : int or str
            number or name of band
            Without : global metadata is set

        Notes
        ---------
        self.vrt.dataset : sets metadata in GDAL current dataset

        """
        # set all metadata to the dataset or to the band
        if band_id is None:
            metadata_receiver = self.vrt.dataset
        else:
            bandNumber = self.get_band_number(band_id)
            metadata_receiver = self.vrt.dataset.GetRasterBand(bandNumber)

        # set metadata from dictionary or from single pair key,value
        if type(key) == dict:
            for k in key:
                metadata_receiver.SetMetadataItem(str(k), str(key[k]))
        else:
            metadata_receiver.SetMetadataItem(str(key), str(value))

    def _get_dataset_metadata(self):
        # open GDAL dataset. It will be parsed to all mappers for testing
        gdal_dataset, metadata = None, dict()
        if not self.filename.startswith('http'):
            try:
                gdal_dataset = gdal.Open(self.filename)
            except RuntimeError:
                self.logger.error('GDAL could not open %s, trying to read with Nansat mappers...'
                                  % self.filename)
        if gdal_dataset is not None:
            # get metadata from the GDAL dataset
            metadata = gdal_dataset.GetMetadata()

        return gdal_dataset, metadata


    def _get_mapper(self, mappername, **kwargs):
        """Create VRT file in memory (VSI-file) with variable mapping

        If mappername is given only this mapper will be used,
        else loop over all availble mappers in mapperList to get the
        matching one.
        In the loop :
            If the specific error appears the mapper is not used
            and the next mapper is tested.
            Otherwise the mapper returns VRT.
        If type of the sensor is identified, add mapping variables.
        If all mappers fail, make simple copy of the input DS into a VSI/VRT

        Parameters
        -----------
        mappername : string, optional (e.g. 'ASAR' or 'merisL2')

        Returns
        --------
        tmp_vrt : VRT object
            tmp_vrt.dataset is a GDAL VRT dataset

        Raises
        --------
        IOError : occurs if the input file does not exist
        ValueError : occurs if given mapper cannot open the input file
        NansatReadError : occurs if no mapper fits the input file

        """
        if os.path.isfile(self.filename):
            # Make sure file exists and can be opened for reading before proceeding
            assert os.access(self.filename, os.R_OK)
        else:
            for f in glob.glob(os.path.join(self.filename, '*.*')):
                assert os.access(f, os.R_OK)
        # lazy import of nansat mappers
        # if nansat mappers were not imported yet
        global nansatMappers
        if nansatMappers is None:
            nansatMappers = _import_mappers()

        # open GDAL dataset. It will be parsed to all mappers for testing
        gdal_dataset, metadata = self._get_dataset_metadata()
        tmp_vrt = None

        # TODO: There seems to be code repetition in this if-test - should be avoided...
        if mappername is not '':
            # If a specific mapper is requested, we test only this one.
            # get the module name
            mappername = 'mapper_' + mappername.replace('mapper_', '').replace('.py', '').lower()
            # check if the mapper is available
            if mappername not in nansatMappers:
                raise ValueError('Mapper ' + mappername + ' not found')

            # check if mapper is importbale or raise an ImportError error
            if isinstance(nansatMappers[mappername], tuple):
                errType, err, traceback = nansatMappers[mappername]
                # self.logger.error(err, exc_info=(errType, err, traceback))
                # TODO: python 3.6 does not support with syntax
                raise EnvironmentError
                #raise errType, err, traceback

            # create VRT using the selected mapper
            tmp_vrt = nansatMappers[mappername](self.filename, gdal_dataset, metadata, **kwargs)
            self.mapper = mappername.replace('mapper_', '')
        else:
            # We test all mappers, import one by one
            import_errors = []
            for iMapper in nansatMappers:
                # skip non-importable mappers
                if isinstance(nansatMappers[iMapper], tuple):
                    # keep errors to show before use of generic mapper
                    import_errors.append(nansatMappers[iMapper][1])
                    continue

                self.logger.debug('Trying %s...' % iMapper)

                # show all ImportError warnings before trying generic_mapper
                if iMapper == 'mapper_generic' and len(import_errors) > 0:
                    self.logger.error('\nWarning! The following mappers failed:')
                    for ie in import_errors:
                        self.logger.error(import_errors)

                # create a Mapper object and get VRT dataset from it
                try:
                    tmp_vrt = nansatMappers[iMapper](self.filename, gdal_dataset, metadata, **kwargs)
                    self.logger.info('Mapper %s - success!' % iMapper)
                    self.mapper = iMapper.replace('mapper_', '')
                    break
                except WrongMapperError:
                    pass

        # if no mapper fits, make simple copy of the input DS into a VSI/VRT
        if tmp_vrt is None and gdal_dataset is not None:
            self.logger.warning('No mapper fits, returning GDAL bands!')
            tmp_vrt = VRT.from_gdal_dataset(gdal_dataset, metadata=metadata)
            for iBand in range(gdal_dataset.RasterCount):
                tmp_vrt.create_band({'SourceFilename': self.filename,
                                     'SourceBand': iBand + 1})
                tmp_vrt.dataset.FlushCache()
            self.mapper = 'gdal_bands'

        # if GDAL cannot open the file, and no mappers exist which can make VRT
        if tmp_vrt is None and gdal_dataset is None:
            # check if given data file exists
            if not os.path.isfile(self.filename):
                raise IOError('%s: File does not exist' % (self.filename))
            raise NansatReadError('%s: File cannot be read with NANSAT - '
                    'consider writing a mapper' % self.filename)

        return tmp_vrt

    def get_band_number(self, band_id):
        """Return absolute band number

        Check if given band_id is valid
        Return absolute number of the band in the VRT

        Parameters
        ----------
        band_id : int or str or dict
            - if int : checks if such band exists and returns band_id
            - if str : finds band with coresponding name
            - if dict : finds first band with given metadata

        Returns
        --------
        absolute band number : int

        """
        band_number = 0
        # if band_id is str: create simple dict with seraching criteria
        if type(band_id) == str:
            band_id = {'name': band_id}

        # if band_id is dict: search self.bands with seraching criteria
        if type(band_id) == dict:
            bands_meta = self.bands()
            for b in bands_meta:
                num_correct_keys = 0
                for key in band_id:
                    if (key in bands_meta[b] and
                            band_id[key] == bands_meta[b][key]):
                        num_correct_keys = num_correct_keys + 1
                    if num_correct_keys == len(band_id):
                        band_number = b
                        break

        # if band_id is int and with bounds: return this number
        if (type(band_id) == int and band_id >= 1 and
                band_id <= self.vrt.dataset.RasterCount):
            band_number = band_id

        # if no band_number found - raise error
        if band_number == 0:
            raise ValueError('Cannot find band %s! '
                              'band_number is from 1 to %s'
                              % (str(band_id), self.vrt.dataset.RasterCount))

        return band_number

    def get_transect(self, points, bands,
                        lonlat=True,
                        smooth_radius=0,
                        smooth_function=nanmedian,
                        data=None,
                        cornersonly=False):
        """Get values from transect from given vector of poins

        Parameters
        ----------
        points : 2xN list or array, N (number of points) >= 1
            coordinates [[x1, x2, y2], [y1, y2, y3]]
        bands : list of int or string
            elements of the list are band number or band Name
        lonlat : bool
            If the points in lat/lon, then True.
            If the points in pixel/line, then False.
        smooth_radius: int
            If smootRadius is greater than 0, smooth every transect
            pixel as the median or mean value in a circule with radius
            equal to the given number.
        smooth_function: func
            function for averaging values collected within smooth radius
        data : ndarray
            alternative array with data to take values from

        Returns
        --------
        transect : numpy record array

        """
        # check if points is 2D array with shape 2xN (N>=1)
        if (len(np.shape(points)) != 2 or
              np.shape(points)[0] != 2 or
              np.shape(points)[1] < 1):
            # points are not 2xN array
            raise ValueError('Input points must be 2xN array with N>0')

        # get names of bands
        band_names = []
        for band in bands:
            try:
                bandN = self.get_band_number(band)
            except ValueError:
                self.logger.error('Wrong band name %s' % band)
            else:
                band_names.append(self.bands()[bandN]['name'])

        if data is not None:
            band_names.append('input')

        pixVector, linVector = self._get_pix_lin_vectors(points, lonlat, cornersonly, smooth_radius)

        # create output transect
        t = np.recarray((len(pixVector)), dtype=[('pixel', int),
                                                ('line', int),
                                                ('lon', float),
                                                ('lat', float), ])

        # add pixel, line, lon, lat values to output
        t['pixel'] = pixVector
        t['line'] = linVector
        t['lon'], t['lat'] = self.transform_points(t['pixel'], t['line'], DstToSrc=0)

        # get values from bands or input data
        for band_name in band_names:
            t = self._extract_transect_data(t, band_name, data, smooth_radius, smooth_function)

        return t

    def _extract_transect_data(self, t, band_name, data, smooth_radius, smooth_function):
        """Extract data along transect from input band"""
        # mask for extraction within circular area
        xgrid, ygrid = np.mgrid[0:smooth_radius * 2 + 1, 0:smooth_radius * 2 + 1]
        distance = ((xgrid - smooth_radius) ** 2 + (ygrid - smooth_radius) ** 2) ** 0.5
        mask = distance <= smooth_radius

        if band_name == 'input':
            bandArray = data
        else:
            bandArray = self[band_name]
        # average values from pixel inside a circle
        bandValues = []
        for r, c in zip(t['line'], t['pixel']):
            subarray = bandArray[r-smooth_radius:r+smooth_radius+1,
                                 c-smooth_radius:c+smooth_radius+1]
            bandValues.append(smooth_function(subarray[mask]))
        return append_fields(t, band_name, bandValues).data

    def digitize_points(self, band=1, **kwargs):
        """Get coordinates of interactively digitized points

        Parameters
        ----------
        band : int or str
            ID of Nansat band
        **kwargs : keyword arguments for imshow

        Returns
        --------
        points : list
            list of 2xN arrays of points to be used in Nansat.get_transect()

        """
        data = self[band]
        browser = PointBrowser(data, **kwargs)
        points = browser.get_points()

        return points

    def crop_interactive(self, band=1, maxwidth=1000, **kwargs):
        """Interactively select boundary and crop Nansat object

        Parameters
        ----------
        band : int or str
            id of the band to show for interactive selection of boundaries
        maxwidth : int
            large input data is downscaled to <maxwidth>
        **kwargs : keyword arguments for imshow

        Notes
        --------
        self.vrt : VRT
            superVRT is created with modified SrcRect and DstRect

        Returns
        -------
        extent : (x_offset, y_offset, x_size, y_size)
            x_offset  - X offset in the original dataset
            y_offset  - Y offset in the original dataset
            x_size - width of the new dataset
            y_size - height of the new dataset

        Examples
        --------
            >>> extent = n.crop_interactive(band=1) # crop a subimage interactively

        """
        resized = False
        if self.shape()[1] > maxwidth:
            factor = self.resize(width=1000)
            resized = True
        else:
            factor = 1
        # use interactive PointBrowser for selecting extent
        try:
           points = self.digitize_points(band=band,**kwargs)[0]
        except:
           if resized:
              self.undo()
           return

        x_offset, x_size = Nansat._get_crop_offset_size(0, points, factor)
        y_offset, y_size = Nansat._get_crop_offset_size(1, points, factor)

        if resized:
            self.undo()

        return self.crop(x_offset, y_offset, x_size, y_size)

    def crop_lonlat(self, lonlim, latlim):
        """Crop Nansat object to fit into given longitude/latitude limit

        Parameters
        ----------
        lonlim : list of 2 float
            min/max of longitude
        latlim : list of 2 float
            min/max of latitude

        Notes
        --------
        self.vrt : VRT
            crops vrt to size that corresponds to lon/lat limits

        Returns
        -------
        extent : (x_offset, y_offset, x_size, y_size)
            x_offset  - X offset in the original dataset
            y_offset  - Y offset in the original dataset
            x_size - width of the new dataset
            y_size - height of the new dataset

        Examples
        --------
            >>> extent = n.crop(lonlim=[-10,10], latlim=[-20,20]) # crop for given lon/lat limits

        """
        # lon/lat lists for four corners
        lon_corners = [lonlim[0], lonlim[0], lonlim[1], lonlim[1]]
        lat_corners = [latlim[0], latlim[1], latlim[0], latlim[1]]
        pix_corners, lin_corners = self.transform_points(lon_corners, lat_corners, 1)

        offset_size = lambda corners: (round(min(corners)), round(max(corners) - min(corners)))
        x_offset, x_size = offset_size(pix_corners)
        y_offset, y_size = offset_size(lin_corners)

        return self.crop(x_offset, y_offset, x_size, y_size)

    @staticmethod
    def _get_crop_offset_size(axis, points, factor):
        """Get offset and size of cropped image"""
        offset = round(points.min(axis=1)[axis] / factor)
        size = round((points.max(axis=1)[axis] - offset) / factor)
        return offset, size

    @staticmethod
    def _fix_crop_offset_size(offset, size, raster_size, allow_larger=False):
        """Check and correct default values of invalid offset or size"""
        # test if crop is totally outside
        if (offset > raster_size or (offset + size) < 0):
            raise ValueError('Cropping region is outside the image! offset: %f size: %f'
                               %(float(offset), float(size)))
        if allow_larger:
            return offset, size

        if offset < 0:
            size += offset
            offset = 0
        if (size + offset) > raster_size:
            size = raster_size - offset
        return offset, size

    def crop(self, x_offset, y_offset, x_size, y_size, allow_larger=False):
        """Crop Nansat object

        Create superVRT, modify the Source Rectangle (SrcRect) and Destination
        Rectangle (DstRect) tags in the VRT file for each band in order
        to take only part of the original image,
        create new GCPs or new GeoTransform for the cropped object.

        Parameters
        ----------
        x_offset : int
            pixel offset of subimage
        y_offset : int
            line offset of subimage
        x_size : int
            width in pixels of subimage
        y_size : int
            height in pizels of subimage
        allow_larger : bool
            Allow resulting extent to be larger than the original image?

        Notes
        --------
        self.vrt : super-VRT is created with modified SrcRect and DstRect

        Returns
        -------
        extent : (x_offset, y_offset, x_size, y_size)
            x_offset  - X offset in the original dataset
            y_offset  - Y offset in the original dataset
            x_size - width of the new dataset
            y_size - height of the new dataset

        Examples
        --------
            >>> extent = n.crop(10, 20, 100, 200)

        """
        x_offset, x_size = Nansat._fix_crop_offset_size(x_offset, x_size,
                                                        self.shape()[1], allow_larger)
        y_offset, y_size = Nansat._fix_crop_offset_size(y_offset, y_size,
                                                        self.shape()[0], allow_larger)

        extent = (int(x_offset), int(y_offset), int(x_size), int(y_size))
        self.logger.debug('x_offset: %d, y_offset: %d, x_size: %d, y_size: %d' % extent)

        # test if crop is larger or equal to image size
        if x_offset == y_offset == 0 and (y_size, x_size) == self.shape():
            self.logger.error(('WARNING! Cropping region is larger or equal to image!'))
            return extent

        # create super VRT and change it
        self.vrt = self.vrt.get_super_vrt()
        self.vrt.set_offset_size('x', x_offset, x_size)
        self.vrt.set_offset_size('y', y_offset, y_size)
        self.vrt.shift_cropped_gcps(x_offset, x_size, y_offset, y_size)
        self.vrt.shift_cropped_geo_transform(x_offset, x_size, y_offset, y_size)
        return extent

    def extend(self, left=0, right=0, top=0, bottom=0):
        """Extend domain from four sides

        Parameters
        ----------
        left : int
            number of pixels to add from left side
        right : int
            number of pixels to add from right side
        top : int
            number of pixels to add from top side
        bottom : int
            number of pixels to add from bottom side

        Notes
        -----
        Canges self.vrt by adding nexgative offset or setting size to be large that original size.

        """
        x_offset, y_offset = -left, -top
        x_size = self.shape()[1] + left + right
        y_size = self.shape()[0] + top + bottom
        self.crop(x_offset, y_offset, x_size, y_size, allow_larger=True)

    def _get_pix_lin_vectors(self, points, lonlat, cornersonly, smooth_radius):
        """Get vectors with pixel and line values for input corner points"""
        # if points in degree, convert them into pix/lin
        if lonlat:
            pix, lin = self.transform_points(points[0], points[1], DstToSrc=1)
        else:
            pix, lin = points[0], points[1]

        if cornersonly:
            return pix, lin

        # full vectors of pixel coordinates based on coordinates of vertices
        pixVector, linVector = [pix[0]], [lin[0]]
        for pn in range(len(pix[1:])):
            px0, px1 = pix[pn], pix[pn+1]
            py0, py1 = lin[pn], lin[pn+1]
            length = np.round(np.hypot(px1-px0, py0-py1))
            pixVector += list(np.linspace(px0, px1, length+1)[1:])
            linVector += list(np.linspace(py0, py1, length+1)[1:])

        # remove out of region points
        pixVector = np.floor(pixVector)
        linVector = np.floor(linVector)
        gpi = ((pixVector >= (0 + smooth_radius)) *
               (linVector >= (0 + smooth_radius)) *
               (pixVector < (self.shape()[1] - smooth_radius)) *
               (linVector < (self.shape()[0] - smooth_radius)))

        return pixVector[gpi], linVector[gpi]


def _import_mappers(log_level=None):
    """Import available mappers into a dictionary

    Returns
    --------
    nansat_mappers : dict
        key  : mapper name
        value: class Mapper(VRT) from the mappers module

    """
    logger = add_logger('import_mappers', logLevel=log_level)
    # import built-in mappers
    import nansat.mappers
    mapper_packages = [nansat.mappers]

    # import user-defined mappers (if any)
    try:
        import nansat_mappers as nansat_mappers_pkg
    except ImportError:
        pass
    else:
        logger.info('User defined mappers found in %s' % nansat_mappers_pkg.__path__)
        mapper_packages = [nansat_mappers_pkg, nansat.mappers]

    # create ordered dict for mappers
    nansat_mappers = OrderedDict()
    for mapper_package in mapper_packages:
        logger.debug('From package: %s' % mapper_package.__path__)
        # scan through modules and load all modules that contain class Mapper
        for finder, name, ispkg in (pkgutil.iter_modules(mapper_package.__path__)):
            logger.debug('Loading mapper %s' % name)
            # Only mappers containing 'mapper' in the module name should be returned
            if not 'mapper' in name: continue
            loader = finder.find_module(name)
            # try to import mapper module
            module = loader.load_module(name)
            # add the imported mapper to nansat_mappers
            if hasattr(module, 'Mapper'):
                nansat_mappers[name] = module.Mapper

        # move netcdfcdf mapper to the end
        if 'mapper_netcdf_cf' in nansat_mappers:
            nansat_mappers['mapper_netcdf_cf'] = nansat_mappers.pop('mapper_netcdf_cf')

        # move generic_mapper to the end
        if 'mapper_generic' in nansat_mappers:
            nansat_mappers['mapper_generic'] = nansat_mappers.pop('mapper_generic')

    return nansat_mappers
