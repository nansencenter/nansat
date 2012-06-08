# Name:    nansat.py
# Purpose: main file of the NANSAT module.
#
# Authors:      Asuka Yamakava, Anton Korosov, Knut-Frode Dagestad
#
# Created:     29.06.2011
# Copyright:   (c) NERSC 2012
# Licence:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details:
# http://www.gnu.org/licenses/

import dateutil.parser

from domain import Domain
from vrt import *
from figure import *
from nansat_tools import add_logger, Node


class GDALError(Error):
    '''Error from GDAL '''
    pass


class DataError(Error):
    '''Error for data.
        e.g. : empty pixel value array in get_pixelValueRange()'''
    pass

class ProjectionError(Error):
    '''Error for reprojection.'''
    pass

class Nansat(Domain):
    '''Main of Nansat

    Construct Nansat object that consist of
        basic dataset information (fileName, dataset, metadata etc..),
        VRT file which points to orignal file with satellite data and
        is saved in an XML format in memory (GDAL VSI).
    '''
    def __init__(self, fileName="", mapperName="", domain=None,
                 array=None, wkv="", parameters={}, logLevel=30):
        '''Construct Nansat object

        Open GDAL dataset,
        Read metadata,
        Generate GDAL VRT file with mapping of variables in memory
        Create logger

        Parameters
        ----------
        fileName : string
            location of the file
        mapperName : string, optional
            "ASAR", "hurlam", "merisL1", "merisL2", "ncep", "radarsat2",
            "seawifsL2" are currently available.  (27.01.2012)
        domain : domain object
        array : numpy array
        parameters: dictionary
            metadata for array band
        logLevel: int, optional, default: logging.DEBUG (30)
            Level of logging. See: http://docs.python.org/howto/logging.html

        Creates
        --------
        self.mapperList: list of file names
            list of available working mappers
        self.fileName : file name
            set file name given by the argument
        self.raw : VRT object
            set VRT object with VRT dataset with mapping of variables
        self.vrt : VRT object
            Copy of self.raw
        self.logger: logging.Logger
            logger for output debugging info
        self.name: string
            name of object (for writing KML)

        Raises
        ------
            GDALError: occurs when the dataset is None or "".

        '''
        # checkt the arguments
        if fileName=="" and (domain is None or array is None):
            raise OptionError("Either fileName or (domain and array) is required.")

        # SET CONSTANTS
        # create logger
        self.logger = add_logger(logName='Nansat', logLevel=logLevel)

        # all available mappers
        self.mapperList = [
              'mapper_ASAR.py',
              'mapper_hirlam.py',
              'mapper_merisL1.py',
              'mapper_merisL2.py',
              'mapper_modisL1.py',
              'mapper_ncep.py',
              'mapper_radarsat2.py',
              'mapper_MOD44W.py',
              'mapper_modisL2NRT.py'
              ]

        self.logger.debug('Mappers: ' + str(self.mapperList))

        # set input file name
        self.fileName = fileName

        # create self from a file using mapper or...
        if fileName != "":
            # Get oroginal VRT object with mapping of variables
            self.raw = self._get_mapper(mapperName)
            # Set current VRT object
            self.vrt = self.raw.copy()
        # create using array, domain, and parameters
        else:
            # Get vrt from domain
            self.vrt = VRT(gdalDataset=domain.vrt.dataset)
            # add a band from array
            self.add_band(array, wkv, parameters)
            # Set raw VRT object
            self.raw = self.vrt.copy()

        # name, for compatibility with some Domain methods
        self.name = self.fileName

        self.logger.debug('Object created from %s ' % self.fileName)

    def __getitem__(self, bandNo):
        ''' Returns the band as a NumPy array, by overloading []
        Parameters:
        -----------
            bandNo: int or string
                If int, array from band with number <bandNo> is returned
                If string, array from band with metadata 'band_name' equal to
                <bandNo> is returned
        Returns
        -------
            self.get_GDALRasterBand(bandNo).ReadAsArray(): NumPy array

        '''
        # get band
        band = self.get_GDALRasterBand(bandNo)
        # get scale and offset
        scale = float(band.GetMetadata().get('scale', '1'))
        offset = float(band.GetMetadata().get('offset', '0'))
        # get data
        bandData = band.ReadAsArray()
        # perform scaling if necessary
        if scale != 1 or offset != 0:
            bandData = bandData * scale + offset

        return bandData

    def __repr__(self):
        '''Creates string basic info about the Nansat object

        '''
        outString = self.fileName + '\n'
        outString += '-' * 40 + '\n'
        outString += self.list_bands(False)
        outString += '-' * 40 + '\n'
        outString += Domain.__repr__(self)
        return outString

    def add_band(self, array, wkv="", p=None):
        '''Add band from the array to self.vrt
        
        Create VRT object which contains VRT and RAW binary file
        Create new band in self.vrt whci points to this vrt

        Parameters
        ----------
            array : Numpy array with band data
            wkv : string with name of WKV
            p: dictionary for defining parameters of the band, name, etc.

        Modifies
        --------
            Creates VRT object with VRT-file and RAW-file
            Adds band to the self.vrt

        '''
        # if p is not given it is empty dictionary
        # NB: default value p={} remembers the value from the last assigment
        if p is None:
            p = {}
        # create a band from array
        arrayVrt = VRT(array=array)
        # add the array band into self.vrt
        self.vrt._create_band(arrayVrt.fileName, 1, wkv, p)
        # delete arrayVrt and unnecessary files
        #delFileList = [arrayVrt.fileName, arrayVrt.fileName.replace(".vrt", ".raw")]
        #arrayVrt._del_gdalMemoryFile(delFileList)

    def export(self, fileName):
        '''Create a netCDF file

        Parameters
        ----------
            fileName: output file name

        Modifies
        --------
            Create a netCDF file

            self.vrt.dataset : VRT dataset of VRT object
                dataType in each VRTRasterBand is modified from GDAL dataType
                to netCDF dataType. The subelements "Offset" and "Scale"
                are added for each VRTRasterBand if need be.

        !! NB !!
        --------
            If number of bands is more than one,
            serial numbers are added at the end of each band name.

            It is possible to fix it by changing
            line.4605 in GDAL/frmts/netcdf/netcdfdataset.cpp :

            "if( nBands > 1 ) sprintf(szBandName,"%s%d",tmpMetadata,iBand);"
            --> "if( nBands > 1 ) sprintf(szBandName,"%s",tmpMetadata);"

        '''
        # Get the node from the XML content
        node0 = Node.create(self.vrt.read_xml())

        # Change the element from GDAL datatype to NetCDF data type
        for iBand in node0.nodeList("VRTRasterBand"):
            dataType = iBand.getAttribute("dataType")
            dataType = {"UInt16": "Int16", "CInt16": "Int16",
                        "UInt32": "Int32",
                        "CFloat32": "Float32", "CFloat64": "Float64"
                        }.get(dataType, dataType)
            iBand.replaceAttribute("dataType", dataType)

        # Copy the vrt and overwrite the modified element
        tmpVrt = self.vrt.copy()
        tmpVrt.write_xml(str(node0.rawxml()))

        # Get projection and add it in the global metadata
        try:
            tmpVrt.dataset.SetMetadataItem("gcpProjection",
                                            self.vrt.dataset.
                                            GetGCPProjection().
                                            replace('"', "&quot;"))
        except:
            self.logger.warning("Unable to get projection for netcdf")

        # Replace the band names
        for iBand in range(tmpVrt.dataset.RasterCount):
            band = tmpVrt.dataset.GetRasterBand(iBand + 1)
            try:
                band.SetMetadataItem('NETCDF_VARNAME',
                                     str(band.GetMetadata_Dict()["band_name"]))
            except:
                self.logger.warning('Unable to set NETCDF_VARNAME for band %d'
                                    % iBand)

        # Create a netCDF file
        dataset = gdal.GetDriverByName("netCDF").CreateCopy(fileName,
                                                            tmpVrt.dataset)

    def resize(self, factor=1, width=None, height=None, method="average"):
        '''Proportional resize of the dataset.

        The dataset is resized as (xSize*factor, ySize*factor) or
        (width, calulated height) or (calculated width, height).
        self.vrt is rewritten to the the downscaled sizes.
        If GCPs are given in a dataset, they are also rewritten.
        If resize() is called without any parameters resizing/reprojection is
        cancelled.


        Parameters
        ----------
        Either factor, or width, or height should be given:
            factor: float, optional, default=1
            width: int, optional
            height: int, optional
            method: "average" (default) or "subsample" (= nearest neighbor),
                    optional
        Modifies
        --------
            self.vrt.dataset : VRT dataset of VRT object
                raster size are modified to downscaled size.
                If GCPs are given in the dataset, they are also overwritten.
        Raises
        ------
            OptionError: occurs when method is not "average" or "subsample"

        '''
        # resize back to original size/setting
        if factor == 1 and width is None and height is None:
            self.vrt = self.raw.copy()
            return

        # estimate factor if width or height is given
        if width is not None:
            factor = float(width) / float(self.vrt.dataset.RasterXSize)
        if height is not None:
            factor = float(height) / float(self.vrt.dataset.RasterYSize)

        if not (method == "average" or method == "subsample"):
            raise OptionError("method should be 'average' or 'subsample'")

        # Get XML content from VRT-file
        vrtXML = self.vrt.read_xml()

        node0 = Node.create(vrtXML)
        rasterXSize = int(float(node0.getAttribute("rasterXSize")) * factor)
        rasterYSize = int(float(node0.getAttribute("rasterYSize")) * factor)
        self.logger.info('New size/factor: (%f, %f)/%f' %
                        (rasterXSize, rasterYSize, factor))
        node0.replaceAttribute("rasterXSize", str(rasterXSize))
        node0.replaceAttribute("rasterYSize", str(rasterYSize))

        for iNode in node0.nodeList("VRTRasterBand"):
            iNode.node("DstRect").replaceAttribute("xSize", str(rasterXSize))
            iNode.node("DstRect").replaceAttribute("ySize", str(rasterYSize))
            # if method="average", overwrite "SimpleSource" to "AveragedSource"
            if method == "average":
                iNode.replaceTag("SimpleSource", "AveragedSource")

        # Edit GCPs to correspond to the downscaled size
        if node0.node("GCPList"):
            for iNode in node0.node("GCPList").nodeList("GCP"):
                pxl = float(iNode.getAttribute("Pixel")) * factor
                if pxl > float(rasterXSize):
                    pxl = rasterXSize
                iNode.replaceAttribute("Pixel", str(pxl))
                lin = float(iNode.getAttribute("Line")) * factor
                if lin > float(rasterYSize):
                    lin = rasterYSize
                iNode.replaceAttribute("Line", str(lin))

        # Write the modified elemements into VRT
        self.vrt.write_xml(str(node0.rawxml()))

    def get_GDALRasterBand(self, bandNo=1, bandID=None):
        ''' Get a GDALRasterBand of a given Nansat object.

        Get a GDALRasterBand specified by the argument.

        If a bandID is given, secify a bandNo based on it.
        Otherwise check if the given bandNo is proper.
        Get a GDALRasterBand from vrt.

        Parameters
        ----------
            bandNo: serial number or string, optional (default is 1)
                if number - a band number of the band to fetch
                if string bandID = {'band_name': bandNo}
            bandID: a dictionary with metadata unique for one band

        Returns
        -------
            self.vrt.dataset.GetRasterBand: a GDAL RasterBand

        Raises
        ------
            OptionError: occurs when the bandNo is not a proper number.

        Example
        -------
            b = get_GDALRasterBand(1)
            b = get_GDALRasterBand('sigma0')
            b = get_GDALRasterBand({"short_name":"radiance",
                                    "wavelength":"1240"})
        See Also
        --------
            _specify_bandNo: specify a band number based on the bandID list

        '''
        # If bandID is given, bandNo is specified here.
        if bandID is not None:
            bandNo = self._specify_bandNo(bandID)
        # if bandNo is given and it is string fetch band which has
        # band_name == bandNo
        elif isinstance(bandNo, str):
            bandNo = self._specify_bandNo({'band_name': bandNo})
        # if given bandNo is over the existing bands, give error message
        elif (bandNo < 1 or bandNo > self.vrt.dataset.RasterCount):
            raise OptionError("Nansat.get_GDALRasterBand(): "
                             "bandNo takes from 1 to",
                             self.vrt.dataset.RasterCount)
        # Based on bandNo,
        # the GDAL RasterBand of the corresponding band is returned
        return self.vrt.dataset.GetRasterBand(bandNo)

    def list_bands(self, doPrint=True):
        ''' Show band information of the given Nansat object

        Show serial number, longName, name and all parameters
        for each band in the metadata of the given Nansat object.

        Parameters:
        -----------
            doPrint: boolean, optional, default=True
                do print, otherwise it is returned as string
        Returns:
        --------
            outString: String
                formatted string with bands info
        '''
        outString = ''
        for iBand in range(self.vrt.dataset.RasterCount):
            band = self.vrt.dataset.GetRasterBand(iBand + 1)
            metadata = band.GetMetadata()
            outString += "Band : %d\n" % (iBand + 1)
            for i in metadata:
                outString += "  %s: %s\n" % (i, band.GetMetadataItem(i))
        if doPrint:
            print outString
        else:
            return outString

    def reproject(self, dstDomain=None, resamplingAlg=0):
        ''' Reproject the object based on the given Domain

        Warp the raw VRT using AutoCreateWarpedVRT() using projection
        from the Domain.
        Modify XML content of the warped vrt using the Domain parameters.
        Generate warpedVRT and replace self.vrt with warpedVRT.

        Parameters
        ----------
            dstDomain: domain
                destination Domain where projection and resolution are set
            resamplingAlg: int
                0, 1 or 2 stand for NearestNeigbour, Bilinear, Cubic

        Modifies
        --------
            self.vrt: VRT object with VRT dataset
                replaced to warpedVRT dataset

        Raises
        ------
            ProjectionError: occurs when the projection of the target data
            is None.
            AttributeError: occurs when it is impossible to get warpedVRT.

        See Also
        --------
            http://www.gdal.org/gdalwarp.html

        '''
        # Get source SRS (either Projection or GCPProjection)
        srcWKT = self.raw.dataset.GetProjection()
        if srcWKT == '':
            srcWKT = self.raw.dataset.GetGCPProjection()

        if srcWKT == '':
            raise ProjectionError("Nansat.reproject(): "
                                  "rawVrt.GetProjection() is None")

        if dstDomain is None:
            # dereproject
            self.vrt = self.raw.copy()
        else:
            # warp the Raw VRT onto the coordinate stystem given by
            # input Domain
            warpedVRT = gdal.AutoCreateWarpedVRT(
                                   self.raw.dataset, srcWKT,
                                   dstDomain.vrt.dataset.GetProjection(),
                                   resamplingAlg)

            # set default vrt to be the warped one
            warpedVRT = VRT(vrtDataset=warpedVRT, logLevel=self.logger.level)

            # modify extent of the created Warped VRT
            warpedVRTXML = warpedVRT.read_xml()
            warpedVRTXML = self._modify_warped_VRT_XML(warpedVRTXML,
                                       dstDomain.vrt.dataset.RasterXSize,
                                       dstDomain.vrt.dataset.RasterYSize,
                                       dstDomain.vrt.dataset.GetGeoTransform())
            warpedVRT.write_xml(warpedVRTXML)

            # check created Warped VRT
            if warpedVRT.dataset is None:
                raise AttributeError("Nansat.reproject():cannot get warpedVRT")

            self.vrt = warpedVRT

    def reproject_on_gcp(self, gcpImage, resamplingAlg=0):
        ''' Reproject the object onto the input object with gcps
        NB! This is a test function required urgently for the open-wind
        project. It is tesed only on NCEP or GLOBAL DEM and
        RADARSAT2 or MERIS images and should be refined and
        added to Nansat.reproject()

        Parameters
        ----------
            gcpImage: Nansat object of an image with GCPs
            resamplingAlg: integer, option for AutoCreateWarpedVRT

        Modifies
        --------
            self.vrt: VRT object with VRT dataset
                replaced with warpedVRT

        '''
        dstDataset = gcpImage.vrt.dataset

        # prepare pure lat/lon WKT
        latlongWKT = self._latlong_srs().ExportToWkt()

        # get source SRS (either Projection or GCPProjection)
        srcWKT = self.vrt.dataset.GetProjection()
        if srcWKT == '':
            srcWKT = self.vrt.dataset.GetGCPProjection()

        # the transformer converts lat/lon to pixel/line of SRC image
        srcTransformer = gdal.Transformer(
                             self.vrt.dataset, None,
                             ['SRC_SRS=' + srcWKT,
                             'DST_SRS=' + latlongWKT])

        # get GCPs from DST image
        gcps = dstDataset.GetGCPs()

        # create 'fake' GCPs
        for g in gcps:
            # transform DST lat/lon to SRC pixel/line
            succ, point = srcTransformer.TransformPoint(1, g.GCPX, g.GCPY)
            srcPixel = point[0]
            srcLine = point[1]

            # swap coordinates in GCPs:
            # pix1/line1 -> lat/lon  =>=>  pix2/line2 -> pix1/line1
            g.GCPX = g.GCPPixel
            g.GCPY = g.GCPLine
            g.GCPPixel = srcPixel
            g.GCPLine = srcLine

        # make copy of the RAW VRT file and replace GCPs
        tmpVRT = self.raw.copy()

        # create 'fake' STEREO projection for 'fake' GCPs of SRC image
        srsString = ("+proj=stere +lon_0=0 +lat_0=0 +k=1 "
                     "+ellps=WGS84 +datum=WGS84 +no_defs ")
        stereoSRS = osr.SpatialReference()
        stereoSRS.ImportFromProj4(srsString)
        stereoSRSWKT = stereoSRS.ExportToWkt()
        tmpVRT.dataset.SetGCPs(gcps, stereoSRSWKT)
        tmpVRT.dataset.SetProjection('')

        # remove GeoTransfomr from SRC image
        # read XML content from VRT
        tmpVRTXML = tmpVRT.read_xml()
        # find and remove GeoTransform
        node0 = Node.create(tmpVRTXML)
        node1 = node0.delNode("GeoTransform")

        # Write the modified elemements back into temporary VRT
        tmpVRT.write_xml(str(node0.rawxml()))

        # create warped vrt out of tmp vrt
        rawWarpedVRT = gdal.AutoCreateWarpedVRT(tmpVRT.dataset, stereoSRSWKT,
                                                stereoSRSWKT,
                                                resamplingAlg)
        self.logger.debug('rawWarpedVRT: %s' % str(rawWarpedVRT))
        warpedVRT = VRT(vrtDataset=rawWarpedVRT, logLevel=self.logger.level)
        warpeVRTXML = warpedVRT.read_xml()

        # change size and geotransform to fit the DST image
        warpeVRTXML = self._modify_warped_VRT_XML(warpeVRTXML,
                                   dstDataset.RasterXSize,
                                   dstDataset.RasterYSize,
                                   (0, 1, 0, 0, 0, 1))
        # update warped VRT with new XML
        warpedVRT.write_xml(warpeVRTXML)

        #append GCPs and Projection from the dstImage
        warpedVRT.dataset.SetGCPs(dstDataset.GetGCPs(),
                                  dstDataset.GetGCPProjection())
        warpedVRT.dataset.SetProjection('')
        # set current VRT to be warped
        self.vrt = warpedVRT

    def watermask(self, mod44path=None, dstDomain=None):
        ''' Create numpy array with watermask (water=1, land=0)

        250 meters resolution watermask from MODIS 44W Product:
        http://www.glcf.umd.edu/data/watermask/

        Watermask is stored as tiles in TIF(LZW) format and a VRT file
        All files are stored in one directory.
        A tarball with compressed TIF and VRT files should be additionally
        downloaded from the Nansat wiki:
        https://svn.nersc.no/nansat/wiki/Nansat/Data/Watermask

        The method:
            Gets the directory either from input parameter or from environment
            variable MOD44WPATH
            Open Nansat object from the VRT file
            Reprojects the watermask onto the current object using reproject()
            or reproject_on_jcps()
            Returns the reprojected Nansat object

        Parameters:
        -----------
            mod44path : string, optional, default=None
                path with MOD44W Products and a VRT file

        Returns:
        --------
            watermask : Nansat object with water mask in current projection

        See also:
        ---------
            250 meters resolution watermask from MODIS 44W Product:
            http://www.glcf.umd.edu/data/watermask/
        '''
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
            # MOD44W data does not exist generate empty matrix
            watermask = np.zeros(self.vrt.dataset.RasterXSize,
                                 self.vrt.dataset.RasterYSize)
        else:
            # MOD44W data does exist: open the VRT file in Nansat
            watermask = Nansat(mod44path + '/MOD44W.vrt',
                               logLevel=self.logger.level)
            # choose reprojection method
            if dstDomain is not None:
                watermask.reproject(dstDomain)
            elif self.vrt.dataset.GetGCPCount() == 0:
                watermask.reproject(Domain(self.vrt.dataset,
                                    logLevel=self.logger.level))
            else:
                watermask.reproject_on_gcp(self)

        return watermask

    def write_figure(self, fileName=None, bands=1, clim=None, **kwargs):

        ''' Save a raster band to a figure in grapfical format.

        Get numpy array from the band(s) and band information specified
        either by given band number or band id.
        -- If three bands are given, merge them and create PIL image.
        -- If one band is given, create indexed image
        Create Figure object and:
        Adjust the array brightness and contrast using the given min/max or
        histogram.
        Apply logarithmic scaling of color tone.
        Generate and append legend.
        Save the PIL output image in PNG or any other graphical fornat.

        Parameters
        ----------
            fileName: string, optional
                Output file name. if one of extensions "png", "PNG", "tif",
                "TIF", "bmp", "BMP", "jpg", "JPG", "jpeg", "JPEG" is included,
                specified file is crated. otherwise, "png" file is created.
                if None, the figure object is returned.
            bands : int or list, default = 1
                the size of the list has to be 1 or 3.
                if the size is 3, RGB image is created based on the
                three bands.
                Then the first element is Red, the second is Green,
                and the third is Blue.
            clim : list with two elements or 'hist' to specify range of
                colormap
                None (default): min/max values are fetched from WKV,
                fallback-'hist'
                [min, max] : min and max are numbers, or
                [[min, min, min], [max, max, max]]: three bands used
                'hist' : a histogram is used to calculate min and max values
            **kwargs : parameters for Figure(). See figure.Figure()

        Modifies
        --------
            if fileName is specified, creates image file

        Returns
        ------
            Figure object

        Example:
        --------
        #write only indexed image, color limits from WKV or from histogram
        n.write_figure('test.jpg')
        #write only RGB image, color limits from histogram
        n.write_figure('test_rgb_hist.jpg', clim='hist', bands=[1, 2, 3])
        #write indexed image, apply log scaling and gamma correction,
        #add legend and type in title 'Title', increase font size and put 15
        tics
        n.write_figure('r09_log3_leg.jpg', logarithm=True, legend=True,
                                gamma=3, titleString='Title', fontSize=30,
                                numOfTicks=15)

        See also
        ------
        Figure()
        http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps

        '''
        # == create 3D ARRAY ==
        if isinstance(bands, list):
            array = np.array([])
            for iBand in bands:
                iArray = self.__getitem__(iBand)
                self.logger.debug(iArray.shape)
                array = np.append(array, iArray)
            array = array.reshape(len(bands),
                                 int(self.vrt.dataset.RasterYSize),
                                 int(self.vrt.dataset.RasterXSize))
        else:
            array = self.__getitem__(bands)
            bands = [bands]

        # == CREATE FIGURE object and parse input parameters ==
        fig = Figure(array, **kwargs)
        array = None

        # == PREPARE cmin/cmax ==
        # try to get clim from WKV if it is not given
        # if failed clim will be evaluated from histogram
        if clim is None:
            clim = [[], []]
            for i, iBand in enumerate(bands):
                try:
                    defValue = (self.vrt.dataset.GetRasterBand(iBand).
                                GetMetadataItem("minmax").split(" "))
                except:
                    clim = 'hist'
                    break
                clim[0].append(float(defValue[0]))
                clim[1].append(float(defValue[1]))

        # Estimate color min/max from histogram
        if clim == 'hist':
            clim = fig.clim_from_histogram()

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
        # get longName and units from vrt
        bandNo = 1
        if isinstance(bands[0], int):
            bandNo = bands[0]
        elif isinstance(bands[0], str):
            bandNo = self._specify_bandNo({'band_name': bands[0]})

        longName = self.vrt.dataset.GetRasterBand(bandNo).GetMetadataItem(
                                                                   "long_name")
        units = self.vrt.dataset.GetRasterBand(bandNo).GetMetadataItem("units")
        # if they don't exist make empty strings
        if longName is None:
            longName = ''
        if units is None:
            units = ''
        caption = longName + ' [' + units + ']'
        self.logger.info('caption: %s ' % caption)

        # == PROCESS figure ==
        fig.process(cmin=clim[0], cmax=clim[1], caption=caption)

        # == finally SAVE to a image file ==
        if fileName is not None:
            fig.save(fileName)

        return fig

    def get_time(self, bandNo=None):
        ''' Get time for dataset and/or its bands

        Parameters
        ----------
            bandNo: integer (default = None)

        If bandNo is given, or if all bands have the same time,
        a single datetime object is returned.
        Othwerwise a list of datetime objects for each band is returned.

        '''
        time = []
        for i in range(self.vrt.dataset.RasterCount):
            band = self.vrt.dataset.GetRasterBand(i + 1)
            try:
                time.append(dateutil.parser.parse(
                                band.GetMetadataItem("time")))
            except:
                self.logger.debug("Band " + str(i + 1) + " has no time")
                time.append(None)

        if bandNo is not None:
            return time[bandNo - 1]
        else:
            if len(set(time)) == 1:
                return time[0]
            else:
                return time

    def get_metadata(self, key=None, bandNo=None):
        ''' Get metadata from self.vrt.dataset

        Parameters:
        -----------
        key: string, optional
            name of the metadata key. If not givem all metadata is returned
        bandNo: int, optional
            number of band to get metadata from. If not given, global metadata
            is returned

        Returns:
            a string with metadata if key is given and found
            an empty string if key is given and not found
            a dictionary with all metadata if key is not given
        '''

        # get all metadata from dataset or from band
        if bandNo is None:
            metadata = self.vrt.dataset.GetMetadata()
        else:
            metadata = self.vrt.dataset.GetRasterBand(bandNo).GetMetadata()

        # get all metadata or from a key
        if key is not None:
            metadata = metadata.get(key, '')

        return metadata

    def set_metadata(self, key='', value='', bandNo=None):
        ''' Set metadata to self.raw.dataset and self.vrt.dataset

        Parameters:
        -----------
        key: string or dictionary with strings
            name of the metadata, or dictionary with metadata names, values
        value: string
            value of metadata
        bandNo: int
            number of band to set metadata to. Without: global metadata is set

        Modifies:
        ---------
            self.raw.dataset: sets metadata in GDAL raw dataset
            self.vrt.dataset: sets metadata in GDAL current dataset
            '''

        # set all metadata from dataset or from band
        if bandNo is None:
            metaReceiverRAW = self.raw.dataset
            metaReceiverVRT = self.vrt.dataset
        else:
            metaReceiverRAW = self.raw.dataset.GetRasterBand(bandNo)
            metaReceiverVRT = self.vrt.dataset.GetRasterBand(bandNo)

        # set metadata from dictionary or from single pair key,value
        if isinstance(key, dict):
            for k in key:
                metaReceiverRAW.SetMetadataItem(k, key[k])
                metaReceiverVRT.SetMetadataItem(k, key[k])
        else:
            metaReceiverRAW.SetMetadataItem(key, value)
            metaReceiverVRT.SetMetadataItem(key, value)

    def _get_mapper(self, mapperName):
        ''' Create VRT file in memory (VSI-file) with variable mapping

        If mapperName is given, it is added as the first in the self.mapperList
        Loop over all availble mappers in mapperList to get the matching one.
        In the loop:
            If the specific error appears the mapper is not used
            and the next mapper is tested.
            Otherwise the mapper returns VRT.
        If type of the sensor is identified, add mapping variables.
        If all mapper do not fit, simply copy the input DS into a VSI/VRT

        Parameters
        ----------
        mapperName : string, optional
            "ASAR", "hurlam", "merisL1", "merisL2", "ncep", "radarsat2",
            "seawifsL2" are currently available.  (27.01.2012)

        Returns
        -------
            tmpVRT : VRT object
                tmpVRT.dataset is a GDAL VRT dataset

        Raises
        --------
            GDALError: occures if the given file cannot be opened with GDAL
            TypeError: occurs when the given driver type is not registarated
                        in the mappers.

        '''
        # open GDAL dataset. It will be parsed to all mappers for testing
        gdalDataset = gdal.Open(self.fileName)
        if (gdalDataset is None) or (gdalDataset == ""):
            raise GDALError("Nansat._get_,apper(): Cannot get the dataset from"
                            + self.fileName)
        # gete metadata from the GDAL dataset
        metadata = gdalDataset.GetMetadata()

        # add the given mapper first
        self.mapperList = ['mapper_' + mapperName] + self.mapperList

        # try to import and get VRT datasaet from all mappers. Break on success
        # if none of the mappers worked - None is returned
        tmpVRT = None
        # For debugging:
        """
        mapper_module = __import__('mapper_modisL2NRT')
        tmpVRT = mapper_module.Mapper(self.fileName, gdalDataset,
                                      metadata, logLevel=self.logger.level)
        """
        # Otherwise
        for iMapper in self.mapperList:
            try:
                #get rid of .py extension
                iMapper = iMapper.replace('.py', '')
                self.logger.debug('Trying %s...' % iMapper)
                #import mapper
                mapper_module = __import__(iMapper)
                #create a Mapper object and get VRT dataset from it
                tmpVRT = mapper_module.Mapper(self.fileName, gdalDataset,
                                              metadata,
                                              logLevel=self.logger.level)
                self.logger.info('Mapper %s - success!' % iMapper)
                break
            except:
                pass
        # """
        # if no mapper fits, make simple copy of the input DS into a VSI/VRT
        if tmpVRT is None:
            self.logger.info('No mapper fits!')
            tmpVRT = VRT(vrtDataset=gdalDataset, logLevel=self.logger.level)

        return tmpVRT

    def _get_pixelValue(self, val, defVal):
        if val == "":
            return defVal
        else:
            return val

    def _modify_warped_VRT_XML(self, vrtXML,
                          rasterXSize, rasterYSize, geoTransform):
        ''' Modify rasterXsize, rasterYsize and geotranforms in the warped VRT

        Parameters
        ----------
            vrtXML: string with XML
                original content of the original VRT file
            rasterXSize: integer
                desired X size of warped image
            rasterYSize: integer
                desired Y size of warped image
            geoTransform: tuple of 6 integers
                desired GeoTransform size of the warped image

        Returns
        --------
            vrtXML : string with XML
                modified content of the original VRT file

        '''
        invGeotransform = gdal.InvGeoTransform(geoTransform)
        node0 = Node.create(vrtXML)

        node0.replaceAttribute("rasterXSize", str(rasterXSize))
        node0.replaceAttribute("rasterYSize", str(rasterYSize))

        # convert proper string style and set to the GeoTransform element
        node0.node("GeoTransform").value = str(geoTransform).strip("()")
        node0.node("DstGeoTransform").value = str(geoTransform).strip("()")
        node0.node("DstInvGeoTransform").value = str(
                                                invGeotransform[1]).strip("()")

        if node0.node("SrcGeoLocTransformer"):
            node0.node("BlockXSize").value = str(rasterXSize)
            node0.node("BlockYSize").value = str(rasterYSize)

        return str(node0.rawxml())

    def _specify_bandNo(self, bandID):
        '''Specify a band number based on bandID {'key_name': 'key_value'}

        Compare the key values of the bandID to the values of the
        metadata in each band.
        If matched, append the band number (iBand) into candidate list.
        If not, go to the next band.
        Iterate these steps until all bands are checked.
        If single band is specified at the end, return the band number.
        Otherwise raise OptionError.

        Parameters
        ----------
            bandID: a dictionary
                Parameters to specify single band
                (e.g. {"name":"radiance", "wavelength":"1234"})

        Returns
        -------
            candidate[0]+1 : a band number

        Raises
        ------
            OptionError: occurs when there is no band which satisfies
            the condition (bandID) or when there are several bands chosen
            by the condition.

        '''
        # search for the specific band based on bandID
        candidate = []
        # loop through all bands
        for iBand in range(self.vrt.dataset.RasterCount):
            # get metadata from band
            bandMetadata = self.raw.dataset.GetRasterBand(
                                                    iBand + 1).GetMetadata()
            allKeysAreGood = True
            # loop through all input keys
            for bandIDKey in bandID:
                # if band doesn't have key from input or value doesn't match
                if (bandIDKey not in bandMetadata or
                        bandMetadata[bandIDKey] != bandID[bandIDKey]):
                    allKeysAreGood = False

            if allKeysAreGood:
                candidate.append(iBand)

        # if a band is specified, set it to bandNo.
        # if some bands are chosen, give an error message and stop.
        if len(candidate) == 1:
            self.logger.info("You chose bandNo: %d" % (candidate[0] + 1))
            return candidate[0] + 1
        elif len(candidate) >= 2:
            raise OptionError("Nansat._specify_bandNo(): "
                              "Cannot specify a single band "
                              "by the given arguments")
        else:
            raise OptionError("Nansat._specify_bandNo(): "
                              "Cannot find any band by the given arguments")
