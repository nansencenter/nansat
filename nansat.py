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

# import standard and additional libraries
from nansat_tools import *


# import nansat parts
try:
    from domain import Domain
except ImportError:
    warnings.warn('Cannot import Domain!'
                  'Nansat will not work.')

try:
    from figure import Figure
except ImportError:
    warnings.warn('Cannot import Figure!'
                  'Nansat will not work.')

try:
    from vrt import VRT
except ImportError:
    warnings.warn('Cannot import VRT!'
                  'Nansat will not work.')

try:
    from nansatmap import Nansatmap
except ImportError:
    warnings.warn('Cannot import Nansatmap!'
                  'Nansat will not work.')


try:
    from nansatshape import Nansatshape
except ImportError:
    warnings.warn('Cannot import NansatOGR!'
                  'Nansat will not work.')

# Force GDAL to raise exceptions
try:
    gdal.UseExceptions()
except:
    warnings.warn('GDAL will not raise exceptions.'
                  'Probably GDAL is not installed')

# Set environment variables, the script directory
nansathome = os.path.dirname(os.path.abspath(inspect.getfile(
                                             inspect.currentframe())))
sys.path.append(nansathome)
sys.path.append(nansathome + '/mappers/')
if not 'GDAL_DRIVER_PATH' in os.environ:
    os.environ['GDAL_DRIVER_PATH'] = nansathome + '/pixelfunctions/'

# Compile pixelfunctions if not already done.
if sys.platform.startswith('win'):
    if not os.path.exists(nansathome + '/pixelfunctions/gdal_PIXFUN.DLL'):
        print 'Cannot find "gdal_PIXFUN.dll". Compile pixelfunctions !!'
else:
    if not os.path.exists(nansathome + '/pixelfunctions/gdal_PIXFUN.so'):
        print 'Cannot find "gdal_PIXFUN.so". Compiling pixelfunctions...'
        os.system('cd ' + nansathome + '/pixelfunctions/; make clean; make')

class Nansat(Domain):
    '''Container for geospatial data, performs all high-level operations

    n = Nansat(fileName) opens the file with satellite or model data for
    reading, adds scientific metadata to bands, and prepares the data for
    further handling.

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
    '''

    def __init__(self, fileName='', mapperName='', domain=None,
                 array=None, parameters=None, logLevel=30, **kwargs):
        '''Create Nansat object

        if <fileName> is given:
            Open GDAL dataset,
            Read metadata,
            Generate GDAL VRT file with mapping of variables in memory
            Create logger
            Create Nansat object for perfroming high-level operations
        if <domain> and <array> are given:
            Create VRT object from data in <array>
            Add geolocation from <domain>

        Parameters
        -----------
        fileName : string
            location of the file
        mapperName : string, optional
            name of the mapper from nansat/mappers dir. E.g.
            'ASAR', 'hirlam', 'merisL1', 'merisL2', etc.
        domain : Domain object
            Geo-reference of a new raster
        array : numpy array
            Firts band of a new raster
        parameters : dictionary
            Metadata for the 1st band of a new raster,e.g. name, wkv, units,...
        logLevel : int, optional, default: logging.DEBUG (30)
            Level of logging. See: http://docs.python.org/howto/logging.html

        Parameters (**kwargs)
        ---------------------
        -- VRT
        eResampleAlg = 0
        use_geolocationArray = True
        use_gcps = True
        use_geotransform = True
        WorkingDataType = None
        tps = False
        blockSize = None
        -- mapper_envisat
        envisat_zoomSize = 500
        envisat_step = 1
        -- mapper_asar
        asar_geolocation = False
        -- mapper_meris
        meris_geolocation = True
        -- mapper_obpg_l2
        obpg_l2_GCP_COUNT = 10
        -- mapper_pathfinder52
        pathfinder52_minQual = 4
        -- mapper_viirs_l1
        viirs_l1_GCP_COUNT0 = 5
        viirs_l1_GCP_COUNT1 = 20
        viirs_l1_pixelStep = 1
        viirs_l1_lineStep = 1
        -- mapper_aster_l1a
        aster_l1a_bandNames = ['VNIR_Band1', 'VNIR_Band2', 'VNIR_Band3N']
        aster_l1a_bandWaves = [560, 660, 820]
        -- mapper_case2reg
        case2regKwargs_wavelengths = [None, 413, 443, 490, 510, 560, 620, 665, 681, 709, 753, None, 778, 864]}
        -- mapper_generic
        generic_rmMetadatas = ['NETCDF_VARNAME', '_Unsigned', 'ScaleRatio',
                               'ScaleOffset', 'dods_variable']

        Creates
        --------
        self.mapperList : list of file names
            list of available working mappers
        self.fileName : file name
            set file name given by the argument
        self.raw : Mapper(VRT) object
            set VRT object with VRT dataset with mapping of variables
        self.vrt : Mapper(VRT) object
            Copy of self.raw
        self.logger : logging.Logger
            logger for output debugging info
        self.name : string
            name of object (for writing KML)

        '''
        # check the arguments
        if fileName == '' and domain is None:
            raise OptionError('Either fileName or domain is required.')

        # create logger
        self.logger = add_logger('Nansat', logLevel)

        # empty dict of VRTs with added bands
        self.addedBands = {}

        # add all available mappers if mapperName is not given
        self.mapper = 'None'
        self.mapperList = []
        if mapperName is '':
            for folder in sys.path:
                for mapper in glob.glob(folder + '/mapper_*.py'):
                    self.mapperList.append(os.path.basename(mapper))

            # pop and append generic mapper to the end
            self.mapperList.pop(self.mapperList.index('mapper_generic.py'))
            self.mapperList.append('mapper_generic.py')

        self.logger.debug('Mappers: ' + str(self.mapperList))

        # set input file name
        self.fileName = fileName
        # name, for compatibility with some Domain methods
        self.name = os.path.basename(fileName)
        self.path = os.path.dirname(fileName)

        # create self.raw from a file using mapper or...
        if fileName != '':
            # Make original VRT object with mapping of variables
            self.raw = self._get_mapper(mapperName, **kwargs)
            # Set current VRT object
            self.vrt = self.raw.copy()
        # ...create using array, domain, and parameters
        else:
            # Get vrt from domain
            self.raw = VRT(gdalDataset=domain.vrt.dataset)
            # Set current VRT object
            self.vrt = VRT(gdalDataset=domain.vrt.dataset)
            if array is not None:
                # add a band from array
                self.add_band(array=array, parameters=parameters)

        self.logger.debug('Object created from %s ' % self.fileName)

    def __getitem__(self, bandID):
        ''' Returns the band as a NumPy array, by overloading []

        Parameters
        -----------
        bandID : int or str
            If int, array from band with number <bandID> is returned
            If string, array from band with metadata 'name' equal to
            <bandID> is returned

        Returns
        --------
        self.get_GDALRasterBand(bandID).ReadAsArray() : NumPy array

        '''
        # get band
        band = self.get_GDALRasterBand(bandID)
        # get expression from metadata
        expression = band.GetMetadata().get('expression', '')
        # get data
        bandData = band.ReadAsArray()
        # execute expression if any
        if expression != '':
            bandData = eval(expression)

        return bandData

    def __repr__(self):
        '''Creates string with basic info about the Nansat object'''

        outString = '-' * 40 + '\n'
        outString += self.fileName + '\n'
        outString += '-' * 40 + '\n'
        outString += 'Mapper: ' + self.mapper + '\n'
        outString += '-' * 40 + '\n'
        outString += self.list_bands(False)
        outString += '-' * 40 + '\n'
        outString += Domain.__repr__(self)
        return outString

    def add_band(self, fileName=None, vrt=None, bandID=1, array=None,
                 parameters=None, resamplingAlg=1):
        '''Add band from the array to self.vrt

        Create VRT object which contains VRT and RAW binary file and append it
        to self.addedBands
        Create new band in self.raw which points to this vrt

        NB : Adding band is possible for raw (nonprojected, nonresized) images
        only. Adding band will cancel any previous reproject() or resize().

        Parameters
        -----------
        fileName : string, name of the file, source of band
        vrt : VRT, source of band
        bandID : int, number of the band in fileName or in vrt
        array : Numpy array with band data
        parameters : dictionary, band metadata: wkv, name, etc.
        resamplingAlg : 0, 1, 2 stands for nearest, bilinear, cubic

        Modifies
        ---------
        Creates VRT object with VRT-file and RAW-file
        Adds band to the self.vrt

        '''
        # None => {} in input p
        if parameters is None:
            parameters = {}

        # default added vrt, source bandNumber and metadata
        vrt2add = None
        bandNumber = None
        p2add = {}

        # get band from input file name
        if fileName is not None:
            # create temporary nansat object
            n = Nansat(fileName)
            # reproject onto current grid
            n.reproject(self)
            # get vrt to be added
            vrt2add = n.vrt
            # get band metadata
            bandNumber = n._get_band_number(bandID)
            p2add = n.get_metadata(bandID=bandID)

        # get band from input VRT
        if vrt is not None:
            # get VRT to be added
            vrt2add = vrt
            # get band metadata
            bandNumber = bandID
            p2add = vrt.dataset.GetRasterBand(bandID).GetMetadata()

        # get band from input array
        if array is not None:
            if array.shape == self.shape():
                # create VRT from array
                vrt2add = VRT(array=array)
            else:
                # create VRT from resized array
                srcVRT = VRT(array=array)
                vrt2add = srcVRT.resized(self.shape()[1],
                                         self.shape()[0],
                                         resamplingAlg)
            # set parameters
            bandNumber = 1

        # add parameters from input
        for pKey in parameters:
            p2add[pKey] = parameters[pKey]

        # add the array band into self.vrt and get bandName
        bandName = self.raw._create_band({'SourceFilename': vrt2add.fileName,
                                          'SourceBand': bandNumber}, p2add)
        # add VRT with the band to the dictionary
        # (not to loose the VRT object and VRT file in memory)
        self.addedBands[bandName] = vrt2add
        self.raw.dataset.FlushCache()  # required after adding bands
        # copy raw VRT object to the current vrt
        self.vrt = self.raw.copy()

    def bands(self):
        ''' Make a dictionary with all bands metadata

        Returns
        --------
        b : dictionary
            key = N, value = dict with all band metadata

        '''
        b = {}
        for iBand in range(self.vrt.dataset.RasterCount):
            b[iBand + 1] = self.get_metadata(bandID=iBand + 1)

        return b

    def export(self, fileName, rmMetadata=[], addGeolocArray=True,
               addGCPs=True, driver='netCDF'):
        '''Export Nansat object into netCDF or GTiff file

        Parameters
        -----------
        fileName : output file name
        rmMetadata : list with metadata names to remove before export.
            e.g. ['name', 'colormap', 'source', 'sourceBands']
        addGeolocArray : Boolean, add geolocation array datasets? [True].
        addGCPs : Boolean, add GCPs? [True]
        driver : Which GDAL driver (format) to use [netCDF]

        Modifies
        ---------
        Create a netCDF file

        !! NB
        ------
        If number of bands is more than one,
        serial numbers are added at the end of each band name.

        It is possible to fix it by changing
        line.4605 in GDAL/frmts/netcdf/netcdfdataset.cpp :
        'if( nBands > 1 ) sprintf(szBandName,"%s%d",tmpMetadata,iBand);'
        --> 'if( nBands > 1 ) sprintf(szBandName,"%s",tmpMetadata);'

        CreateCopy fails in case the band name has special characters,
        e.g. the slash in 'HH/VV'.

        '''
        # temporary VRT for exporting
        exportVRT = self.vrt.copy()
        exportVRT.real = []
        exportVRT.imag = []

        # Find complex data band
        complexBands = []
        node0 = Node.create(exportVRT.read_xml())
        for iBand in node0.nodeList('VRTRasterBand'):
            dataType = iBand.getAttribute('dataType')
            if dataType[0] == 'C':
                complexBands.append(int(iBand.getAttribute('band')))

        # if data includes complex data,
        # create two bands from real and imaginary data arrays
        if len(complexBands) != 0:
            for i in complexBands:
                bandMetadataR = self.get_metadata(bandID=i)
                bandMetadataR.pop('dataType')
                try:
                    bandMetadataR.pop('PixelFunctionType')
                except:
                    pass
                # Copy metadata and modify 'name' for real and imag bands
                bandMetadataI = bandMetadataR.copy()
                bandMetadataR['name'] = bandMetadataR.pop('name')+'_real'
                bandMetadataI['name'] = bandMetadataI.pop('name')+'_imag'
                # Create bands from the real and imaginary numbers
                exportVRT.real.append(VRT(array=self[i].real))
                exportVRT.imag.append(VRT(array=self[i].imag))

                metaDict = [{'src': {'SourceFilename': exportVRT.real[-1].fileName,
                             'SourceBand':  1},
                             'dst': bandMetadataR},
                            {'src': {'SourceFilename': exportVRT.imag[-1].fileName,
                             'SourceBand':  1},
                             'dst': bandMetadataI}]
                exportVRT._create_bands(metaDict)
            # delete the complex bands
            exportVRT.delete_bands(complexBands)

        # add bands with geolocation arrays to the VRT
        if addGeolocArray and len(exportVRT.geolocationArray.d) > 0:
            exportVRT._create_band(
                {'SourceFilename': self.vrt.geolocationArray.d['X_DATASET'],
                 'SourceBand': int(self.vrt.geolocationArray.d['X_BAND'])},
                {'wkv': 'longitude',
                 'name': 'GEOLOCATION_X_DATASET'})
            exportVRT._create_band(
                {'SourceFilename': self.vrt.geolocationArray.d['Y_DATASET'],
                 'SourceBand': int(self.vrt.geolocationArray.d['Y_BAND'])},
                {'wkv': 'latitude',
                 'name': 'GEOLOCATION_Y_DATASET'})

        # add GCPs to VRT metadata
        if addGCPs:
            exportVRT._add_gcp_metadata()

        # add projection metadata
        srs = self.vrt.dataset.GetProjection()
        exportVRT.dataset.SetMetadataItem('NANSAT_Projection',
                                          srs.replace(',',
                                                      '|').replace('"', '&'))

        # add GeoTransform metadata
        geoTransformStr = str(self.vrt.dataset.GetGeoTransform()).replace(',',
                                                                          '|')
        exportVRT.dataset.SetMetadataItem('NANSAT_GeoTransform',
                                          geoTransformStr)

        # manage metadata for each band
        for iBand in range(exportVRT.dataset.RasterCount):
            band = exportVRT.dataset.GetRasterBand(iBand + 1)
            bandMetadata = band.GetMetadata()
            # set NETCDF_VARNAME
            try:
                bandMetadata['NETCDF_VARNAME'] = bandMetadata['name']
            except:
                self.logger.warning('Unable to set NETCDF_VARNAME for band %d'
                                    % (iBand + 1))
            # remove unwanted metadata from bands
            for rmMeta in rmMetadata:
                try:
                    bandMetadata.pop(rmMeta)
                except:
                    self.logger.info('Unable to remove metadata'
                                     '%s from band %d' % (rmMeta, iBand + 1))
            band.SetMetadata(bandMetadata)

        # remove unwanted global metadata
        globMetadata = exportVRT.dataset.GetMetadata()
        for rmMeta in rmMetadata:
            try:
                globMetadata.pop(rmMeta)
            except:
                self.logger.info('Global metadata %s not found' % rmMeta)
        exportVRT.dataset.SetMetadata(globMetadata)

        # Create an output file using GDAL
        self.logger.debug('Exporting to %s using %s...' % (fileName, driver))
        dataset = gdal.GetDriverByName(driver).CreateCopy(fileName,
                                                          exportVRT.dataset)
        self.logger.debug('Export - OK!')

    def resize(self, factor=1, width=None, height=None, eResampleAlg=-1):
        '''Proportional resize of the dataset.

        The dataset is resized as (xSize*factor, ySize*factor) or
        (width, calulated height) or (calculated width, height).
        self.vrt is rewritten to the the downscaled sizes.
        If GCPs are given in a dataset, they are also rewritten.
        If resize() is called without any parameters then previsous
        resizing/reprojection cancelled.

        Parameters
        -----------
        Either factor, or width, or height should be given:
            factor : float, optional, default=1
            width : int, optional
            height : int, optional
            eResampleAlg : int (GDALResampleAlg), optional
                -1 : Average,
                0 : NearestNeighbour,
                1 : Bilinear,
                2 : Cubic,
                3 : CubicSpline,
                4 : Lancoz
                if eResampleAlg > 0 : VRT.resized() is used
                (Although the default is -1 (Average),
                 if fileName start from 'ASA_', the default is 0 (NN).)

        Modifies
        ---------
        self.vrt.dataset : VRT dataset of VRT object
            raster size are modified to downscaled size.
            If GCPs are given in the dataset, they are also overwritten.

        '''
        # if fileName start from 'ASA_' and eResampleAlg is default (Average),
        # then change eResampleAlg to 0 (NearestNeighbour)
        fileName = self.fileName.split('/')[-1].split('\\')[-1]
        if fileName.startswith('ASA_') and eResampleAlg == -1:
            eResampleAlg = 0

        # resize back to original size/setting
        if factor == 1 and width is None and height is None:
            self.vrt = self.raw.copy()
            return

        # get current shape
        rasterYSize = float(self.shape()[0])
        rasterXSize = float(self.shape()[1])

        # estimate factor if width or height is given
        if width is not None:
            factor = float(width) / rasterXSize
        if height is not None:
            factor = float(height) / rasterYSize

        # calculate new size
        newRasterYSize = int(rasterYSize * factor)
        newRasterXSize = int(rasterXSize * factor)

        self.logger.info('New size/factor: (%f, %f)/%f' %
                        (newRasterXSize, newRasterYSize, factor))

        if eResampleAlg > 0:
            # apply affine transformation using reprojection
            self.vrt = self.vrt.resized(newRasterXSize,
                                        newRasterYSize,
                                        eResampleAlg=eResampleAlg)
        else:
            # simply modify VRT rasterX/Ysize and GCPs
            # Get XML content from VRT-file
            vrtXML = self.vrt.read_xml()
            node0 = Node.create(vrtXML)

            # replace rasterXSize in <VRTDataset>
            node0.replaceAttribute('rasterXSize', str(newRasterXSize))
            node0.replaceAttribute('rasterYSize', str(newRasterYSize))

            # replace xSize in <DstRect> of each source
            for iNode1 in node0.nodeList('VRTRasterBand'):
                for sourceName in ['ComplexSource', 'SimpleSource']:
                    for iNode2 in iNode1.nodeList(sourceName):
                        iNodeDstRect = iNode2.node('DstRect')
                        iNodeDstRect.replaceAttribute('xSize',
                                                      str(newRasterXSize))
                        iNodeDstRect.replaceAttribute('ySize',
                                                      str(newRasterYSize))
                # if method=-1, overwrite 'ComplexSource' to 'AveragedSource'
                if eResampleAlg == -1:
                    iNode1.replaceTag('ComplexSource', 'AveragedSource')
                    iNode1.replaceTag('SimpleSource', 'AveragedSource')

            # Edit GCPs to correspond to the downscaled size
            if node0.node('GCPList'):
                for iNode in node0.node('GCPList').nodeList('GCP'):
                    pxl = float(iNode.getAttribute('Pixel')) * factor
                    if pxl > float(rasterXSize):
                        pxl = rasterXSize
                    iNode.replaceAttribute('Pixel', str(pxl))
                    lin = float(iNode.getAttribute('Line')) * factor
                    if lin > float(rasterYSize):
                        lin = rasterYSize
                    iNode.replaceAttribute('Line', str(lin))

            # Write the modified elemements into VRT
            self.vrt.write_xml(str(node0.rawxml()))


    def get_GDALRasterBand(self, bandID=1):
        ''' Get a GDALRasterBand of a given Nansat object

        If str is given find corresponding band number
        If int is given check if band with this number exists.
        Get a GDALRasterBand from vrt.

        Parameters
        -----------
        bandID : serial number or string, optional (default is 1)
            if number - a band number of the band to fetch
            if string bandID = {'name': bandID}

        Returns
        --------
        GDAL RasterBand

        Example
        -------
        b = n.get_GDALRasterBand(1)
        b = n.get_GDALRasterBand('sigma0')

        '''
        # get band number
        bandNumber = self._get_band_number(bandID)
        # the GDAL RasterBand of the corresponding band is returned
        return self.vrt.dataset.GetRasterBand(bandNumber)

    def list_bands(self, doPrint=True):
        ''' Show band information of the given Nansat object

        Show serial number, longName, name and all parameters
        for each band in the metadata of the given Nansat object.

        Parameters
        -----------
        doPrint : boolean, optional, default=True
            do print, otherwise it is returned as string

        Returns
        --------
        outString : String
            formatted string with bands info

        '''
        # get dictionary of bands metadata
        bands = self.bands()
        outString = ''

        for b in bands:
            # print band number, name
            outString += 'Band : %d %s\n' % (b, bands[b].get('name', ''))
            # print band metadata
            for i in bands[b]:
                outString += '  %s: %s\n' % (i, bands[b][i])
        if doPrint:
            # print to screeen
            print outString
        else:
            return outString

    def reproject(self, dstDomain=None, eResampleAlg=0, blockSize=None,
                  WorkingDataType=None, tps=False):
        ''' Reproject the object based on the given Domain

        Warp the raw VRT using AutoCreateWarpedVRT() using projection
        from the Domain.
        Modify XML content of the warped vrt using the Domain parameters.
        Generate warpedVRT and replace self.vrt with warpedVRT.

        Parameters
        -----------
        dstDomain : domain
            destination Domain where projection and resolution are set
        eResampleAlg : int (GDALResampleAlg)
            0 : NearestNeighbour
            1 : Bilinear
            2 : Cubic,
            3 : CubicSpline
            4 : Lancoz

        Modifies
        ---------
        self.vrt : VRT object with VRT dataset
            replaced to warpedVRT dataset

        See Also
        ---------
        http://www.gdal.org/gdalwarp.html

        '''
        # dereproject
        self.vrt = self.raw.copy()

        # if no domain: quit
        if dstDomain is None:
            return

        # get projection of destination dataset
        dstSRS = dstDomain.vrt.dataset.GetProjection()

        # get destination GCPs
        dstGCPs = dstDomain.vrt.dataset.GetGCPs()
        if len(dstGCPs) > 0:
            # get projection of destination GCPs
            dstSRS = dstDomain.vrt.dataset.GetGCPProjection()

        # create Warped VRT
        warpedVRT = self.raw.create_warped_vrt(
                    dstSRS=dstSRS, dstGCPs=dstGCPs,
                    eResampleAlg=eResampleAlg,
                    xSize=dstDomain.vrt.dataset.RasterXSize,
                    ySize=dstDomain.vrt.dataset.RasterYSize,
                    blockSize=blockSize,
                    geoTransform=dstDomain.vrt.dataset.GetGeoTransform(),
                    WorkingDataType=WorkingDataType,
                    tps=tps)

        # set current VRT object
        self.vrt = warpedVRT
        # add metadata from RAW to VRT (except fileName)
        vrtFileName = self.vrt.dataset.GetMetadataItem('fileName')
        rawMetadata = self.raw.dataset.GetMetadata()
        self.vrt.dataset.SetMetadata(rawMetadata)
        self.vrt.dataset.SetMetadataItem('fileName', vrtFileName)

    def watermask(self, mod44path=None, dstDomain=None):
        ''' Create numpy array with watermask (water=1, land=0)

        250 meters resolution watermask from MODIS 44W Product:
        http://www.glcf.umd.edu/data/watermask/

        Watermask is stored as tiles in TIF(LZW) format and a VRT file
        All files are stored in one directory.
        A tarball with compressed TIF and VRT files should be additionally
        downloaded from the Nansat wiki:
        https://svn.nersc.no/nansat/wiki/Nansat/Data/Watermask

        The method :
            Gets the directory either from input parameter or from environment
            variable MOD44WPATH
            Open Nansat object from the VRT file
            Reprojects the watermask onto the current object using reproject()
            or reproject_on_jcps()
            Returns the reprojected Nansat object

        Parameters
        -----------
        mod44path : string, optional, default=None
            path with MOD44W Products and a VRT file

        Returns
        --------
        watermask : Nansat object with water mask in current projection

        See also
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
            watermaskArray = np.zeros([self.vrt.dataset.RasterXSize,
                                      self.vrt.dataset.RasterYSize])
            watermask = Nansat(domain=self, array=watermaskArray)
        else:
            # MOD44W data does exist: open the VRT file in Nansat
            watermask = Nansat(mod44path + '/MOD44W.vrt', mapperName='MOD44W',
                               logLevel=self.logger.level)
            # reproject on self or given Domain
            if dstDomain is None:
                watermask.reproject(self)
            else:
                watermask.reproject(dstDomain)

        return watermask

    def write_figure(self, fileName=None, bands=1, clim=None, addDate=False,
                    **kwargs):
        ''' Save a raster band to a figure in graphical format.

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
        fileName : string, optional
            Output file name. if one of extensions 'png', 'PNG', 'tif',
            'TIF', 'bmp', 'BMP', 'jpg', 'JPG', 'jpeg', 'JPEG' is included,
            specified file is crated. otherwise, 'png' file is created.
            if None, the figure object is returned.
            if True, the figure is shown
        bands : integer or string or list (elements are integer or string),
            default = 1
            the size of the list has to be 1 or 3.
            if the size is 3, RGB image is created based on the three bands.
            Then the first element is Red, the second is Green,
            and the third is Blue.
        clim : list with two elements or 'hist' to specify range of colormap
            None (default) : min/max values are fetched from WKV,fallback-'hist'
            [min, max] : min and max are numbers, or
            [[min, min, min], [max, max, max]]: three bands used
            'hist' : a histogram is used to calculate min and max values
        addDate : boolean
            False (default) : no date will be aded to the caption
            True : the first time of the object will be added to the caption
        **kwargs : parameters for Figure(). See below:
        ---------- Figure.__init__() parameters: -----------
            cmin : number (int ot float) or [number, number, number]
                0, minimum value of varibale in the matrix to be shown
            cmax : number (int ot float) or [number, number, number]
                1, minimum value of varibale in the matrix to be shown
            gamma : float, >0
                2, coefficient for tone curve udjustment
            subsetArraySize : int
                100000, size of the subset array which is used to get histogram
            numOfColor : int
                250, number of colors for use of the palette.
                254th is black and 255th is white.
            cmapName : string
                'jet', name of Matplotlib colormaps
                see --> http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
            ratio : float, [0 1]
                1.0, ratio of pixels which are used to write the figure
            numOfTicks : int
                5, number of ticks on a colorbar
            titleString : string
                '', title of legend (1st line)
            caption : string
                '', caption of the legend (2nd line, e.g. long name and units)
            fontSize : int
                12, size of the font of title, caption and ticks
            logarithm : boolean, defult = False
                If True, tone curve is used to convert pixel values.
                If False, linear.
            legend : boolean, default = False
                if True, information as textString, colorbar, longName and
                units are added in the figure.
            mask_array : 2D numpy array, int, the shape should be equal
                array.shape. If given this array is used for masking land,
                clouds, etc on the output image. Value of the array are
                indeces. LUT from mask_lut is used for coloring upon this
                indeces.
            mask_lut : dictionary
                Look-Up-Table with colors for masking land, clouds etc. Used
                tgether with mask_array:
                {0, [0,0,0], 1, [100,100,100], 2: [150,150,150], 3: [0,0,255]}
                index 0 - will have black color
                      1 - dark gray
                      2 - light gray
                      3 - blue
            logoFileName : string
                name of the file with logo
            logoLocation : list of two int, default = [0,0]
                X and Y offset of the image
                If positive - offset is from left, upper edge
                If Negative - from right, lower edge
                Offset is calculated from the entire image legend inclusive
            logoSize : list of two int
                desired X,Y size of logo. If None - original size is used
            latGrid : numpy array
                full size array with latitudes. For adding lat/lon grid lines
            lonGrid : numpy array
                full size array with longitudes. For adding lat/lon grid lines
            latlonGridSpacing : int
                number of lat/lon grid lines to show
            latlonLabels : int
                number of lat/lon labels to show along each side.
            transparency : int
                transparency of the image background(mask), set for PIL alpha
                mask in Figure.save()
            default : None

            Advanced parameters
            --------------------
            LEGEND_HEIGHT : float, [0 1]
                0.1, legend height relative to image height
            CBAR_HEIGHTMIN : int
                5, minimum colorbar height, pixels
            CBAR_HEIGHT : float, [0 1]
                0.15,  colorbar height relative to image height
            CBAR_WIDTH : float [0 1]
                0.8, colorbar width  relative to legend width
            CBAR_LOCATION_X : float [0 1]
                0.1, colorbar offset X  relative to legend width
            CBAR_LOCATION_Y : float [0 1]
                0.5,  colorbar offset Y  relative to legend height
            CBAR_LOCATION_ADJUST_X : int
                5,  colorbar offset X, pixels
            CBAR_LOCATION_ADJUST_Y : int
                3,  colorbar offset Y, pixels
            TEXT_LOCATION_X : float, [0 1]
                0.1, caption offset X relative to legend width
            TEXT_LOCATION_Y : float, [0 1]
                0.1, caption offset Y relative to legend height
            NAME_LOCATION_X : float, [0 1]
                0.1, title offset X relative to legend width
            NAME_LOCATION_Y :
                0.3, title  offset Y relative to legend height
            DEFAULT_EXTENSION : string
                '.png'
        --------------------------------------------------

        Modifies
        ---------
        if fileName is specified, creates image file

        Returns
        -------
        Figure object

        Example
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
        # write an image to png with transparent Mask set to color
        transparency=[0,0,0], following PIL alpha mask
        n.write_figure(fileName='transparent.png', bands=[3],
               mask_array=wmArray,
               mask_lut={0: [0,0,0]},
               clim=[0,0.15], cmapName='gray', transparency=[0,0,0])

        See also
        --------
        Figure()
        http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps

        '''
        # convert <bands> from integer, or string, or list of strings
        # into list of integers
        if isinstance(bands, list):
            for i, band in enumerate(bands):
                bands[i] = self._get_band_number(band)
        else:
            bands = [self._get_band_number(bands)]

        # if cmin / cmax is scalar, convert to a list
        for ikey in ['cmin', 'cmax']:
            if ikey in kwargs:
                if not isinstance(kwargs[ikey], list):
                    kwargs[ikey] = [kwargs[ikey]]

        # == create 3D ARRAY ==
        array = None
        for band in bands:
            # get array from band and reshape to (1,height,width)
            iArray = self[band]
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
        # try to get clim from WKV if it is not given
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
            caption += self.get_time()[0].strftime(' %Y-%m-%d')

        self.logger.info('caption: %s ' % caption)

        # modify clim if cmin and cmax are given as the arguments
        if 'cmin' in kwargs.keys():
            for i in range(len(kwargs['cmin'])):
                if kwargs['cmin'][i] > clim[0][i]:
                    clim[0][i] = kwargs['cmin'][i]

        if 'cmax' in kwargs.keys():
            for i in range(len(kwargs['cmax'])):
                if kwargs['cmax'][i] < clim[1][i]:
                    clim[1][i] = kwargs['cmax'][i]

        # == PROCESS figure ==
        fig.process(cmin=clim[0], cmax=clim[1], caption=caption)

        # == finally SAVE to a image file or SHOW ==
        if fileName is not None:
            if type(fileName) == bool and fileName:
                try:
                    if __IPYTHON__:
                        from matplotlib.pyplot import imshow, show
                        from numpy import array
                        sz = fig.pilImg.size
                        image = array(fig.pilImg.im)
                        if fig.pilImg.getbands() == ('P',):
                            image.resize(sz[0], sz[1])
                        elif fig.pilImg.getbands() == ('R', 'G', 'B'):
                            image.resize(sz[0], sz[1], 3)
                        imshow(image)
                        show()
                    else:
                        fig.pilImg.show()
                except:
                    fig.pilImg.show()
            elif type(fileName) == str:
                fig.save(fileName, **kwargs)
                # If tiff image, convert to GeoTiff
                if fileName[-3:] == 'tif':
                    self.vrt.copyproj(fileName)

        return fig

    def write_geotiffimage(self, fileName, bandID=1):
        ''' Writes an 8-bit GeoTiff image for a given band.

        The output GeoTiff image is convenient e.g. for display in a GIS tool.
        Colormap is fetched from the metadata item 'colormap'.
            Fallback colormap is 'jet'.
        Color limits are fetched from the metadata item 'minmax'.
            If 'minmax' is not specified, min and max of raster is used.

        The method can be replaced by using nansat.write_figure(),
        however, write_figure uses PIL which does not allow
        Tiff compression, giving much larger files

        Parameters
        -----------
        fileName : string
        bandID : integer or string(default = 1)

        '''
        bandNo = self._get_band_number(bandID)
        band = self.get_GDALRasterBand(bandID)
        minmax = band.GetMetadataItem('minmax')
        # Get min and max from band histogram if not given (from wkv)
        if minmax is None:
            (rmin, rmax) = band.ComputeRasterMinMax(1)
            minmax = str(rmin) + ' ' + str(rmax)

        bMin = float(minmax.split(' ')[0])
        bMax = float(minmax.split(' ')[1])
        # Make colormap from WKV information
        try:
            colormap = band.GetMetadataItem('colormap')
        except:
            colormap = 'jet'
        try:
            cmap = cm.get_cmap(colormap, 256)
            cmap = cmap(arange(256)) * 255
            colorTable = gdal.ColorTable()
            for i in range(cmap.shape[0]):
                colorEntry = (int(cmap[i, 0]), int(cmap[i, 1]),
                              int(cmap[i, 2]), int(cmap[i, 3]))
                colorTable.SetColorEntry(i, colorEntry)
        except:
            print 'Could not add colormap; Matplotlib may not be available.'
        # Write Tiff image, with data scaled to values between 0 and 255
        outDataset = gdal.GetDriverByName('Gtiff').Create(fileName,
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
            print 'Could not set color table'
            print colorTable
        outDataset = None
        self.vrt.copyproj(fileName)

    def write_nansatmap( self, fileName=None, contour=None, contourf=None,
                         quiver=None, mesh=None, color_bar=False,
                         cbar_num_of_ticks=7, cbar_round_decimals=0,
                         grid=False, landmask=True, **kwargs):
        ''' Save a raster band to a figure in graphical format.

        Parameters
        -----------
        fileName : string, optional
            Output file name. if one of extensions 'png', 'emf', 'eps', 'pdf',
            'rgba', 'ps', 'raw', 'svg', 'svgz' is included,
            specified file is crated. otherwise, 'png' file is created.
            if None, the Nansatmap object is returned.
            if True, the Nansatmap is shown
        contour : numpy 2D array, int or string
            input data, band number or band name
        contourf : numpy 2D array, int or string
            input data, band number or band name
        quiver : list of numpy 2D array, int or string
            (e.g. [array, array], [int, str], etc...)
            input data, band number or band name
        mesh : numpy 2D array, int or string
            input data, band number or band name
        color_bar : bool
            add color bar?
        cbar_num_of_ticks : int
            number of ticks on the colorbar
        cbar_round_decimals : int
            decimals of scale on the colorbar
        grid : bool
            draw grid?
        landmask : bool
            draw continents?
        **kwargs : parameters for nansatmap().
            See nansatmap.py

        Modifies
        ---------
        if fileName is specified, creates nansatmap file

        Returns
        -------
        Nansatmap object

        Example
        --------
        # write contour line and save the image
        n.write_nansatmap('test.jpg', contour=n[3])
        # put colors and write quiverplots
        n.write_nansatmap('test.jpg', mesh=1, quiver=['east_wind','north_wind'],
                          grid=True, color_bar=True)

        See also
        --------
        Nansatmap()
        http://matplotlib.org/basemap/api/basemap_api.html#mpl_toolkits.basemap.Basemap

        '''
        # if data is given by band number or name, get the array
        for dataVar in ['contour', 'contourf', 'mesh', 'quiver']:
            if locals()[dataVar] is not None:
                if type(locals()[dataVar])==list:
                    if type(locals()[dataVar][0]) != np.ndarray:
                        listData = []
                        for i, iValue in enumerate(locals()[dataVar]):
                            if type(locals()[dataVar][i])==str:
                                bandNum = self._get_band_number(locals()[dataVar])
                            if type(locals()[dataVar][i])==int:
                                bandNum = locals()[dataVar][i]
                            listData.append(self[bandNum])
                        if listData != []:
                            globals()[dataVar] = listData

                elif type(locals()[dataVar]) == np.ndarray and dataVar+'Valid' in kwargs:
                        globals()[dataVar+'Valid'] = kwargs.pop(dataVar+'Valid')

                elif type(locals()[dataVar]) != np.ndarray:
                    if type(locals()[dataVar])==str:
                        bandNum = self._get_band_number(locals()[dataVar])
                    elif type(locals()[dataVar])==int:
                        bandNum = locals()[dataVar]
                    # set Variable
                    globals()[dataVar] = self[bandNum]
                    # set varlid min and max values for colorbar ticks
                    if dataVar+'Valid' in kwargs:
                        globals()[dataVar+'Valid'] = kwargs.pop(dataVar+'Valid')
                    else:
                        try:
                            globals()[dataVar+'Valid'] = [
                            float(self.vrt.dataset.GetRasterBand(bandNum).GetMetadataItem('valid_min')),
                            float(self.vrt.dataset.GetRasterBand(bandNum).GetMetadataItem('valid_max'))]
                        except:
                            globals()[dataVar+'Valid'] = None
                else:
                    globals()[dataVar+'Valid'] = None

        # Create Nansatmap object
        argKeys = ['lcrnrlon', 'llcrnrlat', 'urcrnrlon', 'urcrnrlat',
                   'llcrnrx', 'llcrnry', 'urcrnrx', 'urcrnry',
                   'width', 'height', 'projection', 'resolution',
                   'area_thresh', 'rsphere', 'lat_ts',
                   'lat_0', 'lat_1', 'lat_2', 'lon_0', 'lon_1', 'lon_2',
                   'k_0', 'no_rot', 'suppress_ticks', 'satellite_height',
                   'boundinglat', 'fix_aspect', 'anchor', 'celestial',
                   'round', 'ax', 'num', 'figsize', 'dpi',
                   'facecolor', 'edgecolor', 'frameon']
        kwargs1 = self._pickup_args(kwargs, argKeys)
        nMap = Nansatmap(self, **kwargs1)

        # draw filled contour plot
        if contourf is not None:
            argKeys = ['smooth', 'mode', 'colors', 'alpha', 'cmap', 'norm',
                       'vmin', 'vmax', 'levels', 'origin', 'extent',
                       'locator', 'extend', 'xunits', 'yunits', 'antialiased',
                       'nchunk', 'hatches']
            kwargs1 = self._pickup_args(kwargs, argKeys)
            if type(contourf) != np.ndarray:
                contourf = globals()['contourf']

            nMap.contourf(contourf, contourfValid, cbar_num_of_ticks,
                          cbar_round_decimals, **kwargs1)

        # draw black smooth contour plot with labels
        if contour is not None:
            argKeys = ['smooth','contourFontsize','contourColors', 'alpha',
                       'cmap', 'norm', 'vmin', 'vmax', 'levels', 'origin',
                       'extent', 'locator', 'extend', 'xunits', 'yunits',
                       'antialiased', 'linewidths', 'linestyles']

            kwargs1 = self._pickup_args(kwargs, argKeys)

            if 'contourFontsize' in kwargs1.keys():
                kwargs1['fontsize'] = kwargs1.pop('contourFontsize')
            if 'contourColors' in kwargs1.keys():
                kwargs1['colors'] = kwargs1.pop('contourColors')
            if type(contour) != np.ndarray:
                contour = globals()['contour']
            nMap.contour(contour, contourValid, cbar_num_of_ticks,
                         cbar_round_decimals, **kwargs1)

        # pseudo-color plot over the map
        if mesh is not None:
            argKeys = ['cmap', 'norm', 'vmin', 'vmax', 'shading',
                       'edgecolors', 'alpha', 'agg_filter', 'animated',
                       'antialiased', 'array', 'axes', 'clim', 'clip_box',
                       'clip_on', 'clip_path', 'cmap', 'meshColor', 'colorbar',
                       'contains', 'edgecolor', 'facecolor', 'figure', 'gid',
                       'hatch', 'label', 'linestyle', 'linewidth', 'lod',
                       'norm', 'offset_position', 'offsets', 'paths',
                       'picker', 'pickradius', 'rasterized', 'snap',
                       'transform', 'url', 'urls', 'visible', 'zorder']
            kwargs1 = self._pickup_args(kwargs, argKeys)
            if 'meshColor' in kwargs1.keys():
                kwargs1['color'] = kwargs1.pop('meshColor')
            if type(mesh) != np.ndarray:
                mesh = globals()['mesh']
            nMap.pcolormesh(mesh, meshValid, **kwargs1)

        # quiver plot
        if quiver is not None:
            if type(quiver)==list and len(quiver)==2:
                argKeys = ['quivectors']
                kwargs1 = self._pickup_args(kwargs, argKeys)
                if type(quiver[0]) != np.ndarray:
                    quiver = globals()['quiver']
                nMap.quiver(quiver[0], quiver[1], **kwargs1)
            else:
                self.logger.warning('"quiver" mast be a list of two numpy arrays.')

        # add colorbar
        if color_bar:
            argKeys = ['orientation', 'pad', 'cbarFontsize']
            kwargs1 = self._pickup_args(kwargs, argKeys)
            if 'cbarFontsize' in kwargs1.keys():
                kwargs1['fontsize'] = kwargs1.pop('cbarFontsize')
            nMap.add_colorbar(**kwargs1)

        # add geocoordinates
        if grid:
            argKeys = ['gridFontsize', 'lat_num', 'lon_num',
                       'lat_labels', 'lon_labels']
            kwargs1 = self._pickup_args(kwargs, argKeys)
            if 'gridFontsize' in kwargs1.keys():
                kwargs1['fontsize'] = kwargs1.pop('gridFontsize')
            nMap.drawgrid(**kwargs1)

        # Save to a image file or Show
        if fileName is not None:
            argKeys = ['color', 'lake_color', 'ax', 'zorder', 'alpha']
            kwargs1 = self._pickup_args(kwargs, argKeys)
            if type(fileName)==bool and fileName:
                if landmask:
                    nMap.draw_continents(**kwargs1)
                plt.show()
            elif type(fileName)==str:
                nMap.save(fileName, landmask, **kwargs1)
        return nMap

    def _pickup_args(self, allkwargs, keys):
        kwargs = {}
        for iArg in keys:
            if iArg in allkwargs.keys():
                kwargs[iArg] = allkwargs[iArg]
        return kwargs

    def get_time(self, bandID=None):
        ''' Get time for dataset and/or its bands

        Parameters
        ----------
        bandID : int or str (default = None)
                band number or name

        Returns
        --------
        time : list with datetime objects for each band.
            If time is the same for all bands, the list contains 1 item

        '''
        time = []
        for i in range(self.vrt.dataset.RasterCount):
            band = self.get_GDALRasterBand(i + 1)
            try:
                time.append(dateutil.parser.parse(
                            band.GetMetadataItem('time')))
            except:
                self.logger.debug('Band ' + str(i + 1) + ' has no time')
                time.append(None)

        if bandID is not None:
            bandNumber = self._get_band_number(bandID)
            return time[bandNumber - 1]
        else:
            return time

    def get_metadata(self, key=None, bandID=None):
        ''' Get metadata from self.vrt.dataset

        Parameters
        ----------
        key : string, optional
            name of the metadata key. If not givem all metadata is returned
        bandID : int or str, optional
            number or name of band to get metadata from.
            If not given, global metadata is returned

        Returns
        --------
        a string with metadata if key is given and found
        an empty string if key is given and not found
        a dictionary with all metadata if key is not given

        '''
        # get all metadata from dataset or from band
        if bandID is None:
            metadata = self.vrt.dataset.GetMetadata()
        else:
            metadata = self.get_GDALRasterBand(bandID).GetMetadata()

        # get all metadata or from a key
        if key is not None:
            metadata = metadata.get(key, None)

        return metadata

    def set_metadata(self, key='', value='', bandID=None):
        ''' Set metadata to self.raw.dataset and self.vrt.dataset

        Parameters
        -----------
        key : string or dictionary with strings
            name of the metadata, or dictionary with metadata names, values
        value : string
            value of metadata
        bandID : int or str
            number or name of band
            Without : global metadata is set

        Modifies
        ---------
        self.raw.dataset : sets metadata in GDAL raw dataset
        self.vrt.dataset : sets metadata in GDAL current dataset

        '''
        # set all metadata to the dataset or to the band
        if bandID is None:
            metaReceiverRAW = self.raw.dataset
            metaReceiverVRT = self.vrt.dataset
        else:
            bandNumber = self._get_band_number(bandID)

            metaReceiverRAW = self.raw.dataset.GetRasterBand(bandNumber)
            metaReceiverVRT = self.vrt.dataset.GetRasterBand(bandNumber)

        # set metadata from dictionary or from single pair key,value
        if type(key) == dict:
            for k in key:
                metaReceiverRAW.SetMetadataItem(k, key[k])
                metaReceiverVRT.SetMetadataItem(k, key[k])
        else:
            metaReceiverRAW.SetMetadataItem(key, value)
            metaReceiverVRT.SetMetadataItem(key, value)

    def _get_mapper(self, mapperName, **kwargs):
        ''' Create VRT file in memory (VSI-file) with variable mapping

        If mapperName is given only this mapper will be used,
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
        mapperName : string, optional (e.g. 'ASAR' or 'merisL2')

        Returns
        --------
        tmpVRT : VRT object
            tmpVRT.dataset is a GDAL VRT dataset

        Raises
        --------
        Error : occurs if given mapper cannot open the input file

        '''
        # open GDAL dataset. It will be parsed to all mappers for testing
        try:
            gdalDataset = gdal.Open(self.fileName)
        except RuntimeError:
            print ('GDAL could not open ' + self.fileName +
                   ', trying to read with Nansat mappers...')
            gdalDataset = None
        if gdalDataset is not None:
            # get metadata from the GDAL dataset
            metadata = gdalDataset.GetMetadata()
        else:
            metadata = None

        tmpVRT = None

        if mapperName is not '':
            # If a specific mapper is requested, we test only this one.
            # Stripping off eventual 'mapper_' and '.py' and converting
            # to lowercase
            mapperName = mapperName.replace('mapper_',
                                            '').replace('.py', '').lower()
            try:
                mapper_module = __import__('mapper_' + mapperName)
            except ImportError:
                raise Error('Mapper ' + mapperName + ' not in PYTHONPATH')
            tmpVRT = mapper_module.Mapper(self.fileName, gdalDataset,
                                          metadata, **kwargs)
            self.mapper = mapperName
        else:
            # We test all mappers, import one by one
            for iMapper in self.mapperList:
                #get rid of .py extension
                iMapper = iMapper.replace('.py', '')
                self.logger.debug('Trying %s...' % iMapper)
                try:
                    mapper_module = __import__(iMapper)
                    #create a Mapper object and get VRT dataset from it
                    tmpVRT = mapper_module.Mapper(self.fileName, gdalDataset,
                                                  metadata, **kwargs)
                    self.logger.info('Mapper %s - success!' % iMapper)
                    self.mapper = iMapper
                    break
                except:
                    pass

        # if no mapper fits, make simple copy of the input DS into a VSI/VRT
        if tmpVRT is None and gdalDataset is not None:
            self.logger.warning('No mapper fits, returning GDAL bands!')
            tmpVRT = VRT(gdalDataset=gdalDataset)
            for iBand in range(gdalDataset.RasterCount):
                tmpVRT._create_band({'SourceFilename': self.fileName,
                                     'SourceBand': iBand + 1})
                tmpVRT.dataset.FlushCache()

        # if GDAL cannot open the file, and no mappers exist which can make VRT
        if tmpVRT is None and gdalDataset is None:
            raise GDALError('NANSAT can not open the file ' + self.fileName)

        return tmpVRT

    def _get_pixelValue(self, val, defVal):
        if val == '':
            return defVal
        else:
            return val

    def _get_band_number(self, bandID):
        '''Return absolute band number

        Check if given bandID is valid
        Return absolute number of the band in the VRT

        Parameters
        ----------
        bandID : int or str or dict
            if int : checks if such band exists and returns band_id
            if str : finds band with coresponding name
            if dict : finds first band with given metadata

        Returns
        --------
        int : absolute band number

        '''
        bandNumber = 0
        # if bandID is str: create simple dict with seraching criteria
        searchDict = None
        if type(bandID) == str:
            bandID = {'name': bandID}

        # if bandID is dict: search self.bands with seraching criteria
        if type(bandID) == dict:
            bandsMeta = self.bands()
            for b in bandsMeta:
                for key in bandID:
                    if (key in bandsMeta[b] and
                            bandID[key] == bandsMeta[b][key]):
                        bandNumber = b
                        break

        # if bandID is int and with bounds: return this number
        if (type(bandID) == int and
            bandID >= 1 and
            bandID <= self.vrt.dataset.RasterCount):
            bandNumber = bandID

        # if no bandNumber found - raise error
        if bandNumber == 0:
            raise OptionError('Cannot find band %s! '
                              'bandNumber is from 1 to %s'
                              % (str(bandID), self.vrt.dataset.RasterCount))

        return bandNumber

    def mosaic(self, files=[], bands=[], doReproject=True, maskName='mask',
               **kwargs):
        '''Mosaic input files. If images overlap, calculate average

        Convert all input files into Nansat objects, reproject onto the
        Domain of the current object, get bands, from each object,
        calculate average and STD, add averaged bands (and STD) to the current
        object.

        mosaic() tries to get band 'mask' from the input files. The mask
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

        mosaic() adds bands to the object, so it works only with empty, or
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
        nClass : child of Nansat, [Nansat]
            This class is used to read input files
        mapperName : str, ['']
            This mapper is used to read input files
        eResampleAlg : int, [0]
            agorithm for reprojection, see Nansat.reproject()

        '''
        # get Nansat child class for opening file
        nClass = kwargs.get('nClass', Nansat)

        # get mapper name for opening file
        mapperName = kwargs.get('mapperName', '')

        # get resampling method for reproject
        eResampleAlg = kwargs.get('eResampleAlg', 0)

        # get desired shape
        dstShape = self.shape()
        self.logger.debug('dstShape: %s' % str(dstShape))

        # preallocate 2D matrices for sum, sum of squares, count of products
        # and mask
        self.logger.debug('Allocating 2D matrices')
        avgMat = {}
        stdMat = {}
        for b in bands:
            avgMat[b] = np.zeros((dstShape[0], dstShape[1]))
            stdMat[b] = np.zeros((dstShape[0], dstShape[1]))

        cntMat = np.zeros((dstShape[0], dstShape[1]), 'float16')
        maskMat = np.zeros((2, dstShape[0], dstShape[1]), 'int8')

        # for all input files
        for i, f in enumerate(files):
            self.logger.info('Processing %s' % f)
            # open file using Nansat or its child class
            # the line below is for debugging
            #n = nClass(f, logLevel=self.logger.level, mapperName=mapperName)
            try:
                n = nClass(f, logLevel=self.logger.level,
                           mapperName=mapperName)
            except:
                self.logger.error('Unable to open %s' % f)
                continue

            # get metadata from the image (only last img metadata is kept)
            bandsMetadata = n.bands()

            # add mask band [0: nodata, 1: cloud, 2: land, 64: data]
            try:
                mask = n[maskName]
            except:
                self.logger.error('Cannot get mask from %s' % f)
                mask = 64 * np.ones(n.shape()).astype('int8')
                n.add_band(array=mask, parameters={'name': maskName})

            if doReproject:
                # reproject image and get reprojected mask
                try:
                    n.reproject(self, eResampleAlg=eResampleAlg)
                    mask = n[maskName]
                except:
                    self.logger.error('Unable to reproject %s' % f)
                    continue
            # if mask was not received from projected image
            # create zeros (out of swath) for blocking this image from
            # averaging
            if mask is None:
                self.logger.error('No mask in reprojected file %s!' % f)
                mask = np.zeros(n.shape()).astype('int8')

            # add data to counting matrix
            cntMatTmp = np.zeros((dstShape[0], dstShape[1]), 'float16')
            cntMatTmp[mask > 2] = 1
            cntMat += cntMatTmp
            # add data to mask matrix (maximum of 0, 1, 2, 64)
            maskMat[0, :, :] = mask
            maskMat[1, :, :] = maskMat.max(0)

            # add data to summation matrix
            for b in bands:
                self.logger.debug('    Adding %s to sum' % b)
                # get projected data from Nansat object
                a = None
                try:
                    a = n[b]
                except:
                    self.logger.error('%s is not in %s' % (b, n.fileName))
                if a is not None:
                    # mask invalid data
                    a[mask <= 2] = 0
                    # sum of valid values and squares
                    avgMat[b] += a
                    stdMat[b] += np.square(a)
            # destroy
            n = None

        # average products
        cntMat[cntMat == 0] = np.nan
        for b in bands:
            self.logger.debug('    Averaging %s' % b)
            # get average
            avg = avgMat[b] / cntMat
            # calculate STD
            # STD = sqrt(sum((x-M)^2)/n) = (sqrt((sum(x^2) -
            #                                2*mean(x)*sum(x) +
            #                                sum(mean(x)^2))/n))
            stdMat[b] = np.sqrt((stdMat[b] - 2.0 * avg * avgMat[b] +
                                np.square(avg) * cntMat) / cntMat)
            # set std
            avgMat[b] = avg

        # calculate mask (max of 0, 1, 2, 64)
        maskMat = maskMat.max(0)
        # if old 'valid' mask was applied in files, replace with new mask
        maskMat[maskMat == 128] = 64

        self.logger.debug('Adding bands')
        # add mask band
        self.logger.debug('    mask')
        self.add_band(array=maskMat, parameters={'name': maskName, 'long_name': 'L2-mask', 'standard_name': 'mask'})
        # add averaged bands with metadata
        for b in bands:
            self.logger.debug('    %s' % b)

            # get metadata of this band
            for bm in bandsMetadata:
                if bandsMetadata[bm]['name'] == b:
                    parameters = bandsMetadata[bm]

            self.add_band(array=avgMat[b], parameters=parameters)
            parameters['name'] = b + '_std'
            self.add_band(array=stdMat[b], parameters=parameters)

    def process(self, opts=None):
        '''Default L2 processing of Nansat object. Overloaded in childs.'''

    def export_band(self, fileName, bandID=1, driver='netCDF'):
        '''Export only one band of the Nansat object
        Get array from the required band
        Create temporary Nansat from the array
        Export temporary Nansat to file

        Parameters
        ----------
        fileName : str
            name of the output file
        bandID : int or str, [1]
            number of name of the band
        driver : str, ['netCDF']
            name of the GDAL Driver (format) to use

        '''
        # get array from self
        bandArray = self[bandID]
        # get root, band metadata
        rootMetadata = self.get_metadata()
        bandMetadata = self.get_metadata(bandID=bandID)
        # create temporary nansat
        tmpNansat = Nansat(domain=self, array=bandArray)
        # set metadata
        tmpNansat.set_metadata(rootMetadata)
        tmpNansat.set_metadata(bandMetadata, bandID=1)
        # export
        tmpNansat.export(fileName, driver=driver)

    def get_transect(self, points=None, bandList=[1], latlon=True, transect=True, returnOGR=False, layerNum=0):
        '''Get transect from two poins and retun the values by numpy array

        Parameters
        ----------
        points : tuple with one or more points or shape file name
            i.e. ((lon1, lat1),(lon2, lat2),(lon3, lat3), ...) or
                 ((col1, row1),(col2, row2),(col3, row3), ...)
        bandID : list of int or string
            elements of the list are band number or band Name
        latlon : bool
            If the points in lat/lon, then True.
            If the points in pixel/line, then False.
        transect : bool
            If True, get all transact values
            If False, get values of points
        returnOGR: bool
            If True, then return numpy array
            If False, return OGR object
        layerNum: int
            If shapefile is given as points, it is the number of the layer

        Returns
        --------
        if returnOGR:
            transect : OGR object with points coordinates and values
        else:
            transect : list or 
                values of the transect or OGR object with the transect values
            [lonVector, latVector] : list with longitudes, latitudes
            pixlinCoord : numpy array with pixels and lines coordinates

        '''
        # if shapefile is given, get corner points from it
        if type(points) == str:
            nansatOGR = Nansatshape(fileName=points)
            #nansatOGR = NansatOGR(fileName=points, layerNum=layerNum)
            points, latlon = nansatOGR.get_corner_points(latlon)

        # if points is not given, get points from GUI ...
        if points is None:
            firstBand =bandList[0]
            if type(firstBand) == str:
                firstBand = self._get_band_number(firstBand)
            browser = PointBrowser(self[firstBand])
            browser.get_points()
            points = tuple(browser.coordinates)
            latlon = False

        # get wkt
        wkt = self._get_projection(self.vrt.dataset)

        pixlinCoord = np.array([[],[]])
        for iPoint in range(len(points)):
            # if one point is given
            if type(points[iPoint]) != tuple:
                point0 = (points[0], points[1])
                point1 = (points[0], points[1])
            # if we want points ...
            elif not transect:
                point0 = (points[iPoint][0], points[iPoint][1])
                point1 = (points[iPoint][0], points[iPoint][1])
            # if we want a transect ...
            else:
                try:
                    point0 = points[iPoint]
                    point1 = points[iPoint+1]
                except:
                    break
            # if points in degree, convert them into pix/lin
            if latlon:
                pix, lin = self._transform_points([point0[0], point1[0]],[point0[1], point1[1]], DstToSrc=1)
                point0 = (pix[0], lin[0])
                point1 = (pix[1], lin[1])
            # compute Euclidean distance between point0 and point1
            length = int(np.hypot(point0[0]-point1[0], point0[1]-point1[1]))
            # if a point is given
            if length == 0:
                length = 1
            # get sequential coordinates on pix/lin between two points
            pixVector = list(np.linspace(point0[0], point1[0], length).astype(int))
            linVector = list(np.linspace(point0[1], point1[1], length).astype(int))
            pixlinCoord = np.append(pixlinCoord, [pixVector, linVector], axis=1)

        # convert pix/lin into lon/lat
        lonVector, latVector = self._transform_points(pixlinCoord[0], pixlinCoord[1], DstToSrc=0)

        transect = []
        # get data
        for iBand in bandList:
            if type(iBand) == str:
                iBand = self._get_band_number(iBand)
            data = self[iBand]
            # extract values
            transect.append(data[list(pixlinCoord[1]), list(pixlinCoord[0])].tolist())
        if returnOGR:
            NansatOGR = Nansatshape(wkt=wkt)
            NansatOGR.set_layer(lonlatCoord=[lonVector, latVector], pixlinCoord=pixlinCoord, fieldNames=map(str, bandList), fieldValues=transect)
            return NansatOGR
        else:
            return transect, [lonVector, latVector], pixlinCoord

