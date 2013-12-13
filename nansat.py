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
import scipy

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
        kwargs : additional arguments for mappers

        Creates
        --------
        self.mapperList : list of file names
            list of available working mappers
        self.fileName : file name
            set file name given by the argument
        self.vrt : VRT object
            Wrapper around VRT file and GDAL dataset with satellite raster data
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

        # create self.vrt from a file using mapper or...
        if fileName != '':
            # Make original VRT object with mapping of variables
            self.vrt = self._get_mapper(mapperName, **kwargs)
        # ...create using array, domain, and parameters
        else:
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

        # Set invalid and missing data to np.nan
        if '_FillValue' in band.GetMetadata():
            fillValue = float(band.GetMetadata()['_FillValue'])
            try:
                bandData[bandData == fillValue] = np.nan
            except:
                self.logger.info('Cannot replace _FillValue values with np.NAN!')
        try:
            bandData[np.isinf(bandData)] = np.nan
        except:
            self.logger.info('Cannot replace inf values with np.NAN!')

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

    def add_band(self, array, parameters=None, nomem=False):
        '''Add band from the array to self.vrt

        Create VRT object which contains VRT and RAW binary file and append it
        to self.vrt.subVRTs

        Parameters
        -----------
        array : Numpy array with band data
        parameters : dictionary, band metadata: wkv, name, etc.
        nomem : boolean, saves the vrt to a tempfile if nomem is True

        Modifies
        ---------
        Creates VRT object with VRT-file and RAW-file
        Adds band to the self.vrt

        '''
        # None => {} in input p
        if parameters is None:
            parameters = {}

        # create VRT from array
        bandVRT = VRT(array=array, nomem=nomem)

        self.vrt = self.vrt.get_super_vrt()

        # add the array band into self.vrt and get bandName
        bandName = self.vrt._create_band({'SourceFilename': bandVRT.fileName,
                                          'SourceBand': 1}, parameters)
        self.vrt.subVRTs[bandName] = bandVRT
        self.vrt.dataset.FlushCache()  # required after adding bands

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

    def has_band(self, band):
        '''Check if self has band with name <band>
        Parameters
        ----------
            band : str
                name of the band to check

        Returns
        -------
            True/False if band exists or not
                
        '''
        bandExists = False
        for b in self.bands():
            if self.bands()[b]['name'] == band:
                bandExists = True

        return bandExists

    def export(self, fileName, rmMetadata=[], addGeolocArray=True,
               addGCPs=True, driver='netCDF', bottomup=False):
        '''Export Nansat object into netCDF or GTiff file

        Parameters
        -----------
        fileName : str
            output file name
        rmMetadata : list
            metadata names for removal before export.
            e.g. ['name', 'colormap', 'source', 'sourceBands']
        addGeolocArray : bool
            add geolocation array datasets to exported file?
        addGCPs : bool
            add GCPs?  to exported file?
        driver : str
            Name of GDAL driver (format)
        bottomup : bool
            Write swath-projected data bottomup?

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
                bandMetadataR['name'] = bandMetadataR.pop('name') + '_real'
                bandMetadataI['name'] = bandMetadataI.pop('name') + '_imag'
                # Create bands from the real and imaginary numbers
                exportVRT.real.append(VRT(array=self[i].real))
                exportVRT.imag.append(VRT(array=self[i].imag))

                metaDict = [{'src': {
                             'SourceFilename': exportVRT.real[-1].fileName,
                             'SourceBand':  1},
                             'dst': bandMetadataR},
                            {'src': {
                             'SourceFilename': exportVRT.imag[-1].fileName,
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

        gcps = exportVRT.dataset.GetGCPs()
        if addGCPs and len(gcps) > 0:
            # add GCPs in VRT metadata and remove geotransform
            exportVRT._add_gcp_metadata(bottomup)
            exportVRT._remove_geotransform()

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

        # if output filename is same as input one...
        if self.fileName == fileName:
            numOfBands = self.vrt.dataset.RasterCount
            # create VRT from each band and add it
            for iBand in range(numOfBands):
                vrt = VRT(array=self[iBand + 1])
                self.add_band(vrt=vrt)
                metadata = self.get_metadata(bandID=iBand + 1)
                self.set_metadata(key=metadata,
                                  bandID=numOfBands + iBand + 1)

            # remove source bands
            self.vrt.delete_bands(range(1, numOfBands))

        # set CreateCopy() options
        if bottomup:
            options = 'WRITE_BOTTOMUP=NO'
        else:
            options = 'WRITE_BOTTOMUP=YES'

        # Create an output file using GDAL
        self.logger.debug('Exporting to %s using %s...' % (fileName, driver))
        dataset = gdal.GetDriverByName(driver).CreateCopy(fileName,
                                                          exportVRT.dataset,
                                                          options=[options])
        self.logger.debug('Export - OK!')

    def resize(self, factor=1, width=None, height=None, eResampleAlg=-1):
        '''Proportional resize of the dataset.

        The dataset is resized as (xSize*factor, ySize*factor) or
        (width, calulated height) or (calculated width, height).
        self.vrt is rewritten to the the downscaled sizes.
        Georeference is stored in the object. Useful e.g. for export.
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
                0 : NearestNeighbour
                1 : Bilinear,
                2 : Cubic,
                3 : CubicSpline,
                4 : Lancoz

        Modifies
        ---------
        self.vrt.dataset : VRT dataset of VRT object
            raster size are modified to downscaled size.
            If GCPs are given in the dataset, they are also overwritten.

        '''
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

        if eResampleAlg <= 0:
            self.vrt = self.vrt.get_subsampled_vrt(newRasterXSize,
                                                   newRasterYSize,
                                                   factor,
                                                   eResampleAlg)
        else:
            # update size and GeoTranform in XML of the warped VRT object
            self.vrt = self.vrt.get_resized_vrt(newRasterXSize,
                                                newRasterYSize,
                                                eResampleAlg=eResampleAlg)

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
            geoTransform = map(float, geoTransform)
            self.vrt.dataset.SetGeoTransform(geoTransform)

        # set global metadata
        subMetaData = self.vrt.vrt.dataset.GetMetadata()
        subMetaData.pop('fileName')
        self.set_metadata(subMetaData)
        
        return factor

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
                  WorkingDataType=None, tps=False, **kwargs):
        ''' Change projection of the object based on the given Domain

        Create superVRT from self.vrt with AutoCreateWarpedVRT() using
        projection from the dstDomain.
        Modify XML content of the warped vrt using the Domain parameters.
        Generate warpedVRT and replace self.vrt with warpedVRT.
        If current object spans from 0 to 360 and dstDomain is west of 0,
        the object is shifted by 180 westwards.

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
        blockSize : int
            size of blocks for resampling. Large value decrease speed
            but increase accuracy at the edge
        WorkingDataType : int (GDT_int, ...)
            type of data in bands. Shuold be integer for int32 bands
        tps : boolean
            Use thin-spline trasnformation or not

        Modifies
        ---------
        self.vrt : VRT object with dataset replaced to warpedVRT dataset

        See Also
        ---------
        http://www.gdal.org/gdalwarp.html
        '''
        # if no domain: quit
        if dstDomain is None:
            return

        # if self spans from 0 to 360 and dstDomain is west of 0:
        #     shift self westwards by 180 degrees
        # check span
        srcCorners = self.get_corners()
        if round(min(srcCorners[0])) == 0 and round(max(srcCorners[0])) == 360:
            # check intersection of src and dst
            dstCorners = dstDomain.get_corners()
            if min(dstCorners[0]) < 0:
                # shift
                self.vrt = self.vrt.get_shifted_vrt(-180)

        # get projection of destination dataset
        dstSRS = dstDomain.vrt.dataset.GetProjection()

        # get destination GCPs
        dstGCPs = dstDomain.vrt.dataset.GetGCPs()
        if len(dstGCPs) > 0:
            # get projection of destination GCPs
            dstSRS = dstDomain.vrt.dataset.GetGCPProjection()

        xSize = dstDomain.vrt.dataset.RasterXSize
        ySize = dstDomain.vrt.dataset.RasterYSize

        # get geoTransform
        if 'use_gcps' in kwargs.keys() and not (kwargs['use_gcps']):
            corners = dstDomain.get_corners()
            ext = '-lle %0.3f %0.3f %0.3f %0.3f -ts %d %d' % (min(corners[0]),
                                                              min(corners[1]),
                                                              max(corners[0]),
                                                              max(corners[1]),
                                                              xSize, ySize)
            d = Domain(srs=dstSRS, ext=ext)
            geoTransform = d.vrt.dataset.GetGeoTransform()
        else:
            geoTransform = dstDomain.vrt.dataset.GetGeoTransform()

        # create Warped VRT
        self.vrt = self.vrt.get_warped_vrt(dstSRS=dstSRS,
                                            dstGCPs=dstGCPs,
                                            eResampleAlg=eResampleAlg,
                                            xSize=xSize, ySize=ySize,
                                            blockSize=blockSize,
                                            geoTransform=geoTransform,
                                            WorkingDataType=WorkingDataType,
                                            tps=tps, **kwargs)

        # set global metadata from subVRT
        subMetaData = self.vrt.vrt.dataset.GetMetadata()
        subMetaData.pop('fileName')
        self.set_metadata(subMetaData)

    def undo(self, steps=1):
        '''Undo reproject, resize, add_band or crop of Nansat object

        Restore the self.vrt from self.vrt.vrt

        Parameters
        -----------
        steps : int
            How many steps back to undo

        Modifies
        --------
        self.vrt

        '''

        self.vrt = self.vrt.get_sub_vrt(steps)

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
            None (default) : min/max values are fetched from WKV,
            fallback-'hist'
            [min, max] : min and max are numbers, or
            [[min, min, min], [max, max, max]]: three bands used
            'hist' : a histogram is used to calculate min and max values
        addDate : boolean
            False (default) : no date will be aded to the caption
            True : the first time of the object will be added to the caption
        **kwargs : parameters for Figure().

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
            caption += self.get_time()[0].strftime(' %Y-%m-%d')

        self.logger.info('caption: %s ' % caption)

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
        ''' Set metadata to self.vrt.dataset

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
        self.vrt.dataset : sets metadata in GDAL current dataset

        '''
        # set all metadata to the dataset or to the band
        if bandID is None:
            metaReceiverVRT = self.vrt.dataset
        else:
            bandNumber = self._get_band_number(bandID)
            metaReceiverVRT = self.vrt.dataset.GetRasterBand(bandNumber)

        # set metadata from dictionary or from single pair key,value
        if type(key) == dict:
            for k in key:
                metaReceiverVRT.SetMetadataItem(k, key[k])
        else:
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
            # create VRT
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
                # get rid of .py extension
                iMapper = iMapper.replace('.py', '')
                self.logger.debug('Trying %s...' % iMapper)
                try:
                    mapper_module = __import__(iMapper)
                    # create a Mapper object and get VRT dataset from it
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
        if type(bandID) == str:
            bandID = {'name': bandID}

        # if bandID is dict: search self.bands with seraching criteria
        if type(bandID) == dict:
            bandsMeta = self.bands()
            for b in bandsMeta:
                numCorrectKeys = 0
                for key in bandID:
                    if (key in bandsMeta[b] and
                            bandID[key] == bandsMeta[b][key]):
                        numCorrectKeys = numCorrectKeys + 1
                    if numCorrectKeys == len(bandID):
                        bandNumber = b
                        break

        # if bandID is int and with bounds: return this number
        if (type(bandID) == int and bandID >= 1 and
                bandID <= self.vrt.dataset.RasterCount):
            bandNumber = bandID

        # if no bandNumber found - raise error
        if bandNumber == 0:
            raise OptionError('Cannot find band %s! '
                              'bandNumber is from 1 to %s'
                              % (str(bandID), self.vrt.dataset.RasterCount))

        return bandNumber

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

    def get_transect(self, points=None, bandList=[1], latlon=True,
                     transect=True, returnOGR=False, layerNum=0,
                     smooth=0, **kwargs):
        '''Get transect from two poins and retun the values by numpy array

        Parameters
        ----------
        points : tuple with one or more points or shape file name
            i.e. ((lon1, lat1),(lon2, lat2),(lon3, lat3), ...) or
                 ((col1, row1),(col2, row2),(col3, row3), ...)
        bandList : list of int or string
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
        smooth: int or [int, int]
            If smooth or smooth[0] is greater than 0, smooth every transect
            pixel as the median (default, smooth[1]=0) or mean (smooth[1]=1)
            value in a box with sides equal to the given number.
            smooth or smooth[0] must be 0 or a positive odd number
            smooth[1] can be 0 or 1 for median or mean

        vmin, vmax : int (optional)
            minimum and maximum pixel values of an image shown
            in case points is None.

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
        smooth_function = scipy.stats.nanmedian
        if type(smooth) is list:
            if smooth[0] < 0:
                raise ValueError("smooth[0] must be 0 or a positive odd number.")
            if smooth[0] > 0 and (smooth[0] % 2) != 1:
                raise ValueError("The kernel size smooth[0] should be odd.")
            if not smooth[1] == 0 and not smooth[1] == 1:
                raise ValueError("smooth[1] must be 0 or 1 for median or mean filter.")
            if smooth[1] == 1:
                smooth_function = scipy.stats.nanmean
        else:
            if smooth < 0:
                raise ValueError("smooth must be 0 or a positive odd number.")
            if smooth > 0 and (smooth % 2) != 1:
                raise ValueError("The kernel size smooth should be odd.")
            tmp = smooth
            smooth = []
            smooth.append(tmp)

        data = None
        # if shapefile is given, get corner points from it
        if type(points) == str:
            nansatOGR = Nansatshape(fileName=points)
            #nansatOGR = NansatOGR(fileName=points, layerNum=layerNum)
            points, latlon = nansatOGR.get_corner_points(latlon)

        # if points is not given, get points from GUI ...
        if points is None:
            firstBand = bandList[0]
            if type(firstBand) == str:
                firstBand = self._get_band_number(firstBand)
            data = self[firstBand]

            browser = PointBrowser(data, **kwargs)
            browser.get_points()
            points = tuple(browser.coordinates)
            latlon = False

        # get wkt
        wkt = self.vrt.get_projection()

        pixlinCoord = np.array([[], []])
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
                    point1 = points[iPoint + 1]
                except:
                    break
            # if points in degree, convert them into pix/lin
            if latlon:
                pix, lin = self.transform_points([point0[0], point1[0]],
                                                  [point0[1], point1[1]],
                                                  DstToSrc=1)
                point0 = (pix[0], lin[0])
                point1 = (pix[1], lin[1])
            # compute Euclidean distance between point0 and point1
            length = int(np.hypot(point0[0] - point1[0],
                                  point0[1] - point1[1]))
            # if a point is given
            if length == 0:
                length = 1
            # get sequential coordinates on pix/lin between two points
            pixVector = list(np.linspace(point0[0],
                                         point1[0],
                                         length).astype(int))
            linVector = list(np.linspace(point0[1],
                                         point1[1],
                                         length).astype(int))
            pixlinCoord = np.append(pixlinCoord,
                                    [pixVector, linVector],
                                    axis=1)
            if smooth[0]:
                pixlinCoord0 = pixlinCoord - int(smooth[0]) / 2
                pixlinCoord1 = pixlinCoord + int(smooth[0]) / 2

        # convert pix/lin into lon/lat
        lonVector, latVector = self.transform_points(pixlinCoord[0],
                                                      pixlinCoord[1],
                                                      DstToSrc=0)
        transect = []
        # get data
        for iBand in bandList:
            if type(iBand) == str:
                iBand = self._get_band_number(iBand)
            if data is None:
                data = self[iBand]
            # extract values
            if smooth[0]:
                transect0 = []
                for xmin, xmax, ymin, ymax in zip(pixlinCoord0[1],
                                                  pixlinCoord1[1],
                                                  pixlinCoord0[0],
                                                  pixlinCoord1[0]):
                    transect0.append(smooth_function(data[xmin:xmax,
                                                          ymin:ymax],
                                                     axis=None))
                transect.append(transect0)
            else:
                transect.append(data[list(pixlinCoord[1]),
                                list(pixlinCoord[0])].tolist())
            data = None
        if returnOGR:
            # Lists for field names and datatype
            names = ['X (pixel)', 'Y (line)']
            formats = ['i4', 'i4']
            for iBand in bandList:
                names.append('transect_' + str(iBand))
                formats.append('f8')
            # Create zeros structured numpy array
            fieldValues = np.zeros(len(pixlinCoord[1]),
                                   dtype={'names': names,
                                          'formats': formats})
            # Set values into the structured numpy array
            fieldValues['X (pixel)'] = pixlinCoord[0]
            fieldValues['Y (line)'] = pixlinCoord[1]
            for i, iBand in enumerate(bandList):
                fieldValues['transect_' + str(iBand)] = transect[i]
            # Create Nansatshape object
            srs = osr.SpatialReference()
            srs.ImportFromWkt(wkt)
            NansatOGR = Nansatshape(srs=srs)
            # Set features and geometries into the Nansatshape
            NansatOGR.add_features(coordinates=[lonVector, latVector],
                                   values=fieldValues,
                                   AddPixLine=False)
            # Return Nansatshape object
            return NansatOGR
        else:
            return transect, [lonVector, latVector], pixlinCoord.astype(int)

    def crop(self, xOff=0, yOff=0, xSize=None, ySize=None, lonlim=None, latlim=None):
        '''Crop Nansat object
        
        Create superVRT, modify the Source Rectangle (SrcRect) and Destination
        Rectangle (DstRect) tags in the VRT file for each band in order
        to take only part of the original image,
        create new GCPs or new GeoTransform for the cropped object.
        
        Parameters
        ----------
        xOff : int
            pixel offset of subimage
        yOff : int
            line offset of subimage
        xSize : int
            width in pixels of subimage
        ySize : int
            height in pizels of subimage
        
        Modifies:
            self.vrt : VRT
                superVRT is created with modified SrcRect and DstRect
        '''
        # use interactive PointBrowser for selecting extent
        if      (xOff==0 and yOff==0 and
                 xSize is None and ySize is None and
                 lonlim is None and latlim is None):
            factor = self.resize(width=1000)
            data = self[1]
            browser = PointBrowser(data)
            browser.get_points()
            points = np.array(browser.coordinates)
            xOff = round(points.min(axis=0)[0] / factor)
            yOff = round(points.min(axis=0)[1] / factor)
            xSize = round((points.max(axis=0)[0] - points.min(axis=0)[0]) / factor)
            ySize = round((points.max(axis=0)[1] - points.min(axis=0)[1]) / factor)
            print xOff, yOff, xSize, ySize
            self.undo()

        # get xOff, yOff, xSize and ySize from lonlim and latlim
        if       (xOff==0 and yOff==0 and
                 xSize is None and ySize is None and
                 type(lonlim) is list and type(latlim) is list):
            crnPix, crnLin = self.transform_points([lonlim[0], lonlim[0], lonlim[1], lonlim[1]],
                                                   [latlim[0], latlim[1], latlim[0], latlim[1]],
                                                    1)
            xOff = round(min(crnPix))
            yOff = round(min(crnLin))
            xSize = round(max(crnPix) - min(crnPix))
            ySize = round(max(crnLin) - min(crnLin))
        
        self.vrt.dataset.RasterYSize, self.vrt.dataset.RasterXSize
        
        RasterXSize = self.vrt.dataset.RasterXSize
        RasterYSize = self.vrt.dataset.RasterYSize
        
        if xSize is None:
            xSize =  RasterXSize - xOff
        if ySize is None:
            ySize =  RasterYSize - yOff
        
        self.logger.debug('xOff: %d, yOff: %d, xSize: %d, ySize: %d' % xOff,
                                                                       yOff,
                                                                       xSize,
                                                                       ySize)
        if    ((xSize + xOff) > RasterXSize or
               (ySize + yOff) > RasterYSize or
               xOff > RasterXSize or
               yOff > RasterXSize or
               xOff < 0 or xSize < 0 or
               yOff < 0 or ySize < 0):
            self.logger.error('WARNING! Cropping region is outside the image!')
            return
        
        self.vrt = self.vrt.get_super_vrt()
        xml = self.vrt.read_xml()
        node0 = Node.create(xml)
        
        node0.node('VRTDataset').replaceAttribute('rasterXSize', str(xSize))
        node0.node('VRTDataset').replaceAttribute('rasterYSize', str(ySize))
        
        # replace xOff and xSize in <SrcRect> and <DstRect> of each source
        for iNode1 in node0.nodeList('VRTRasterBand'):
            iNode2 = iNode1.node('ComplexSource')

            iNode3 = iNode2.node('SrcRect')
            iNode3.replaceAttribute('xOff', str(xOff))
            iNode3.replaceAttribute('yOff', str(yOff))
            iNode3.replaceAttribute('xSize', str(xSize))
            iNode3.replaceAttribute('ySize', str(ySize))

            iNode3 = iNode2.node('DstRect')
            iNode3.replaceAttribute('xSize', str(xSize))
            iNode3.replaceAttribute('ySize', str(ySize))
            

        xml = node0.rawxml()
        self.vrt.write_xml(xml)
        
        if xOff > 0 or yOff > 0:
            gcps = self.vrt.dataset.GetGCPs()
            if len(gcps) > 0:
                newGCPs = []
                gcpProjection = self.vrt.dataset.GetGCPProjection()
                for newPix in np.r_[0:xSize:10j]:
                    for newLin in np.r_[0:ySize:10j]:
                        newLon, newLat = self.vrt.transform_points([newPix+xOff], [newLin+yOff])
                        newGCPs.append(gdal.GCP(newLon[0], newLat[0], 0, newPix, newLin))
                self.vrt.dataset.SetGCPs(newGCPs, gcpProjection)
                self.vrt._remove_geotransform()
            else:
                geoTransfrom = self.vrt.dataset.GetGeoTransform()
                geoTransfrom = map(float, geoTransfrom)
                geoTransfrom[0] += geoTransfrom[1] * xOff
                geoTransfrom[3] += geoTransfrom[5] * yOff
                self.vrt.dataset.SetGeoTransform(geoTransfrom)

        # set global metadata
        subMetaData = self.vrt.vrt.dataset.GetMetadata()
        subMetaData.pop('fileName')
        self.set_metadata(subMetaData)
