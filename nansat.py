#-------------------------------------------------------------------------------
# Name:    nansat
# Purpose: main of nansat module.
#          Reference nansat_open, nansat_transform and nansat_write
#
# Author:      asumak
#
# Created:     29.06.2011
# Copyright:   (c) asumak 2011
# Licence:
#-------------------------------------------------------------------------------
import os
import os.path
import sys
import time
import fnmatch

try:
    from osgeo import gdal
except ImportError:
    import gdal

try:
    from osgeo import osr
except ImportError:
    import osr

from xml.etree.ElementTree import *
from scipy.misc import toimage
from scipy.misc import pilutil
import string
import re
import numpy as np
from scipy.stats import cumfreq

from domain import Domain
from vrt import *

class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass;

class GDALError(Error):
    '''Error from GDAL '''
    pass;

class ProjectionError(Error):
    '''Cannot get the projection'''
    pass;

class DataError(Error):
    '''Error for data.
        e.g. : empty pixel value array in get_pixelValueRange()'''
    pass;

class OptionError(Error):
    '''Error for unproper options (arguments) '''
    pass;

class Nansat():
    '''Main of Nansat

    Construct Nansat object that consist of
        basic dataset information (file location, metadata etc..),
        well know variavles that are defined in NANSAT and
        VRT format which is saved in an XML format.
    Show information of the bands in a given object.
    Give GDALRasterBand of the given band number.
    Export in-memory VRT dataset to a physical file.

    '''

    def __init__(self, fileName, mapperName = '', bandList = None):
        '''Construct Nansat object

        Open GDAL dataset,
        Read metadata,
        Identify type of the sensor,
        Map the variables into internal format,
        Generate GDAL VRT file in memory

        Args:
            fileName: location of the file

        Side effects:
            set attributes: fileName, ds, metadata, rawVrt,
            warpedVrt and vrt

        '''
        # location of the data
        self.fileName = fileName;

        # dataset
        self.ds = gdal.Open(self.fileName);
        if (self.ds == None) or (self.ds == ""):
            raise GDALError("Nansat._init_(): Cannot get the dataset from " + self.fileName);

        # metadata
        self.metadata = self.ds.GetMetadata();
        if (self.metadata == None) or (self.metadata == ""):
            raise GDALError("Nansat._init_(): Cannot get the metdadata");

        # names of raw and warped VRT files in memory
        self.rawVRTName = '/vsimem/vsiFile.vrt';
        self.warpedVRTName = '/vsimem/vsi_warped.vrt';
        self.vrtDriver = gdal.GetDriverByName("VRT");

        # VRT with mapping of variables
        self.rawVrt = self._get_mapper(mapperName, bandList);
        # Warped VRT
        self.warpedVrt = None;
        # Current VRT
        self.vrt = self.rawVrt;

    def list_bands(self):
        '''Show band information of the given Nansat object

        Show serial number, longName, name and all parameters
        for each band in the metadata of the given Nansat object.

        Side effects:
            show serial number, longName, name and parameters

        '''
        print self.fileName;
        for i in range(self.rawVrt.RasterCount):
            metadata = self.rawVrt.GetRasterBand(i+1).GetMetadata();
            print "Band :", i+1;
            for j in metadata:
                if j != "units":
                    print "    ", j, " : ",\
                            self.rawVrt.GetRasterBand(i+1).GetMetadataItem(j);

    def downscale(self, factor = 1, method = "average"):
        '''Downscale the size of the data.

        The size of data is downscaled as (xSize/factor, ySize/factor).
        self.vrt is rewritten to the the downscaled sizes.

        Args:
            factor: int
            method: "average" (default) or
                    "subsample" ( = nearest neighbor)

        Side effects:
            edit self.vrt

        Raises:
            OptionError: occur when method is neither "average" nor "subsample".

        '''
        if not (method == "average" or method == "subsample"):
            raise OptionError("method should be 'average' or 'subsample'");

        # Write the vrt to a VSI-file
        vrtDsCopy = self.vrtDriver.CreateCopy(self.rawVRTName, self.vrt);

        # Get XML content from VSI-file
        # open
        vsiFile = gdal.VSIFOpenL(self.rawVRTName, "r")
        # get file size
        gdal.VSIFSeekL(vsiFile, 0, 2)
        vsiFileSize = gdal.VSIFTellL(vsiFile)
        gdal.VSIFSeekL(vsiFile, 0, 0) #fseek to start
        # read
        vsiFileContent = gdal.VSIFReadL(vsiFileSize, 1, vsiFile);
        gdal.VSIFCloseL(vsiFile);

        # Get element from the XML content and modify some elements
        # using Domain object parameters
        element = XML(vsiFileContent);
        XSize = int(float(element.get("rasterXSize")) / factor);
        YSize = int(float(element.get("rasterYSize")) / factor);
        element.set("rasterXSize", str(XSize));
        element.set("rasterYSize", str(YSize));

        for elem in element.iter("DstRect"):
            elem.set("xSize", str(XSize));
            elem.set("ySize", str(YSize));

        # if method = "average", overwrite "SimpleSource" to "AveragedSource"
        if method == "average":
            for elem in element.iter("SimpleSource"):
                elem.tag = "AveragedSource";

        # Edit GCPs to correspond to the downscaled size
        for elem in element.iter("GCP"):
            pxl = float(elem.get("Pixel")) / factor;
            if pxl > float(XSize):
                pxl = XSize;
            lin = float(elem.get("Line")) / factor;
            if lin > float(YSize):
                lin = YSize;
            elem.set("Pixel", str(pxl));
            elem.set("Line", str(lin));

        # Overwrite element
        # Write the modified elemements into VSI-file
        vsiFile = gdal.VSIFOpenL(self.rawVRTName, 'w')
        gdal.VSIFWriteL(tostring(element), \
                            len(tostring(element)), 1, vsiFile);
        gdal.VSIFCloseL(vsiFile);

        self.vrt = gdal.Open(self.rawVRTName);


    def get_GDALRasterBand(self, bandNo = 1, bandID = None):
        '''Get a GDALRasterBand of a given Nansat object.

        Get a GDALRasterBand that is specified by the given argument.

        If a bandID is given, secify a bandNo based on it.
        Check if the given bandNo is proper.
        Get a GDALRasterBand from vrt.

        Args:
            bandNo: a serial number of the band to fetch.
                    default setting is 1
            bandID: a dictionary for parameters to specify a band
                    (example: bandIdList = {"ShortName":"radiance",
                                             "Wavelength":"1240"})
            bandID is prior to bandNo

        Returns :
            a GDAL RasterBand

        Raises:
            OptionError: An error occur when the bandNo is not a proper number.

        '''
        # If bandID is given, bandNo is specified here.
        if bandID != None:
            bandNo = self._specify_bandNo(bandID);
        # if given bandNo is over the existing bands, give error message
        elif 1 > bandNo or bandNo > self.rawVrt.RasterCount :
           raise OptionError("Nansat.get_GDALRasterBand(): bandNo takes from 1 to",\
                            self.rawVrt.RasterCount);

        # Based on bandNo,
        # the GDAL RasterBand of the corresponding band is returned
        return self.vrt.GetRasterBand(bandNo);

    def write_figure(self, fileName, bandNo = 1, bandName = None,\
                     pixelValMin = None, pixelValMax = None,\
                     imageDatatype = None, thresholdRatio = 1.0, useFullMatrix = False, extension = 'png'):
        '''Get proper pixel value range for writing a figure in PNG

        Save a raster band to a figure in PNG format.
        If bandName is used as the argument,
        it is converted to the bandNo in get_GDALRasterBand().
        The proper pixel value range is calculated by setting thresholdRatio.

        Args:
            fileName : file name
            bandNo  : int
            bandName (option) : a list
                (e.g.: bandIdList = {"name":"radiance", "wavelength":"645"})
            thresholdRatio (option) : float (0.0 - 1.0).
                e.g. : thresholdRatio = 0.95 means to round off 5%
                        form the both sides (upper and lower sides).
            useFullMatrix (option): boolean
                if true, the full matrix is used for estimating min/max,
                otherwise only image scaled down to 100x100 (MUCH FASTER)

        Raises:
            DataError: occurs when the array of the band is empty


        Side effects:
            write a band in PNG file

        '''
        #fetch band from the object
        if bandName != None:
            band = self.get_GDALRasterBand(bandID = bandName);
        else:
            band = self.get_GDALRasterBand(bandNo);

        #read NumPy array from band
        tic = time.clock();
        print "Writing figure (%d x %d) " %  (band.XSize, band.YSize), ;
        rawArray = band.ReadAsArray();
        if rawArray == None:
            raise DataError("Nansat.write_figure(): array of the band is empty");
        toc = time.clock();
        print "(%3.1f sec) " % (toc-tic),

        # if value < pixelValMin then replace as value = pixelValMin
        # if value > pixelValMax then replace as value = pixelValMax
        if pixelValMin is None:

            #reduce input matrix to the size 100 x 100 for calculating histogram
            if not useFullMatrix:
                step1 = max(rawArray.shape[0] / 100, 1);
                step2 = max(rawArray.shape[1] / 100, 1);
                histArray = rawArray[::step1, ::step2];
            else:
                histArray = rawArray;

            #get minmax from histogram analysis
            pixelValMin, pixelValMax = self._get_pixelValueRange\
                                           (histArray, thresholdRatio);
        print "[%f %f]" % (pixelValMin, pixelValMax),
        toc = time.clock();
        print "(%3.1f sec) " % (toc-tic),

        #cut away values over limits and save to a PNG
        np.clip(rawArray, pixelValMin, pixelValMax, out=rawArray);
        toimage(rawArray).save(fileName + "." + extension);
        toc = time.clock();
        print "(%3.1f sec) " % (toc-tic)

    def _get_pixelValueRange(self, array, ratio):
        '''Get proper pixel value range for writing a figure in PNG

        Return a proper pixel value range (cmin, cmax)
        to wrige a figure with toimage.
        the argument "ratio" is used to specify the threshold of a pixel value
        that should be counted.

        Args:
            array : array of a band
            ratio  : float (0.0 - 1.0)

        Returns:
            edge_min : float
            edge_max : float
        '''
        #exclude zeros from array (wich spoil histo)
        array.flatten();
        array = array[array != 0];

        #try to make histogram
        tic = time.clock();
        try:
            hist, lowerreallimit, binsize, extrapoint = cumfreq(array, numbins = 15);
        except:
            hist = None
            
        if hist is None:
            edge_min = 0
            edge_max = 1
        else:
            toc = time.clock();
            #print "hist : ", hist;
            #print "lowerreallimit : ", lowerreallimit, "binsize : ", binsize;
            #print "array : ", np.histogram(array, bins=15);
    
            hist_eq = hist / max(hist);
            #print "hist_eq : ", hist_eq;
            hist_min = hist_eq[hist_eq < 1 - ratio];
            hist_max = hist_eq[hist_eq > ratio];
    
            if len(hist_min) == len(hist_eq):
                edge_min = lowerreallimit + (len(hist_eq) - 1.5) * binsize;
            elif len(hist_min) == 0:
                edge_min = lowerreallimit + 0.5 * binsize;
            else:
                edge_min = lowerreallimit + (len(hist_min) - 0.5) * binsize;
    
            if len(hist_max) == len(hist_eq):
                edge_max = lowerreallimit + (1.0 + 0.5) * binsize;
            elif len(hist_eq) == 0:
                edge_max = lowerreallimit + (len(hist_eq) - 0.5) * binsize;
            else:
                edge_max = lowerreallimit + (len(hist_eq) - len(hist_max) + 0.5) * binsize;

        return edge_min, edge_max;

    def export_VRT(self, fileName = None):
        '''Export in-memory VRT dataset to a physical file

        Check if fileName is given.
        If fileName is None, this method is skipped.
        Otherwise, open VSI-file.
        Copy it to a physical file whose location is given by the argument.

        Args:
            fileName: location for an output VRT file

        Side effects:
            "Create a physical VRT file and write down the data
             from memory (VSI-file).

        '''
        if fileName == None:
            fileName = "/vsistdout/" # GDAL special name to flush output to console.
            #Unfortunately an error message about non-existing file
            # is still reported; this must be a bug in GDAL.
        vrtDsCopy = self.vrtDriver.CreateCopy(fileName, self.vrt);

    def _get_mapper(self, mapperName, bandList):
        '''write VSI-file

        Loop over all availble drivers to get the corresponding one.
        All mappers are in separate '*.py' files in the subdir 'mappers'.
        All mappers are imported in the vrt.py file.
        In the loop:
            If the specific error appears the mapper is not used
            and the next mapper is tested.
            Otherwise the mapper returns VRT.

        Returns:
            vsiDs

        Side effects:
            Create VSI-file

        Raises:
            TypeError: occurs when the given driver type is not registarated
                        in the mappers.

        '''
        # create a mapper list based on the files in the folder "mappers"
        nansatdir = os.path.dirname(os.path.realpath( __file__ ))

        allMapperFiles = os.listdir(os.path.join(nansatdir, "mappers"))
        allMapperFiles = fnmatch.filter(allMapperFiles, 'mapper_*.py')

        #add the given mapper first
        mapperList = ['mapper_'+ mapperName];

        #loop through appropriate files and add to the list
        for ifile in allMapperFiles:
            ifile = ifile.replace(".py", "");
            mapperList.append(ifile);

        #try to add path for windows, add for linux otherwise
        try:
            sys.path.append(os.path.join(unicode(nansatdir, "mbcs"), "mappers"));
        except:
            sys.path.append(os.path.join(nansatdir, "mappers"));

        #try to import and get VRT datasaet from all mappers. Break on success
        #if none of the mappers worked - None is returned
        vrtDs = None;
        for iMapper in mapperList:
            try:
                mapper_module = __import__(iMapper);
                vrtDs = mapper_module.Mapper(self.ds, self.fileName, \
                                             self.metadata, \
                                             vrtBandList = bandList, \
                                             rawVRTName = self.rawVRTName).vsiDs;
                break;
            except:
                pass;

        #if no mapper fits, make simple copy of the input DS into a VSI/VRT
        if vrtDs is None:
            print 'No mapper fits!'
            vrtDs = self.vrtDriver.CreateCopy(self.rawVRTName, self.ds);

        return vrtDs;

    def dereproject(self):
        '''Cancel reprojection'''
        self.vrt = self.rawVrt;

    def _specify_bandNo(self, bandID):
        '''Specify a band number based on bandID (shortName + parameters)

        Check if the keys given by the argument(bandID)
          are in metadata keys.
        Compare the key values of the bandID
            to the values of the metadata dictionary.
        If they are matched, append the band number (iBand) into candidate list.
        If not, go to the next band.
        Iterate these steps until all bands are checked.
        If single band is specified at the end, return the band number.
        Otherwise show SpecificationError.

        Args:
            bandID: dictionary that shows parameters and the values
                    to specify single band.
                    (e.g.  {"ShortName":"radiance", "Wavelength":" 1234"})

        Returns:
            candidate[0]+1 : a band number

        Side effects:
            show the specified band number

        Raises:
            KeyError: occurs when the given key is not in the
                      rawVrt.GetRasterBand(1).GetMetadata_Dict().keys().

        '''
        metaItemKeys = self.rawVrt.GetRasterBand(1).GetMetadata_Dict().keys();
        bandIDkeys = bandID.keys();
        bandIDvalues = bandID.values();

        # check if the keys in the bandID exist
        for i in range(len(bandIDkeys)):
            if (bandIDkeys[i] not in metaItemKeys):
                raise KeyError("Nansat.specify_bandNo(): Cannot find a such key: ", \
                                          bandIDkeys[i]);

        # search for the specific band based on bandID
        candidate = [];
        for iBand in range(self.rawVrt.RasterCount):
            counter = 0;
            for iItemKey in bandIDkeys:
                counter += 1;
                if bandID[iItemKey] != self.rawVrt.GetRasterBand(iBand+1).\
                                            GetMetadataItem(iItemKey):
                    break;
                else:
                    if counter == len(bandIDkeys):
                        candidate.append(iBand);

        # if a band is specified, set it to bandNo.
        # if some bands are chosen, give an error message and stop.
        if len(candidate) == 1:
            print "You chose bandNo:", candidate[0]+1;
            return candidate[0]+1;
        elif len(candidate) >= 2:
            raise OptionError("Nansat._specify_bandNo(): Cannot specify a single band by the given arguments");
        else:
            raise OptionError("Nansat._specify_bandNo(): Cannot find any band by the given arguments");

    def reproject(self, proj4string = \
                            "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs",\
                             extentOption = None, trgDomain = None, \
                             resamplingAlg = 0):
        '''Reproject the object based on the given arguments proj4string
        and extentString

        Get srcWkt from the raw VRT
        Create Domain object from pro4String and extentString.
        Warp the raw VRT using AutoCreateWarpedVRT() using projection from the Domain
        Modify XML content of the warped vrt using the Domain parameters.
        Generate self.warpedVrt and replace self.vrt to warpedVrt.

        Args:
            proj4string : proj4string
            extentOption (option): string that shows extent.
                                   "te", "ts", "tr", "lle"
            trgDomain : target Domain

        Side effects:
            self.warpedVrt : warped vrt
            self.vrt : replace current self.vrt to self.warpedVrt

        Raises:
            ProjectionError: occurs when the projection of the source data is None.
            ProjectionError: occurs when the projection of the target data is None.
            OptionError: occures when the option combination is not proper.
            AttributeError: occurs when it is impossible to get warpedVRT.

        '''
        #Generate source WKT
        srcWKT = self.rawVrt.GetProjection()
        if srcWKT == None:
            raise ProjectionError("Nansat.reproject(): rawVrt.GetProjection() is None");

        #check input options
        if proj4string is not None and extentOption is None and trgDomain is None:
            #generate destination WKT
            dstSRS = osr.SpatialReference();
            dstSRS.ImportFromProj4(proj4string);
            dstWKT = dstSRS.ExportToWkt();
            if dstWKT == "":
                raise ProjectionError("Nansat.reproject(): Projection of the target data is empty." \
                                       "Is the 'proj4string' correct?");
            #warp the Raw Vrt onto the coordinate stystem given by proj4string
            rawWarpedVRT = gdal.AutoCreateWarpedVRT(self.rawVrt, \
                              srcWKT, dstWKT, resamplingAlg);

            #generate Domain from the warped VRT
            trgDomain = Domain(rawWarpedVRT);

        elif proj4string is not None and extentOption is not None and trgDomain is None:

            #Generate Domain from srs and extent strings
            trgDomain = Domain(ds = None, srsString = proj4string, extentString = extentOption);

            #warp the Raw Vrt onto the coordinate stystem given by Domain
            rawWarpedVRT = gdal.AutoCreateWarpedVRT(self.rawVrt, \
                              srcWKT, trgDomain.memDs.GetProjection(), resamplingAlg);
            
        elif proj4string is None and extentOption is None and trgDomain is not None:
            print 'Reprojection with given Domain'

            #warp the Raw Vrt onto the coordinate stystem given by input Domain
            rawWarpedVRT = gdal.AutoCreateWarpedVRT(self.rawVrt, \
                              srcWKT, trgDomain.memDs.GetProjection(), resamplingAlg);

        else:
            #Potentially erroneous input options
            raise OptionError("Nansat.reproject():wrong combination of input options");
        
        #modify extent of the created Warped VRT
        self.warpedVRT = self._modify_warped_vrt(rawWarpedVRT,
                                trgDomain.memDs.RasterXSize, \
                                trgDomain.memDs.RasterYSize, \
                                trgDomain.memDs.GetGeoTransform())

        #test created Warped VRT
        if self.warpedVRT == None:
            raise AttributeError("Nansat.reproject():cannot get warpedVRT");

        #set default vrt to be the warped one
        self.vrt = self.warpedVRT;
        
        
    def _modify_warped_vrt(self, rawWarpedVRT, rasterXSize, rasterYSize, geoTransform):
        ''' Modify rasterXsize, rasterYsize and geotranforms in the warped VRT
        Args:
            rasterXSize: integer, desired X size of warped image
            rasterYSize: integer, desired Y size of warped image
            rasterYSize: tuple of 6 integers, desired GeoTransform size of the warped image

        Side effects:
            the VRT file which keepes warped vrt is modified
        '''
        # Write the warpedVrt to a VSI-file
        vrtDsCopy = self.vrtDriver.CreateCopy(self.warpedVRTName, \
                                          rawWarpedVRT);
        # Get XML content from VSI-file
        #open
        vsiFile = gdal.VSIFOpenL(self.warpedVRTName, "r")
        # get file size
        gdal.VSIFSeekL(vsiFile, 0, 2)
        vsiFileSize = gdal.VSIFTellL(vsiFile)
        gdal.VSIFSeekL(vsiFile, 0, 0) #fseek to start
        # read
        vsiFileContent = gdal.VSIFReadL(vsiFileSize, 1, vsiFile);
        gdal.VSIFCloseL(vsiFile);

        # Get element from the XML content and modify some elements
        # using Domain object parameters
        #print vsiFileContent
        element = XML(vsiFileContent);
        element.set("rasterXSize", str(rasterXSize));
        element.set("rasterYSize", str(rasterYSize));
        tree = ElementTree(element);

        elem = tree.find("GeoTransform");
        #print "539 GeoTarmsform:", str(d.memDs.GetGeoTransform()).translate(string.maketrans("", ""), "()");
        # convert proper string style and set to the GeoTransform element
        elem.text = str(geoTransform).\
                        translate(string.maketrans("", ""), "()");

        elem = tree.find("GDALWarpOptions/Transformer/GenImgProjTransformer/DstGeoTransform");
        #print "545 DstGeoTransform:",str(d.memDs.GetGeoTransform()).translate(string.maketrans("", ""), "()");
        # convert proper string style and set to the DstGeoTransform element
        elem.text = str(geoTransform).\
                        translate(string.maketrans("", ""), "()");

        elem = tree.find("GDALWarpOptions/Transformer/GenImgProjTransformer/DstInvGeoTransform");
        # get inverse geotransform
        invgeotransform = gdal.InvGeoTransform(geoTransform);
        #print "551 DstInvGeoTransform:", str(invgeotransform[1]).translate(string.maketrans("", ""), "()");
        # convert proper string style and set to the DstInvGeoTransform element
        elem.text = str(invgeotransform[1]).translate(\
                                                string.maketrans("", ""), "()");

        # Overwrite element
        element = tree.getroot()
        # Write the modified elemements into VSI-file
        vsiFile = gdal.VSIFOpenL(self.warpedVRTName, 'w')
        #print tostring(element)
        gdal.VSIFWriteL(tostring(element), \
                            len(tostring(element)), 1, vsiFile);
        gdal.VSIFCloseL(vsiFile);
        
        newWarpedVRT = gdal.Open(self.warpedVRTName)
        #print newWarpedVRT.RasterXSize
        
        return newWarpedVRT


    def export(self, fileName, bandsList, dataType = gdal.GDT_Int16):
        copyFileName = "/vsimem/vrtCopy"
        vrtDsCopy = self.vrtDriver.CreateCopy(copyFileName, self.vrt);

        xsize = self.vrt.RasterXSize;
        ysize = self.vrt.RasterYSize;

        #create empty dataset with N bands
        vrtDrv = gdal.GetDriverByName("VRT");
        vrtDs = vrtDrv.Create("/vsimem/export_vrt.vrt", xsize, ysize, \
                                len(bandsList), dataType);
        vrtDs.SetGCPs(self.vrt.GetGCPs(), vrtDs.GetGCPProjection());
        vrtDs.SetGeoTransform(self.vrt.GetGeoTransform());
        vrtDs.SetMetadata(self.vrt.GetMetadata());
        vrtDs.SetProjection(self.vrt.GetProjection());

    	#populate the bands with source metadata
        for bn in range(len(bandsList)):
            metaItemKeys = self.rawVrt.GetRasterBand(bandsList[bn]).\
                            GetMetadata_Dict().keys();
            for iItemKey in metaItemKeys:
                vrtDs.GetRasterBand(bn+1).SetMetadataItem\
                (iItemKey,
                    self.rawVrt.GetRasterBand(bandsList[bn]).\
                    GetMetadataItem(iItemKey));

            BlockSize = vrtDsCopy.GetRasterBand(bandsList[bn]).GetBlockSize();

            bandSourceXML = '\
       	    <SimpleSource>\
              <SourceFilename relativeToVRT="0">%s</SourceFilename>\
              <SourceBand>%d</SourceBand>\
              <SourceProperties RasterXSize="%d" RasterYSize="%d" DataType="UInt16" BlockXSize="%d" BlockYSize="%d" />\
              <SrcRect xOff="0" yOff="0" xSize="%d" ySize="%d" />\
              <DstRect xOff="0" yOff="0" xSize="%d" ySize="%d" />\
            </SimpleSource>' % (copyFileName, bandsList[bn], xsize, ysize,\
                                 BlockSize[0], BlockSize[1], xsize, ysize, \
                                 xsize, ysize);

            vrtDs.GetRasterBand(bn+1).SetMetadataItem("source_0", \
                                                       bandSourceXML, \
                                                       "new_vrt_sources");

        vrtDs.FlushCache();

        tiffDrv = gdal.GetDriverByName("GTiff");
        copyDs = tiffDrv.CreateCopy(fileName + ".tif" ,vrtDs, 0);
        copyDs = None;
        vrtDsCopy = None;


    def get_domain(self):
        ''' Returns: Domain of the Nansat object '''
        return Domain(self.vrt)

    def __repr__(self):
        '''Prints basic info about the Nansat object to the terminal'''
        print '-'*40
        print 'Nansat object:'
        print 'Size : ' + str(self.vrt.RasterXSize) + ' x ' + str(self.vrt.RasterYSize)
        print 'Projection: ' + self.vrt.GetProjection()
        print 'GCP Projection: ' + self.vrt.GetGCPProjection()
        if self.vrt.GetProjection():
            print 'Corner Coordinates:'
            from gdalinfo import GDALInfoReportCorner
            hProj = osr.SpatialReference(self.vrt.GetProjection())
            hLatLong = hProj.CloneGeogCS()
            hTransform = osr.CoordinateTransformation(hProj, hLatLong)
            GDALInfoReportCorner( self.vrt, hTransform, "Upper Left", 0.0, 0.0 )
            GDALInfoReportCorner( self.vrt, hTransform, "Lower Left", 0.0, self.vrt.RasterYSize);
            GDALInfoReportCorner( self.vrt, hTransform, "Upper Right", self.vrt.RasterXSize, 0.0 );
            GDALInfoReportCorner( self.vrt, hTransform, "Lower Right", self.vrt.RasterXSize, \
                          self.vrt.RasterYSize );
        print '-'*20
        self.list_bands()
        print '-'*40
        return ''

    def __getitem__(self, bandNo):
        ''' Returns the band as a NumPy array, by overloading [] '''
        return self.get_GDALRasterBand(bandNo).ReadAsArray()



    def reproject_on_gcp(self, gcpImage, resamplingAlg = 0):
        ''' Reproject the object onto the input object with gcps
        NB! This is a test function required urgently for the open-wind
        project. It is tesed only on NCEP or GLOBAL DEM and 
        RADARSAT2 or MERIS images and should be refined and
        added to Nansat.reproject()

        Args:
            gcpImage : Nansat object of an image with GCPs
            resamplingAlg : integer, option for AutoCreateWarpedVRT

        Side effects:
            self.warpedVrt : new warped vrt
            self.vrt : replace current self.vrt to self.warpedVrt
        '''
        #name of VRT with 'fake' GCPs
        tmpVrtName = '/vsimem/vsiFileFakeGCP.vrt';
        
        #prepare pure lat/lon WKT
        proj4string = "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs"
        latlongSRS = osr.SpatialReference();
        latlongSRS.ImportFromProj4(proj4string);
        latlongWkt = latlongSRS.ExportToWkt();
        
        #get source SRS (either Projection or GCPProjection)
        srcWkt = self.vrt.GetProjection()
        if srcWkt == '':
            srcWkt = self.vrt.GetGCPProjection()

        #the transformer converts lat/lon to pixel/line of SRC image
        srcTransformer = gdal.Transformer(self.vrt, None, ['SRC_SRS='+srcWkt, 'DST_SRS='+latlongWkt])
        
        #get GCPs from DST image
        gcps = gcpImage.vrt.GetGCPs();
        
        #create 'fake' GCPs
        for g in gcps:
            #transform DST lat/lon to SRC pixel/line
            succ,point = srcTransformer.TransformPoint(1, g.GCPX, g.GCPY);
            srcPixel = point[0];
            srcLine = point[1];

            #swap coordinates in GCPs:
            #pix1/line1 -> lat/lon  =>=>  pix2/line2 -> pix1/line1
            g.GCPX = g.GCPPixel
            g.GCPY = g.GCPLine
            g.GCPPixel = srcPixel
            g.GCPLine = srcLine
        
        #make copy of the RAW VRT file and replace GCPs
        tmpVrt = self.vrtDriver.CreateCopy(tmpVrtName, self.rawVrt);

        #create 'fake' STEREO projection for 'fake' GCPs of SRC image
        srsString = "+proj=stere +lon_0=0 +lat_0=0 +k=1 +ellps=WGS84 +datum=WGS84 +no_defs ";
        stereoSRS = osr.SpatialReference();
        stereoSRS.ImportFromProj4(srsString);
        stereoSRSWKT = stereoSRS.ExportToWkt()
        tmpVrt.SetGCPs(gcps, stereoSRSWKT)
        tmpVrt.SetProjection('');
        tmpVrt = None
        
        #remove GeoTransfomr from SRC image
        # open XML content from VSI-file
        vsiFile = gdal.VSIFOpenL(tmpVrtName, "r")
        # get file size
        gdal.VSIFSeekL(vsiFile, 0, 2)
        vsiFileSize = gdal.VSIFTellL(vsiFile)
        gdal.VSIFSeekL(vsiFile, 0, 0) #fseek to start
        # read
        vsiFileContent = gdal.VSIFReadL(vsiFileSize, 1, vsiFile);
        gdal.VSIFCloseL(vsiFile);

        #find and remove GeoTransform
        tree = XML(vsiFileContent);
        elemGT = tree.find("GeoTransform");
        tree.remove(elemGT);

        # Write the modified elemements back into VSI-file
        vsiFile = gdal.VSIFOpenL(tmpVrtName, 'w')
        gdal.VSIFWriteL(tostring(tree), \
                            len(tostring(tree)), 1, vsiFile);
        gdal.VSIFCloseL(vsiFile);
        
        #create warped vrt out of tmp vrt
        tmpVrt = gdal.Open(tmpVrtName)
        rawWarpedVRT = gdal.AutoCreateWarpedVRT(tmpVrt, \
                          stereoSRSWKT, stereoSRSWKT, resamplingAlg);

        #change size and geotransform to fit the DST image
        self.warpedVRT = self._modify_warped_vrt(rawWarpedVRT,
                                gcpImage.vrt.RasterXSize, \
                                gcpImage.vrt.RasterYSize, \
                                (0, 1, 0, 0, 0, 1))
        self.vrt = self.warpedVRT;
