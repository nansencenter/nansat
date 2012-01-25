#------------------------------------------------------------------------------#
# Name:    wkv
# Purpose:
#           Constaruct wkv from XML files.
#
# Author:      asumak
#
# Created:     25.07.2011
# Copyright:   (c) asumak 2011
# Licence:
#------------------------------------------------------------------------------#

import os
import sys

try:
    from osgeo import gdal
except ImportError:
    import gdal

import numpy as np

from xml.etree.ElementTree import *

class MatrixSizeError(Exception):
    '''Error for improper matrix size'''
    pass;

class NoFlagError(Exception):
    '''Error occures when the flag is not recorded'''
    pass;


class WKV():
    '''Set WKV

    Set proper WKV, which is defined in nansat_filetypes.xml and wkv.xml,
    based on the given fileType.

    The fileType should be in the list of nansat_filetypes.xml,
        otherwise get an error.
    The given fileType corresponds to fileType-name in nansat_filetypes.xml.
    The band locations, wkv, and each parameter are fetched
        from nansat_filetypes.xml.
    Some parameter values are obtained from
        not nansat_filetypes.xml but the dataset.
    name, longname and unit corrsepond to the wkv in nansat_filetypes.xml
        are fetched from wkv.xml.

    '''
    def __init__(self, ds, fileType):
        '''Construct WKV object

        Set the locations of two XMLfiles (nansat_filetypes.xml and wkv.xml)
        Currently suppose these files are in the same folder as this file.
        Read wkv in wkv.xml
        Set wkv_list

        Args:
            ds : dataset
            fileType :"Radarsat2_singlechannel","Radarsat2_dualpol",
                      "Radarsat2_quadpol",
                      "ASAR_1channel", "ASAR_dualpol",
                      "NCEP_GRIB", "HIRLAM_GRIB",
                      "MYD02QKM", "MYD02HKM", "MOD021KM"

        Side effects:
            set attributes; fileName_filetypes, fileName_wkv,
            wkvXML and wkv_list

        '''
        self.fileName_filetypes = os.path.dirname(__file__) \
                                + "/nansat_filetypes.xml";
        self.fileName_wkv = os.path.dirname(__file__) + "/wkv.xml";

        fd = file(self.fileName_wkv, "rb");
        elem = ElementTree(file=fd).getroot();
        self.wkvXML = elem;
        self.wkv_list = self.get_wkvList(ds, fileType);

    def get_wkvList(self, ds, fileType):
        '''Set wkv_list

        Fetch band information from nansat_filetypes.xml
            which corresponds to a given fileType.
        Set the information in a numpyarray.
        Two numpyarraies are constructed, default and option.
        "default" is common band information for any sensors.
        The keys of "option" depends on sensors.
        Finally these two arrays are combined and create an array "wkv_list".

        Args:
            ds : dataset
            fileType :"Radarsat2_singlechannel","Radarsat2_dualpol",
                      "Radarsat2_quadpol",
                      "MERIS_L1", "MERIS_L2",
                      "ASAR_1channel", "ASAR_dualpol",
                      "NCEP_GRIB", "HIRLAM_GRIB",
                      "MYD02QKM", "MYD02HKM", "MOD021KM"

        Returns:
            wkv_list : numpyarray
                    wkv_list["index"]
                    wkv_list["name"]
                    wkv_list["longname"]
                    wkv_list["sub"]
                    wkv_list["src"]
                    wkv_list["wkv"]
                    wkv_list["unit"]
                        + parameters (e.g. wkv_list["wavelength"])

        Raises:
            AttributeError : occurs when cannot get bandXML[i].wkv
            NoFlagError: occurs when the flag has no method to fetch the value

        '''
        # fetch band information corresponding to the fileType

        fd = file(self.fileName_filetypes, "rb");
        element = ElementTree(file=fd).getroot();

        for e1 in list(element):
            keyList1 =  e1.keys();
            for i in range(len(keyList1)):
                if e1.get(keyList1[i]) == fileType:
                    # create a nampyarray for default (common) wkv
                    # with len(bandXML) bands
                    default = np.ndarray(len(e1.findall("band")), dtype = \
                            {"names": ["index", "name", "longname",\
                                        "sub", "src", "wkv", "unit"], \
                            "formats": ["a2", "a20", "a100", \
                                        "a3", "a2", "a40", "a10"]});
                    # create a nampyarray for optional wkv with len(bandXML) bands
                    prmNames = [];
                    prmFormats = [];
                    for e2 in list(e1):
                        for e3 in list(e2):
                            prmNameValue = e3.attrib["name"]
                            if prmNameValue in prmNames:
                                pass;
                            else:
                                prmNames.append(prmNameValue);
                                prmFormats.append("a20");

                    option = np.ndarray(len(e1.findall("band")), \
                                            dtype = {"names":prmNames, \
                                                     "formats":prmFormats});

                    # Clear option array
                    for j in range(len(e1.findall("band"))):
                        for k in range(len(prmNames)):
                            option[j][k] = None;

                    j = 0;
                    for e2 in list(e1):
                        # set common wkv
                        default["index"][j] = j + 1;
                        default["name"][j] = fileType;
                        default["sub"][j] = str(e2.attrib["sub"]);
                        default["src"][j] = str(e2.attrib["src"]);
                        default["wkv"][j] = str(e2.attrib["wkv"]);
                        # get corresponding longname and unit from wkv.xml
                        longname, unit = self.get_wkv(default["wkv"][j]);
                        default["longname"][j] = longname;
                        default["unit"][j] = unit;

                        # get corresponding longname and unit from wkv.xml
                        k = 0;
                        for e3 in list(e2):
                            # set optional wkv (parameters)
                            for iPar in range(len(prmNames)):
                                if prmNames[iPar] == e3.attrib["name"]:
                                    option[prmNames[iPar]][j] = e3.attrib["value"];
                            # Set optional parameters
                            #   that is fetched from dataset
                            if option[prmNames[k]][j] =="":
                                flag = e3.attrib["flag"];
                                # get Radarsat Polarization
                                if flag == "RadarsatPolarization":
                                    option[prmNames[k]][j] = \
                                                self.get_RadarsatPolarization\
                                                (ds, int(default["src"][j]));
                                # get Radarsat DataType
                                elif flag == "RadarsatDataType":
                                    option[prmNames[k]][j] = \
                                                   self.get_RadarsatDataType\
                                                   (ds, int(default["src"][j]));
                                # get ASAR Polarization
                                elif flag == "ASAR_Polarization":
                                    option[prmNames[k]][j] = \
                                                self.get_ASARPolarization(ds);
                                # get ASAR DataType
                                elif flag == "ASARDataType":
                                    option[prmNames[k]][j] = \
                                                self.get_ASARDataType(ds);
                                # show error if the flag is not recognized
                                else:
                                    raise NoFlagError("get_wkv(): Flag is not recorded.");
                            k += 1;
                        j += 1;
                break;

        # combine default wkv and optional wkv
        wkv_list = self.add_keys(default, option);
        return wkv_list;

    def get_wkv(self, wkvName):
        '''Search for an item number, that corresponds to wkvName, from wkv.xml

        Args:
            wkvName: "wkv" value in nansat_filetype.xml

        Return:
            longname : string
            unit : string

        '''
        for element in list(self.wkvXML):
            if element[0].text == wkvName:
                longname = element[1].text;
                unit = element[2].text;
                break;

        return longname, unit;

    def add_keys(self, vctr1, vctr2):
        '''Create a composite numpyarray from two numpyarrays

        Check the matrices sizes.
        Get names and datatypes from the matrices and create composite lists.
        Create an empty composite numpyarray based on the compsite lists.
        Matrix(n1 x 1 x m) + Matrix(n2 x 1 x m) --> Matrix((n1 + n2) x 1 x m)
        Set the corresponding values into the comopsite array.

        Args:
            vctr1: numpyarray whose size is "n1 x 1 x m"
            vctr2: numpyarray whose size is "n2 x 1 x m"

        Return:
            vctrComp: numpy array whose size is "(n1 + n2) x 1 x m"

        Raises:
            MatrixSizeError : occurs when the two vector sizes
                              (numbers of the bands) are different.

        '''
        # Check the sizes of two matrices
        if len(vctr1) != len(vctr2):
            raise MatrixSizeError("add_keys(): \
                                   the sizes of two vectors are different.");

        # Create composite lists for keys and datatypes
        keyNames = [];
        typeFormats = [];

        for i in range(len(vctr1.dtype.names)):
            keyNames.append(vctr1.dtype.names[i]);
            dt = vctr1.dtype[i];
            typeFormats.append(self.convert_datatype(str(dt)));

        for i in range(len(vctr2.dtype.names)):
            keyNames.append(vctr2.dtype.names[i]);
            dt = vctr2.dtype[i];
            typeFormats.append(self.convert_datatype(str(dt)));

        # create a numpyarray
        # row :  7 (number of common wkv )
        #        + number of optional wkv (len(vctr2.dtype.names))
        # data(band) number : len(vctr1)
        vctrComp = np.ndarray(len(vctr1), dtype = {"names": keyNames, \
                                                   "formats":typeFormats})

        # Set the corresponding values from vctr 1 and vctr2
        for j in range(len(vctr1.dtype.names)):
            for i in range(len(vctr1)):
                vctrComp[vctr1.dtype.names[j]][i] \
                        = vctr1[vctr1.dtype.names[j]][i];

        for j in range(len(vctr2.dtype.names)):
            for i in range(len(vctr1)):
                vctrComp[vctr2.dtype.names[j]][i] \
                        = vctr2[vctr2.dtype.names[j]][i];
        return vctrComp;

    def convert_datatype(self, datatype):
        ''' Convert String data type from "|S10" to "a10"

        Args :
            datatype: datatype

        Returns :
            dt: datatype

        Raises:
            TypeError: occurs when the argunment (datatype) is not string

        '''
        if datatype[0:2] == "|S" :
            dt = datatype.replace("|S", "a", 1);
        else:
            raise TypeError("convert_datatype():the data type is not string.");
        return dt;

    def get_RadarsatPolarization(self, ds, bn):
        ''' Get Polarization for a given Radarsat band

        Parameters:
            ds: dataset
            bn: band number

        Returns:
            value: HH / VV / HV / VH

        Raises:
            KeyError: occurs when band.GetMetadata() does not have
                        "POLARIMETRIC_INTERP" key

        '''
        band = ds.GetRasterBand(bn);
        bandMetadata = band.GetMetadata();
        if bandMetadata.has_key("POLARIMETRIC_INTERP"):
            value = bandMetadata["POLARIMETRIC_INTERP"];
        else:
            raise KeyError("get_RadarsatPolarization():\
                    Cannot find 'POLARIMETRIC_INTERP' key in the band metadata");
        return value;

    def get_RadarsatDataType(self, ds, bn):
        ''' Get DataType for a given Radarsat band

        Parameters :
            ds: dataset
            bn: band number

        Returns :
            value: Detected / Complex

        Raises:
            AttributeError: occurs when band.DataType does not have
                            attribute "DataTypeIsComplex"

        '''
        band = ds.GetRasterBand(bn);
        if gdal.DataTypeIsComplex (band.DataType) == 0:
            value = "Detected";
        elif gdal.DataTypeIsComplex (band.DataType) != 0:
            value = "Complex";
        else:
            raise AttributeError("get_type(): no attribute 'DataTypeIsComplex'");
        return value;

    def get_ASARPolarization(self, ds):
        ''' Get Polarization for a given ASAR dataset

        Parameters :
            ds: dataset

        Returns :
            value: H/H, V/V, H/V or V/H

        Raises:
            KeyError: occurs when ds.GetMetadata() does not have
                            the key "SPH_MDS1_TX_RX_POLAR"

        '''
        dsMetadata = ds.GetMetadata();
        if dsMetadata.has_key("SPH_MDS1_TX_RX_POLAR"):
            value = dsMetadata["SPH_MDS1_TX_RX_POLAR"];
        else:
            raise KeyError("get_ASAR_Polarization():\
                    Cannot find 'SPH_MDS1_TX_RX_POLAR' key in the dataset metadata");
        return value;

    def get_ASARDataType(self, ds):
        ''' Get DataType for a given ASAR dataset

        Parameters :
            ds: dataset

        Returns :
            value: ???

        Raises:
            KeyError: occurs when ds.GetMetadata() does not have
                           the key "SPH_DATA_TYPE"

        '''
        dsMetadata = ds.GetMetadata();
        ##for i in dsMetadata:
        ##    print i, "-- ", ds.GetMetadataItem(i);
        if dsMetadata.has_key("SPH_DATA_TYPE"):
            value = dsMetadata["SPH_DATA_TYPE"];
        else:
            raise KeyError("get_ASARDataType():\
                    Cannot find 'SPH_DATA_TYPE' key in the DSmetadata");
        return value;







