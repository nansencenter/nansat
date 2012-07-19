from dateutil.parser import parse
from struct import unpack

class Envisat():
    '''Methods shared between Envisat mappers'''
    def _set_envisat_time(self, gdalMetadata):
        ''' Get time from metadata, set time to VRT'''
        # set time
        productTime = gdalMetadata["SPH_FIRST_LINE_TIME"]
        self._set_time(parse(productTime))

    def read_scaling_gads(self, fileName, indeces):
        ''' Read Scaling Factor GADS to get scalings of MERIS L1/L2'''

        maxGADS = max(indeces) + 1

        # open file, find offset
        f = file(fileName, 'rt')
        headerLines = f.readlines(100)
        gadsDSNameString = 'DS_NAME="Scaling Factor GADS         "\n'
        if gadsDSNameString in headerLines:
            i1 = headerLines.index(gadsDSNameString)
            iGadsOffset = i1 + 3
            scalingGADSOffset = int(headerLines[iGadsOffset].\
                               replace('DS_OFFSET=', '').replace('<bytes>', ''))
        f.close()

        # fseek to gads, read all into a list
        f = file(fileName, 'rb')
        f.seek(scalingGADSOffset, 0)
        allGADSValues = []
        for i in range(maxGADS):
            fbString = f.read(4)
            fbVal = unpack('>f', fbString)[0]
            allGADSValues.append(fbVal)
        f.close()

        #get only values required for the mapper
        return [allGADSValues[i] for i in indeces]

    def read_GeolocationGrid_Offset(self, fileName):
        ''' Read ADSR offsets for ASAR

        See Also
        --------
            http://envisat.esa.int/handbooks/asar/CNTR6-6-9.htm#eph.asar.asardf.asarrec.ASAR_Geo_Grid_ADSR

        '''
        # open file, find offset
        f = file(fileName, 'rt')
        headerLines = f.readlines(150)
        gadsDSNameString = 'DS_NAME="GEOLOCATION GRID ADS        "\n'

        offsetDict = {}
        if gadsDSNameString in headerLines:
            gridOffset = headerLines.index(gadsDSNameString)
            offset = gridOffset + 3
            offsetDict["geolocationOffset"] = int(headerLines[offset].replace('DS_OFFSET=', '').replace('<bytes>', ''))
            offsetDict["longitudeOffset"] = offsetDict["geolocationOffset"] + 25 + (4*4)*11
            offsetDict["latitudeOffset"] = offsetDict["geolocationOffset"] + 25 + (4*3)*11
            offsetDict["incidentAngleOffset"] = offsetDict["geolocationOffset"] + 25 + (4*2)*11
            offset = gridOffset + 5
            offsetDict["numOfDSR"] = int(headerLines[offset].replace('NUM_DSR=', ''))
            offset = gridOffset + 6
            offsetDict["DSRsize"] = int(headerLines[offset].replace('DSR_SIZE=', '').replace('<bytes>', ''))
        else:
            print "Cannot find GEOLOCATION GRID ADS for Envisat"
            #for iText in headerLines:
            #    print iText
        f.close()
        return offsetDict

#m = MERIS();
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_1PNPDK20110817_110451_000004053105_00339_49491_7010.N1', range(7, 22))
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(7, 20) + [20, 21, 22, 20])
#print m.read_scaling_gads('/Data/sat/GDAL_test/MER_FRS_2CNPDK20110503_105820_000000813102_00109_47968_7906.N1', range(33, 46) + [46, 47, 48, 46])
