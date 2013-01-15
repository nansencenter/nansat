# Name:        mapper_modisL1
# Purpose:     Mapping for MODIS-L1 data
# Authors:      Anton Korosov
# Licence:      This file is part of NANSAT. You can redistribute it or modify
#               under the terms of GNU General Public License, v.3
#               http://www.gnu.org/licenses/gpl-3.0.html

from vrt import VRT, gdal, parse

class Mapper(VRT):
    ''' VRT with mapping of WKV for MODIS Level 1 (QKM, HKM, 1KM) '''

    def __init__(self, fileName, gdalDataset, gdalMetadata):
        ''' Create MODIS_L1 VRT '''
        #get 1st subdataset and parse to VRT.__init__() for retrieving geo-metadata
        gdalSubDataset = gdal.Open(gdalDataset.GetSubDatasets()[0][0])
       
        #list of available modis names:resolutions
        modisResolutions = {'MYD02QKM':250, 'MOD02QKM':250,
                            'MYD02HKM':500, 'MOD02HKM':500,
                            'MYD021KM':1000, 'MOD021KM':1000};
        
        #should raise error in case of not MODIS_L1
        mResolution = modisResolutions[gdalMetadata["SHORTNAME"]];

        # create empty VRT dataset with geolocation only
        VRT.__init__(self, gdalSubDataset)
       
        subDsString = 'HDF4_EOS:EOS_SWATH:"%s":MODIS_SWATH_Type_L1B:%s'
                
        #provide all mappings
        metaDict250SF = ['EV_250_RefSB']
        metaDict250 = [
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_250_RefSB'), 'SourceBand': 1}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '645'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_250_RefSB'), 'SourceBand': 2}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '858'}}
        ];
        
        metaDict500SF = ['EV_250_Aggr500_RefSB', 'EV_500_RefSB']
        metaDict500 = [
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_250_Aggr500_RefSB'), 'SourceBand': 1}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '645'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_250_Aggr500_RefSB'), 'SourceBand': 2}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '858'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_RefSB'), 'SourceBand': 1}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '469'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_RefSB'), 'SourceBand': 2}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '555'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_RefSB'), 'SourceBand': 3}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '1240'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_RefSB'), 'SourceBand': 4}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '1640'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_RefSB'), 'SourceBand': 5}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '2130'}}
        ];

        metaDict1000SF = ['EV_250_Aggr1km_RefSB', 'EV_500_Aggr1km_RefSB', 'EV_1KM_RefSB', 'EV_1KM_Emissive']
        metaDict1000 = [
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_250_Aggr1km_RefSB'), 'SourceBand': 1}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '645'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_250_Aggr1km_RefSB'), 'SourceBand': 2}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '858'}},

        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_Aggr1km_RefSB'), 'SourceBand': 1}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '469'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_Aggr1km_RefSB'), 'SourceBand': 2}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '555'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_Aggr1km_RefSB'), 'SourceBand': 3}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '1240'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_Aggr1km_RefSB'), 'SourceBand': 4}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '1640'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_500_Aggr1km_RefSB'), 'SourceBand': 5}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '2130'}},

        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 1}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '412'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 2}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '443'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 3}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '488'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 4}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '531'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 5}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '551'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 6}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '667'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 7}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '667'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 8}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '678'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 9}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '678'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 10}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '748'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 11}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '869'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 12}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '905'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 13}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '936'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 14}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '940'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_RefSB'), 'SourceBand': 15}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '1375'}},

        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 1}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '3750'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 2}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '3959'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 3}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '3959'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 4}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '4050'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 5}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '4465'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 6}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '4515'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 7}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '6715'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 8}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '7325'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 9}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '8550'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 10}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '9730'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 11}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '11030'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 12}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '12020'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 13}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '13335'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 14}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '13635'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 15}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '13935'}},
        {'src': {'SourceFilename': subDsString % (fileName, 'EV_1KM_Emissive'), 'SourceBand': 16}, 'dst': {'wkv': 'toa_outgoing_spectral_radiance', 'wavelength': '14235'}}
        ];

        # get proper mapping depending on resolution
        metaDict = {
            250: metaDict250,
            500: metaDict500,
            1000: metaDict1000,
        }[mResolution];
        # get proper mapping depending on resolution
        metaDictSF = {
            250: metaDict250SF,
            500: metaDict500SF,
            1000: metaDict1000SF,
        }[mResolution];
        
        # read all scales/offsets
        rScales = {}
        rOffsets = {}
        for sf in metaDictSF:
            dsName = subDsString % (fileName, sf)
            ds = gdal.Open(dsName)
            rScales[dsName] = map(float, ds.GetMetadataItem('radiance_scales').split(','))
            rOffsets[dsName] = map(float, ds.GetMetadataItem('radiance_offsets').split(','))
            self.logger.debug('radiance_scales: %s' % str(rScales))
            
        # add 'band_name' to 'parameters'
        for bandDict in metaDict:
            SourceFilename = bandDict['src']['SourceFilename']
            SourceBand = bandDict['src']['SourceBand']
            bandDict['dst']['suffix'] = bandDict['dst']['wavelength']
            scale = rScales[SourceFilename][SourceBand-1]
            offset = rOffsets[SourceFilename][SourceBand-1]
            self.logger.debug('band, scale, offset: %s_%d %s %s' % (SourceFilename, SourceBand, scale, offset))
            bandDict['src']['ScaleRatio'] = scale
            bandDict['src']['ScaleOffset'] = offset
            
        # add bands with metadata and corresponding values to the empty VRT
        self._create_bands(metaDict)

        productDate = gdalMetadata["RANGEBEGINNINGDATE"]
        productTime = gdalMetadata["RANGEENDINGTIME"]
        self._set_time(parse(productDate+' '+productTime))
