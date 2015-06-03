#-------------------------------------------------------------------------------
# Name:
# Purpose:      
#
# Author:       Morten Wergeland Hansen
# Modified:	Morten Wergeland Hansen
#
# Created:	02.06.2015
# Last modified:03.06.2015 13:23
# Copyright:    (c) NERSC
# License:      
#-------------------------------------------------------------------------------
import os.path
import json
import re
import numpy as np
import datetime
from collections import OrderedDict

from nansat.vrt import VRT
from nansat.tools import WrongMapperError

class Mapper(VRT):
    def __init__(self, fileName, *args, **kwargs):

        try:
            rvl = gsar(fileName)
        except ValueError:
            raise WrongMapperError

        channels_info = [rvl.getinfo(channel = 0)]
        try:
            channels_info.append(rvl.getinfo(channel = 1))
        except KeyError:
            pass # Only one channel in the dataset

        # Read data into numpy record array (in memory, violates lazy
        # operations approach...)
        channels_data = [rvl.getdata(channel = 0)]
        if len(channels_info)==2:
            channels_data.append(rvl.getdata(channel = 1))
            np.testing.assert_array_equal(channels_data[0]['LONGITUDE'],
                    channels_data[1]['LONGITUDE'])
            np.testing.assert_array_equal(channels_data[0]['LATITUDE'],
                    channels_data[1]['LATITUDE'])
            np.testing.assert_array_equal(channels_data[0]['INCANGLE'],
                    channels_data[1]['INCANGLE'])
            np.testing.assert_array_equal(channels_data[0]['HEADING'],
                    channels_data[1]['HEADING'])
            np.testing.assert_array_equal(channels_data[0]['PREDDOPFREQ'],
                    channels_data[1]['PREDDOPFREQ'])

        longitude = channels_data[0]['LONGITUDE']
        latitude = channels_data[0]['LATITUDE']
        VRT.__init__(self, lon = longitude, lat = latitude)

        incVRT = VRT(array = channels_data[0]['INCANGLE'], lon = longitude,
                lat = latitude)
        azVRT = VRT(array = np.mod(channels_data[0]['HEADING'] + 90., 360.),
                lon = longitude, lat = latitude)
        dcpVRT = VRT(array = channels_data[0]['PREDDOPFREQ'], lon = longitude,
                lat = latitude)

        metaDict = []
        self.bandVRTs['incVRT'] = incVRT
        self.bandVRTs['azVRT'] = azVRT
        self.bandVRTs['dcpVRT'] = dcpVRT
        metaDict.append({
            'src': {
                'SourceFilename': self.bandVRTs['incVRT'].fileName,
                'SourceBand': 1
            },
            'dst': {
                'wkv': 'angle_of_incidence',
                'name': 'incidence_angle'
            }
        })
        metaDict.append({
            'src': {
                'SourceFilename': self.bandVRTs['azVRT'].fileName,
                'SourceBand': 1
            },
            'dst': {
                'wkv': 'sensor_azimuth_angle',
                #'name': 'SAR_look_direction'
            }
        })
        metaDict.append({
            'src': {
                'SourceFilename': self.bandVRTs['dcpVRT'].fileName,
                'SourceBand': 1
            },
            'dst': {
                'wkv': 'predicted_surface_backwards_doppler_centroid_frequency_shift_of_radar_wave',
            }
        })

        for i, channel in enumerate(channels_data):
            nrcsVRT = VRT(array = channel['NRCS'], lon = longitude,
                lat = latitude)
            dcVRT = VRT(array = channel['DOPFREQ'], lon = longitude,
                lat = latitude)
            dc_std_VRT = VRT(array = channel['DOPFREQSTD'], lon = longitude,
                lat = latitude)

            pol = channels_info[i].polarization

            self.bandVRTs[u'nrcs_'+pol+u'_VRT'] = nrcsVRT
            self.bandVRTs[u'dc_'+pol+u'_VRT'] = dcVRT
            self.bandVRTs[u'dc_std_'+pol+u'_VRT'] = dc_std_VRT

            metaDict.append({
                'src': {
                    'SourceFilename': self.bandVRTs['nrcs_'+pol+'_VRT'].fileName,
                    'SourceBand': 1,
                },
                'dst': {
                    'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                    'suffix': pol,
                    'polarization': pol,
                }
            })
            metaDict.append({
                'src': {
                    'SourceFilename': self.bandVRTs['dc_'+pol+'_VRT'].fileName,
                    'SourceBand': 1
                },
                'dst': {
                    'wkv': 'surface_backwards_doppler_centroid_frequency_shift_of_radar_wave',
                    'suffix': pol,
                    'polarization': pol,
                }
            })
            metaDict.append({
                'src': {
                    'SourceFilename': self.bandVRTs['dc_std_'+pol+'_VRT'].fileName,
                    'SourceBand': 1
                },
                'dst': {
                    'wkv': 'standard_deviation_of_surface_backwards_doppler_centroid_frequency_shift_of_radar_wave',
                    'suffix': pol,
                    'polarization': pol,
                }
            })

        self._create_bands(metaDict)

        # Add pixelfunction for getting the Doppler anomaly
        for i, channel in enumerate(channels_data):
            pol = channels_info[i].polarization
            src = [{
                    'SourceFilename': self.bandVRTs['dc_'+pol+'_VRT'].fileName,
                    'SourceBand': 1
                },{
                    'SourceFilename': self.bandVRTs['dcpVRT'].fileName,
                    'SourceBand': 1
                }]
            dst = {
                    'wkv': 'surface_backwards_doppler_frequency_shift_of_radar_wave_due_to_surface_velocity',
                    'PixelFunctionType': 'diff',
                    'suffix': pol,
                    'polarization': pol,
                }

            self._create_band(src, dst)
            self.dataset.FlushCache()


        self.dataset.SetMetadataItem('mapper', 'gsar')

def cohist(zdata):
    return np.histogram(abs(zdata.ravel()), bins = np.linspace(0, 1, 101))

def phase(zz):
    return np.arctan2(zz.imag, zz.real)

def dtemplate(s):
    def ttrans(typespec):
        if type(typespec) is OrderedDict:
            return [ (str(key), ttrans(val)) for key, val in typespec.iteritems() ]
        if typespec == 'float': return 'float32'
        if typespec == 'complex': return 'complex64'

    # Dirty hack: Strip off '(1)' array indicators. non-singleton arrays must be dealt with
    etemp = re.sub('(\w+)', r'"\1"', s.replace('(1)', ''))
    dtemp = json.loads(etemp, object_pairs_hook=OrderedDict)

    return [ (str(key), ttrans(val)) for key, val in dtemp.iteritems() ]

class ytime:
    pattern = re.compile(r'(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})(\.\d+)?Z?')

    def __init__(self, s):
        m = self.pattern.match(s)
        if m is None:
            raise ValueError('Not a valid ytime string')
        groups = m.groups()
        if groups[-1] is None:
            usec = 0
        else:
            usec = int(round(1e6 * float(groups[-1])))  # Loses sub-usec digits

        self.dtime = datetime.datetime(*map(int, groups[:-1]) + [usec])

    def __repr__(self):
        return repr(self.dtime)

    def __str__(self):
        s = self.dtime.strftime('%FT%T')
        if self.dtime.microsecond != 0:
            s += ".%06d" % self.dtime.microsecond
        s += 'Z'

        return s



class struct:
    t_map = (   None,           # 0 used for no-data files
        np.uint8,               # 1
        np.int16,               # 2 signed short
        np.int32,               # 3 signed long
        np.float32,             # 4 single
        np.float64,             # 5 double
        np.complex64,           # 6 single-precision complex
        None,                   # 7
        None,                   # 8
        np.complex128,          # 9 double-precision complex
        None,                   # 10
        None,                   # 11
        np.uint16,              # 12 unsigned short
        np.uint32,              # 13 unsigned long
        np.int64,               # 14 signed long long
        np.uint64,              # 15 unsigned long long
    )

    def __init__(self, d):
        for k, v in d.iteritems():
            # print "k: {}, v: {}".format(k, v)
            if k.lower() == 'datatype':
                v = self.t_map[v]
            elif type(v) is dict:
                v = struct(v)
            elif type(v) in [unicode, str]:
                try:
                    v = ytime(v)
                except ValueError:
                    pass
            setattr(self, k.lower(), v)

        if hasattr(self, 'datatemplate'):
            self.datatype = np.dtype(dtemplate(self.datatemplate))
        elif hasattr(self, 'datatype'):
            self.datatype = np.dtype(self.datatype)

    def __repr__(self, level=4):
        return "%s(%s)" % (self.__class__.__name__,
            (",\n" + " " * level).join("%s=%s" % (k, (v.__repr__(level+4) \
                    if isinstance(v, struct) else v)) \
                            for k, v in self.__dict__.iteritems()))

class gsar:

    def __init__(self, fname):

        fnbase, ext = os.path.splitext(fname)

        if ext.lower() in ['.xml', '.json']:
            raise ValueError("Open GSAR file, not metadata file")

        if ext.lower() != '.gsar':
            raise ValueError("Not a GSAR file")

        open(fname)     # not used, but will cause error if permissions are denied

        self.datafile = fname

        if os.path.exists(fname + '.json'):
            self.metafile = fname + '.json'
        elif os.path.exists(fnbase + '.json'):
            self.metafile = fnbase + '.json'
        else:
            raise IOError("JSON Metadata file not found")

        self.__info = json.load(open(self.metafile))

    def __getitem__(self, ix):
        if len(ix) != 2:
            raise ValueError("need y and x indices")
        
        if type(ix[0]) is slice:
            ypos = ix[0].start or 0
            ysize = ix[0].stop - ypos
        else:
            # ypos, ysize = 0, ix[0]
            ypos, ysize = ix[0], 1
        
        if type(ix[1]) is slice:
            xpos = ix[1].start or 0
            xsize = ix[1].stop - xpos
        else:
            # xpos, xsize = 0, ix[1]
            xpos, xsize = ix[1], 1

        return self.getdata(ysize, xsize, ypos, xpos)


    def getinfo(self, channel=0):
        key = 'CHANNEL_' + str(channel)
        info = struct(self.__info[key])
        # info.datatype = self.t_map[info.datatype]

        return info

    def getdata(self, ysize=None, xsize=None, ypos=0, xpos=0, channel=0):

        # last/inner dimension is x (x indices run faster)
        info = self.getinfo(channel=channel)

        if ysize is None: ysize = info.ysize
        if xsize is None: xsize = info.xsize

        assert ypos >= 0 and xpos >= 0, 'FIXME: allow negative xpos/ypos in getdata'
        assert ysize >= 0 and xsize >= 0, 'x/ysize must be positive nonzero'

        data = np.zeros((ysize, xsize), dtype=info.datatype)

        if ypos >= info.ysize or xpos >= info.xsize: return data

        nrows_toread = min(info.ysize, ypos + ysize) - ypos

        # Now, read up an array with all the rows from the file which are necessary
        # (wasteful, I know). Will fix later
        df = open(self.datafile, 'rb')
        if ypos > 0: df.seek(ypos * info.xsize * info.datatype.itemsize)

        buf = np.fromfile(df, info.datatype, info.xsize * nrows_toread)
        buf.shape = (nrows_toread, info.xsize)

        if xpos == 0 and xsize == info.xsize:
            assert ypos >= 0, 'FIXME before allowing negative ypos'
            if ypos+ysize < info.ysize:
                # Data read are exactly data requested
                return buf
            data[:nrows_toread,:] = buf
        elif xpos+xsize < info.xsize:
            data[:nrows_toread,:] = buf[:,xpos:xpos+xsize]
        else:
            data[:nrows_toread,:info.xsize-xpos] = buf[:,xpos:]

        return data



