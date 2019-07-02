from nansat.mappers.sentinel1 import Sentinel1
from nansat.mappers.mapper_netcdf_cf import Mapper as NetCDF_CF_Mapper

class Mapper(Sentinel1, NetCDF_CF_Mapper):

    def __init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs):
        NetCDF_CF_Mapper.__init__(self, filename, gdal_dataset, gdal_metadata, *args, **kwargs)
        Sentinel1.__init__(self, filename)
        self.add_calibrated_nrcs()

    def add_calibrated_nrcs(self):
        polarizations = [self.ds.polarisation[i:i+2] for i in range(0,len(self.ds.polarisation),2)]
        for pol in polarizations:
            amp_fn = 'NETCDF:"' + self.input_filename + '":Amplitude_%s' %pol
            bdict_amp = self._get_band_from_subfile(amp_fn)
            s0_fn = 'NETCDF:"' + self.input_filename + '":sigmaNought_%s' %pol
            bdict_s0 = self._get_band_from_subfile(s0_fn)
            src = [
                    bdict_amp['src'],
                    bdict_s0['src']
                ]
            dst = {
                    'wkv': 'surface_backwards_scattering_coefficient_of_radar_wave',
                    'PixelFunctionType': 'Sentinel1Calibration',
                    'polarization': pol,
                    'suffix': pol,
                }
            self.create_band(src, dst)
            self.dataset.FlushCache()




