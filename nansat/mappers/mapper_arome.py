from nansat.mappers.mapper_netcdfcf import Mapper as NetcdfCF

class Mapper(NetcdfCF):

    def __init__(self, *args, **kwargs):
    
        super(Mapper, self).__init__(*args, **kwargs)
