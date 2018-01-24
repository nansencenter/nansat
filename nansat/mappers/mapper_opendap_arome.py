from nansat.mappers.mapper_arome import Mapper as MapperArome
from nansat.mappers.opendap import Opendap
from nansat.exceptions import WrongMapperError

class Mapper(Opendap, MapperArome):

    baseURLs = ['http://thredds.met.no/thredds/catalog/arome25/catalog.html']

    def __init__(self, *args, **kwargs):

        raise WrongMapperError('Mapper is under development...')

        fn = args[0]
        #ds = args[1] - None
        #mm = args[2] - None

        self.test_mapper(fn)
