#-----------------------------------------------------------------------------
# Name:        setup.py
# Purpose:
#
# Author:      asumak
#
# Created:     17.06.2013
# Copyright:   (c) asumak 2012
# Licence:     <your licence>
# =========  !! NB !! HOW TO DO FOR MAC USERS??  ==========
#-----------------------------------------------------------------------------

from subprocess import Popen
import subprocess
import sys
import errno
import os

# Check if numpy, gdal, and basemap packages are installed
try:
    import numpy
except ImportError:
    raise ImportError("Nansat requires numpy, which should be installed separately")
try:
    from osgeo import gdal, osr, ogr
except ImportError:
    try:
        import gdal
        import osr
        import ogr
    except ImportError:
        raise ImportError("Nansat requires gdal, which should be installed separately")
try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
        raise ImportError("Nansat requires Basemap, which should be installed separately")

NAME                = 'nansat'
MAINTAINER          = "Nansat Developers"
MAINTAINER_EMAIL    = "nansat-dev@googlegroups.com"
DESCRIPTION         = "A scientist friendly Python toolbox for processing 2D satellite Earth observation data"
LONG_DESCRIPTION    = "A scientist friendly Python toolbox for processing 2D satellite Earth observation data"
URL                 = "https://github.com/nansencenter/nansat"
DOWNLOAD_URL        = "https://github.com/nansencenter/nansat"
LICENSE             = "GNU General Public License"
CLASSIFIERS         = '***'  # filter(None, CLASSIFIERS.split('\n'))
AUTHOR              = ("Asuka Yamakawa, Anton Korosov, Morten W. Hansen, Kunt-Frode Dagestad")
AUTHOR_EMAIL        = "nansat-dev@googlegroups.com"
PLATFORMS           = ["UNKNOWN"]
MAJOR               = 0
MINOR               = 6
MICRO               = 0
ISRELEASED          = False
VERSION             = '%d.%d-dev.%d' % (MAJOR, MINOR, MICRO) # Remember to remove "dev" when releasing
REQS                = [
                        "scipy",
                        "matplotlib",
                        "Pillow",
                    ]

#----------------------------------------------------------------------------#
#                   Prepare compilation of C pixel functions
#----------------------------------------------------------------------------#
skip_compile = False
libraries = []
include_dirs = []
library_dirs = []
if sys.platform == 'win32':
    extra_compile_args = ['-nologo', '-DLL']
    path = os.environ['LIB'].split(';')
    for iFolder in path:
        try:
            files = os.listdir(iFolder)
        except:
            extra_link_args = []
        else:
            if 'gdal_i.lib' in files:
                extra_link_args = [iFolder + '/gdal_i.lib']
                break
else:
    extra_compile_args = ['-fPIC', '-Wall', '-Wno-long-long', '-pedantic', '-O3']
    extra_link_args = [] # not used currently

def _ask_gdal_config(resultlist, option, result_prefix):
    try:
        p = Popen(['gdal-config', option], stdout=subprocess.PIPE)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise
    else:
        t = p.stdout.read().decode().strip()
        if p.wait() != 0:
            return
        res = t.split()
        res = filter(lambda x: x.startswith(result_prefix), res)
        # '-I/usr/...' -> '/usr/...'
        res = [x[len(result_prefix):] for x in res]
        resultlist[:] = res

def use_gdal_config():
    _ask_gdal_config(include_dirs, '--cflags', '-I')
    _ask_gdal_config(library_dirs, '--libs', '-L')
    _ask_gdal_config(libraries,    '--libs', '-l')

try:
    use_gdal_config()
except Exception as e:
    print 'WARNING: gdal-config could not be called, ' +\
          'pixel functions will not be available.'
    print 'Error details follow:'
    print repr(e)
    skip_compile = True

#----------------------------------------------------------------------------#
#                               Install package
#----------------------------------------------------------------------------#
from setuptools import setup
from distutils.extension import Extension
from distutils.errors import CCompilerError, DistutilsExecError,\
    DistutilsPlatformError

# the following is adapted from simplejson's setup.py
if sys.platform == 'win32' and sys.version_info > (2, 6):
    # 2.6's distutils.msvc9compiler can raise an IOError when failing to
    # find the compiler
    # It can also raise ValueError http://bugs.python.org/issue7511
    ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError,
                  IOError, ValueError)
else:
    ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError)

def run_setup(skip_compile):
    if skip_compile:
        kw = dict()
    else:
        kw = dict(
            ext_modules = [
                Extension(NAME + '._pixfun',
                          [NAME + '/pixelfunctions/pixelfunctions.c',
                           NAME + '/pixelfunctions/_pixfun.c'],
                          include_dirs=include_dirs,
                          libraries=libraries,
                          library_dirs=library_dirs,
                          extra_compile_args=extra_compile_args,
                          extra_link_args=extra_link_args)
            ])

    setup(
        name=NAME,
        version=VERSION,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        url=URL,
        download_url=DOWNLOAD_URL,
        license=LICENSE,
        classifiers=CLASSIFIERS,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        platforms=PLATFORMS,
        packages=[NAME, NAME + '.mappers', NAME + '.tests'],
        package_data={NAME: ['wkv.xml', "fonts/*.ttf"]},
        install_requires=REQS,
        **kw
        )

try:
    run_setup(skip_compile)
except ext_errors:
    BUILD_EXT_WARNING = ("WARNING: The C extension could not be compiled, "
                         "pixel functions will not be available.")
    print('*' * 75)
    print(BUILD_EXT_WARNING)
    print("Failure information, if any, is above.")
    print("I'm retrying the build without the C extension now.")
    print('*' * 75)

    run_setup(True)

    print('*' * 75)
    print(BUILD_EXT_WARNING)
    print("Plain-Python installation succeeded.")
    print('*' * 75)
