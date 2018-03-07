#-----------------------------------------------------------------------------
# Name:        setup.py
# Purpose:
#
# Author:      asumak
#
# Created:     17.06.2013
# Copyright:   (c) asumak 2012
# Licence:     <your licence>
#-----------------------------------------------------------------------------

from __future__ import print_function
from subprocess import Popen
import subprocess
import sys
import errno
import os

import_error_msg = "Nansat requires %s, which should be installed separately"

# Check if required packages are installed
try:
    import numpy
except ImportError:
    raise ImportError(import_error_msg % 'numpy')

try:
    from osgeo import gdal, osr, ogr
except ImportError:
    try:
        import gdal
        import osr
        import ogr
    except ImportError:
        raise ImportError(import_error_msg % 'gdal')

NAME                = 'nansat'
MAINTAINER          = "Nansat Developers"
MAINTAINER_EMAIL    = "nansat-dev@googlegroups.com"
DESCRIPTION         = "A scientist friendly Python toolbox for processing 2D satellite Earth observation data"
LONG_DESCRIPTION    = "A scientist friendly Python toolbox for processing 2D satellite Earth observation data"
URL                 = "https://github.com/nansencenter/nansat"
DOWNLOAD_URL        = "https://github.com/nansencenter/nansat"
LICENSE             = "GNU General Public License"
CLASSIFIERS         = [
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Utilities'
    ]
AUTHOR              = ("Anton Korosov, Morten W. Hansen, Kunt-Frode Dagestad, Aleksander Vines, Asuka Yamakawa")
AUTHOR_EMAIL        = "nansat-dev@googlegroups.com"
PLATFORMS           = ["Linux", "OS X", "Windows"]
MAJOR               = 1
MINOR               = 0
MICRO               = 11
ISRELEASED          = True
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO) # Remember to remove "dev" when releasing
REQS                = [
                        "Pillow",
                        "pythesint",
                        "cfunits",
                        "urllib3",
                    ]

#----------------------------------------------------------------------------#
#                   Prepare compilation of C pixel functions
#----------------------------------------------------------------------------#
pixfun_module_name = '_pixfun_py{0}'.format(sys.version_info[0])
skip_compile = False
libraries = []
include_dirs = []
library_dirs = []
if os.name == 'nt':
    extra_compile_args = ['-nologo', '-DLL']
    try:
        path = os.environ['CONDA_PREFIX'] + '\Library\lib'
    except:
        # TODO: figure out how this works if you are not using conda
        path = os.environ['LIB'].split(';')
    #for iFolder in path:
    try:
        files = os.listdir(path)
    except:
        extra_link_args = []
    else:
        if 'gdal_i.lib' in files:
            extra_link_args = [path + '\gdal_i.lib']
else:
    extra_compile_args = ['-fPIC', '-Wall', '-Wno-long-long', '-pedantic', '-O3']
    extra_link_args = [] # not used currently


def _ask_gdal_config(resultlist, option, result_prefix):
    p = Popen(['gdal-config', option], stdout=subprocess.PIPE)
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

def use_win_config():
    #try:
    path = os.environ['CONDA_PREFIX']
    include_dirs[:] = [path+'\Library\include']
    #except:
        # TODO: figure out how this works if you are not using conda
        # path = os.environ['LIB'].split(';')

try:
    if os.name == 'nt':
        use_win_config()
    else:
        use_gdal_config()
except Exception as e:
    if os.name == 'nt':
        print('WARNING: CONDA_PREFIX could not be found pixel functions will not be available.')
    else:
        print('WARNING: gdal-config could not be called, pixel functions will not be available.')
    print('Error details follow:')
    print(e)
    skip_compile = True

#----------------------------------------------------------------------------#
#                               Install package
#----------------------------------------------------------------------------#
from setuptools import setup, find_packages
from setuptools.command.install_scripts import install_scripts
from distutils import log
from distutils.extension import Extension
from distutils.errors import CCompilerError, DistutilsExecError,\
    DistutilsPlatformError

# Windows batch file handling
# from https://matthew-brett.github.io/pydagogue/installing_scripts.html
BAT_TEMPLATE = \
r"""@echo off
set mypath=%~dp0
set pyscript="%mypath%{FNAME}"
set /p line1=<%pyscript%
if "%line1:~0,2%" == "#!" (goto :goodstart)
echo First line of %pyscript% does not start with "#!"
exit /b 1
:goodstart
set py_exe=%line1:~2%
call %py_exe% %pyscript% %*
"""


class my_install_scripts(install_scripts):
    def run(self):
        install_scripts.run(self)
        if not os.name == "nt":
            return
        for filepath in self.get_outputs():
            # If we can find an executable name in the #! top line of the script
            # file, make .bat wrapper for script.
            with open(filepath, 'rt') as fobj:
                first_line = fobj.readline()
            if not (first_line.startswith('#!') and
                    'python' in first_line.lower()):
                log.info("No #!python executable found, skipping .bat "
                            "wrapper")
                continue
            pth, fname = os.path.split(filepath)
            froot, _ = os.path.splitext(fname)
            bat_file = os.path.join(pth, froot + '.bat')
            bat_contents = BAT_TEMPLATE.replace('{FNAME}', fname)
            log.info("Making %s wrapper for %s" % (bat_file, filepath))
            if self.dry_run:
                continue
            with open(bat_file, 'wt') as fobj:
                fobj.write(bat_contents)

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
                Extension('{0}.{1}'.format(NAME, pixfun_module_name),
                          ['{0}/pixelfunctions/pixelfunctions.c'.format(NAME),
                           '{0}/pixelfunctions/{1}.c'.format(NAME, pixfun_module_name)],
                          include_dirs=include_dirs,
                          libraries=libraries,
                          library_dirs=library_dirs,
                          extra_compile_args=extra_compile_args,
                          extra_link_args=extra_link_args)
            ])


    # remove mapper_tests from installed packages
    packages = find_packages()
    #if 'mapper_tests' in packages:
    #    packages.remove('mapper_tests')

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
        packages=packages,
        package_data={NAME:["fonts/*.ttf", 'mappers/*.pl', 'tests/data/*.*']},
        scripts=[os.path.join('utilities', name) for name in
                    ['nansatinfo',
                     'nansat_add_coastline',
                     'nansat_geotiffimage',
                     'nansat_show',
                     'nansat_translate',
                     ]],
        cmdclass = {'install_scripts': my_install_scripts},
        install_requires=REQS,
        test_suite="nansat.tests",
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
