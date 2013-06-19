#-------------------------------------------------------------------------------
# Name:        setup.py
# Purpose:
#
# Author:      asumak
#
# Created:     17.06.2013
# Copyright:   (c) asumak 2012
# Licence:     <your licence>
# =========  !! NB !! HOW TO DO FOR MAC USERS??  ==========
#-------------------------------------------------------------------------------

import os
from subprocess import Popen
import subprocess
import sys

NAME                = 'nansat'
MAINTAINER          = "Nansat Developers"
MAINTAINER_EMAIL    = "numpy-discussion@scipy.org"
DESCRIPTION         = "***" #DOCLINES[0]
LONG_DESCRIPTION    = "***" #"\n".join(DOCLINES[2:])
URL                 = "http://normap.nersc.no/"
DOWNLOAD_URL        = "http://normap.nersc.no/"
LICENSE             = '***'
CLASSIFIERS         = '***' #filter(None, CLASSIFIERS.split('\n'))
AUTHOR              = "Asuka Yamakawa, Anton Korosov, Morten W. Hansen, Kunt-Frode Dagestad"
AUTHOR_EMAIL        = "asuka.yamakawa@nersc.no"
PLATFORMS           = ["UNKNOWN"]
MAJOR               = 1
MINOR               = 0
MICRO               = 0
ISRELEASED          = True
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

osName = sys.platform
myNansatDir = os.getcwd()
""" !! NB: How to do for Mac ??  """
if not('win' in osName):
    myHomeDir = os.environ.get("HOME")
    index = myNansatDir.rfind(myHomeDir)
    myNansatDir = myNansatDir[index:]


#------------------------------------------------------------------------------#
#                       Set environment variables
#------------------------------------------------------------------------------#
if 'win' in osName:
    dicDir = {'PYTHONPATH':"\\mappers", 'GDAL_DRIVER_PATH':"\\pixelfunctions"}
    for iKey in dicDir.keys():
        # check if iKey (environment variable) exist
        command = ("set %s" %iKey)
        val = Popen(command, shell=True, stdout=subprocess.PIPE)
        stdout_value = val.communicate()[0]
        # if iKey does not exists
        if stdout_value == "":
            # command to add new environment variable
            myNansatDir = myNansatDir.replace("/","\\")
            command = ("setx %s %s%s" % ( iKey, myNansatDir, dicDir[iKey]))
        # if iKey exist
        else:
            # if the folder is not registered yet
            if stdout_value.find(myNansatDir + dicDir[iKey]) == -1:
                oldPath = stdout_value.replace("/","\\").rstrip().split('=')[1]
                newPath = (myNansatDir + dicDir[iKey]).replace("/","\\")
                # command to replace oldfolder to (oldfolder+addFolder)
                command = ('setx %s %s;%s'% (iKey, oldPath, newPath))
            else:
                command = ''
        if command != '':
            process = Popen(command, shell=True, stdout=subprocess.PIPE)
            process.stdout.close()
    """ !! NB: How to do for Mac ??  """
else:
    dicDir = {'PYTHONPATH':"/mappers", 'GDAL_DRIVER_PATH':"/pixelfunctions"}
    for iKey in dicDir.keys():
        # check if iKey (environment variable) exist
        command = ("grep '%s' .bashrc" %iKey)
        val = Popen(command, cwd=myHomeDir, shell=True, stdout=subprocess.PIPE)
        stdout_value = val.communicate()[0]
        # if iKey does not exists
        if stdout_value =="":
            # command to add new environment variable
            command = ("echo 'export %s=%s%s' >> .bashrc" % ( iKey, myNansatDir, dicDir[iKey]))
        # if iKey exist
        else:
            # if the folder is not registered yet
            if stdout_value.find(myNansatDir + dicDir[iKey]) == -1:
                command = ('echo "export %s=\$%s:%s" >> .bashrc' % (iKey, iKey, myNansatDir + dicDir[iKey]))
            else:
                command = ""
        if command != "":
            process = Popen(command, cwd=myHomeDir, shell=True, stdout=subprocess.PIPE)
            process.stdout.close()

#------------------------------------------------------------------------------#
#                               Copy files
#------------------------------------------------------------------------------#
from distutils.core import setup
setup(
    name=NAME,
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
    package_dir={NAME : ''},
    packages = {NAME, NAME + '.mappers'},
    package_data={NAME: ['wkv.xml', "fonts/*.ttf", "pixelfunctions/*"]},
    )
