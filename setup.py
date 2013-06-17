#-------------------------------------------------------------------------------
# Name:        setup.py
# Purpose:
#
# Author:      asumak
#
# Created:     17.06.2013
# Copyright:   (c) asumak 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os
from subprocess import Popen
import subprocess

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

#myDir = os.environ.get("HOME")
#myNansatDir = str(os).split(" ")[-1][1:].replace("os.pyc'>","site-packages\\"+NAME)
myNansatDir = os.getcwd()

#------------------------------------------------------------------------------#
#                       Set environment variables
#------------------------------------------------------------------------------#
dicDir = {'PYTHONPATH':"\\mappers", 'GDAL_DRIVER_PATH':"\\pixelfunctions"}
'''
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
        print ''
        print iKey, ' : does not exist'
    # if iKey exist
    else:
        # if the folder is not registered yet
        if stdout_value.find(myNansatDir + dicDir[iKey]) == -1:
            oldPath = stdout_value.replace("/","\\").rstrip().split('=')[1]
            newPath = (myNansatDir + dicDir[iKey]).replace("/","\\")
            # command to replace oldfolder to (oldfolder+addFolder)
            command = ('setx %s %s;%s'% (iKey, oldPath, newPath))
            print ''
            print 'oldpath : ', oldPath
            print iKey, ' : not resgistered'
        else:
            command = ''
    if command != '':
        print command
        #process = Popen(command, cwd=path, shell=True, stdout=subprocess.PIPE)
        process = Popen(command, shell=True, stdout=subprocess.PIPE)
        process.stdout.close()
'''
'''
path = "/Home/asumak/package/"
#path = myDir

dicDir = {'PYTHONPATH':"/mappers", 'GDAL_DRIVER_PATH':"/pixelfunctions"}
for iKey in dicDir.keys():
    # check if iKey (environment variable) exist
    command = ("grep '%s' .bashrc" %iKey)
    val = Popen(command, cwd=path, shell=True, stdout=subprocess.PIPE)
    stdout_value = val.communicate()[0]
    # if iKey does not exists
    if stdout_value =="":
        # command to add new environment variable
        command = ("echo 'export %s = %s%s' >> .bashrc" % ( iKey, myNansatDir, dicDir[iKey]))
    # if iKey exist
    else:
        # if the folder is not registered yet
        if stdout_value.find(myNansatDir + dicDir[iKey]) == -1:
            oldFolder = stdout_value.replace("/","\/").rstrip()
            addFolder = (myNansatDir + dicDir[iKey]).replace("/","\/")
            # command to replace oldfolder to (oldfolder+addFolder)
            command = ('sed -i -e "s/%s/%s:%s/g" .bashrc' % (oldFolder, oldFolder, addFolder))
        else:
            command = ""
    if command != "":
        Popen(command, cwd=path, shell=True)
'''
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
