import os
import inspect
def testio():
    # input and output file names
    iPath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    iFileName = os.path.join(iPath, 'gcps.tif')
    print 'Input file: ', iFileName
    oPath = os.path.join(iPath, 'tmpdata')
    print 'Output path:', oPath
    if not os.path.exists(oPath):
        os.mkdir(oPath)
    oFileName = os.path.join(oPath, 'output_')
    print 'Output file:', oFileName

    return iPath, iFileName, oPath, oFileName
