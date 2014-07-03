import os
import inspect


def testio():
    # input and output pathes
    iPath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    oPath = os.path.join(iPath, 'tmpdata')
    if not os.path.exists(oPath):
        os.mkdir(oPath)

    return iPath, oPath
