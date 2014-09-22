#!/usr/bin/env python
#
# Utility to quickly view one band of a file supported by Nansat/GDAL

# Could be extended to use more of the features of the Figure-class

import sys
import os
import Image
from nansat import Nansat

def Usage():
    sys.exit('Usage: nansat_show <band_number> <input filename> [<output_file>]')

def main():
    if (len(sys.argv) < 3):
        Usage()
    
    bandNo = int(sys.argv[1])
    try:
        n = Nansat(sys.argv[2])
    except:
        Usage()
    
    try:
        outfile = sys.argv[3]
        delete = False
    except:
        outfile = 'tmp.png'
        delete = True
    
    n.write_figure(outfile, bandNo, legend=True)
    Image.open(outfile).show()
    
    if delete:
        os.remove(outfile)

if __name__ == '__main__':
    main()