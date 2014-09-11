#!/usr/bin/env python
#
# Analogue to gdalinfo, but for Nansat datasets 
# Refers to Nansat band numbers

import sys
from nansat import Nansat

def main():
    if (len(sys.argv) != 2):
        sys.exit('Usage: nansatinfo <filename>')
    
    n = Nansat(sys.argv[1])
    print n
    
if __name__ == '__main__':
    main()