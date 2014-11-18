import scipy.ndimage as nd
import matplotlib.pyplot as plt

import numpy as np
from nansat import Nansat, Domain

###############################################################################
####   Increase AQUARIUS data product resolution from 1 deg to 0.1 deg   ######
###############################################################################

# Download Aquarius data and unzip
# wget http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/Q20113052013334.L3m_MC_SCISM_V3.0_SSS_bias_adj_1deg -O AQUARIUS.hdf.bz2
# bunzip2 AQUARIUS.hdf.bz2

# Open Aquarius file with Nansat
nOrig = Nansat('AQUARIUS.hdf', mapperName='obpg_l3')

# Get SSS from Nansat
sssOrig = nOrig['SSS']

# Mask invalid vales with nan
sssOrig[sssOrig < 0] = np.nan

# Fill gaps: extrapolate valid values into 2 pixels border using nearest neighbour
#     get distance and indices of nearest neighbours
dst, ind = nd.distance_transform_edt(np.isnan(sssOrig),
                                     return_distances=True,
                                     return_indices=True)
#     erase row,col indeces further than 2 pixels
ind[0][dst > 2] = 0
ind[1][dst > 2] = 0
#    fill gaps
sssExtra = sssOrig[tuple(ind)]

# Create a Nansat object from matrix with extrapolated product
nExtra = Nansat(domain=nOrig, array=sssExtra)

# Increase resolution
nExtra.resize(10, eResampleAlg=1)

# Get upscaled SSS
sssUpscaled = nExtra[1]

# Plot all SSS for comparison
for sss, sssName in [(sssOrig, 'sss0.png'),
                     (sssExtra, 'sss1.png'),
                     (sssUpscaled, 'sss2.png')]:
    f = plt.figure(figsize=(10,5))
    plt.imshow(sss, vmin=20, vmax=38, interpolation='nearest')
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.savefig(sssName, bbox_inches='tight', pad_inches=0.0, dpi=300)
    plt.close('all')
