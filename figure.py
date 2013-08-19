# Name:    figure.py
# Purpose: Container of Figure class
# Authors:      Asuka Yamakawa, Anton Korosov, Knut-Frode Dagestad,
#               Morten W. Hansen, Alexander Myasoyedov,
#               Dmitry Petrenko, Evgeny Morozov
# Created:      29.06.2011
# Copyright:    (c) NERSC 2011 - 2013
# Licence:
# This file is part of NANSAT.
# NANSAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

from nansat_tools import *

class Figure():
    '''Perform opeartions with graphical files: create, append legend, save.

    Figure instance is created in the Nansat.write_figure method
    The methods below are applied consequently in order to generate a figure
    from one or three bands, estimate min/max, apply logarithmic scaling,
    convert to uint8, append legend, save to a file
    '''

    # default values of ALL params of Figure
    cmin = [0.]
    cmax = [1.]
    gamma = 2.
    subsetArraySize = 100000
    numOfColor = 250
    cmapName = 'jet'
    ratio = 1.0
    numOfTicks = 5
    titleString = ''
    caption = ''
    fontSize = 12
    logarithm = False
    legend = False
    mask_array = None
    mask_lut = None

    logoFileName = None
    logoLocation = [0, 0]
    logoSize = None

    latGrid = None
    lonGrid = None
    nGridLines = 10
    latlonLabels = 0

    transparency = None

    LEGEND_HEIGHT = 0.1
    CBAR_HEIGHTMIN = 5
    CBAR_HEIGHT = 0.15
    CBAR_WIDTH = 0.8
    CBAR_LOCATION_X = 0.1
    CBAR_LOCATION_Y = 0.5
    CBAR_LOCATION_ADJUST_X = 5
    CBAR_LOCATION_ADJUST_Y = 3
    TEXT_LOCATION_X = 0.1
    TEXT_LOCATION_Y = 0.1
    NAME_LOCATION_X = 0.1
    NAME_LOCATION_Y = 0.3
    DEFAULT_EXTENSION = '.png'

    palette = None
    pilImg = None
    pilImgLegend = None

    extensionList = ['png', 'PNG', 'tif', 'TIF', 'bmp',
                          'BMP', 'jpg', 'JPG', 'jpeg', 'JPEG']

    _cmapName = 'jet'
        
    def __init__(self, nparray, **kwargs):
        ''' Set attributes

        Parameters
        -----------
        array : numpy array (2D or 3D)
            dataset from Nansat
        kwargs : dictionary
            parameters that are used for all operations.
            See Nansat.write_figure()

        Modifies
        ---------
        self.sizeX, self.sizeY : int
            width and height of the image
        self.pilImg : PIL image
            figure
        self.pilImgLegend : PIL image
            if pilImgLegend is None, legend is not added to the figure
            if it is replaced, pilImgLegend includes text string, color-bar,
            longName and units.

        '''
        from nansat_tools import add_logger

        # make a copy of nparray (otherwise a new reference to the same data is
        # created and the original input data is destroyed at process())
        array = np.array(nparray)

        self.logger = add_logger('Nansat')

        # if 2D array is given, reshape to 3D
        if array.ndim == 2:
            self.array = array.reshape(1, array.shape[0], array.shape[1])
        else:
            self.array = array

        # note swaping of axis by PIL
        self.width = self.array.shape[2]
        self.height = self.array.shape[1]

        # modify the default values using input values
        self._set_defaults(kwargs)


        # set fonts for Legend
        self.fontFileName = os.path.join(os.path.dirname(
                                         os.path.realpath(__file__)),
                                         'fonts/DejaVuSans.ttf')

    def apply_logarithm(self, **kwargs):
        '''Apply a tone curve to the array

        After the normalization of the values from 0 to 1, logarithm is applied
        Then the values are converted to the normal scale.

        Parameters
        -----------
        Any of Figure__init__() parameters

        Modifies
        ---------
        self.array : numpy array

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        # apply logarithm/gamme correction to pixel values
        for iBand in range(self.array.shape[0]):
            self.array[iBand, :, :] = (
                (np.power((self.array[iBand, :, :] - self.cmin[iBand]) /
                         (self.cmax[iBand] - self.cmin[iBand]),
                          (1.0 / self.gamma))) *
                (self.cmax[iBand] - self.cmin[iBand]) +
                self.cmin[iBand])

    def apply_mask(self, **kwargs):
        '''Apply mask for coloring land, clouds, etc

        If mask_array and mask_lut are provided as input parameters
        The pixels in self.array which have index equal to mask_lut kay
        in mask_array will have color equal to mask_lut value

        apply_mask should be called only after convert_palettesize
        (i.e. to uint8 data)

        Parameters
        -----------
        Any of Figure__init__() parameters

        Modifies
        ---------
        self.array : numpy array

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        # get values of free indeces in the palette
        availIndeces = range(self.numOfColor, 255 - 1)

        # for all lut color indeces
        for i, maskValue in enumerate(self.mask_lut):
            if i < len(availIndeces):
                # get color for that index
                maskColor = self.mask_lut[maskValue]
                # get indeces for that index
                maskIndeces = self.mask_array == maskValue
                # exchange colors
                if self.array.shape[0] == 1:
                    # in a indexed image
                    self.array[0][maskIndeces] = availIndeces[i]
                elif self.array.shape[0] == 3:
                    # in RGB image
                    for c in range(0, 3):
                        self.array[c][maskIndeces] = maskColor[c]

                # exchage palette
                self.palette[(availIndeces[i] * 3):
                             (availIndeces[i] * 3 + 3)] = maskColor

    def add_logo(self, **kwargs):
        '''Insert logo into the PIL image

        Read logo from file as PIL
        Resize to the given size
        Pan using the given location
        Paste into pilImg

        Parameters
        ----------
        Any of Figure__init__() parameters

        Modifies
        ---------
        self.pilImg

        '''
        # set/get default parameters
        self._set_defaults(kwargs)
        logoFileName = self.logoFileName
        logoLocation = self.logoLocation
        logoSize = self.logoSize

        # check if pilImg was created already
        if self.pilImg is None:
            self.logger.warning('Create PIL image first')
            return
        # check if file is available
        try:
            logoImg = Image.open(logoFileName)
        except:
            self.logger.warning('No logo file %s' % logoFileName)
            return
        # resize if required
        if logoSize is None:
            logoSize = logoImg.size
        else:
            logoImg = logoImg.resize(logoSize)
        # get location of the logo w.r.t. sign of logoLocation
        box = [0, 0, logoSize[0], logoSize[1]]
        for dim in range(2):
            if logoLocation[dim] >= 0:
                box[dim + 0] = box[dim + 0] + logoLocation[dim + 0]
                box[dim + 2] = box[dim + 2] + logoLocation[dim + 0]
            else:
                box[dim + 0] = (self.pilImg.size[dim + 0] +
                                logoLocation[dim + 0] -
                                logoSize[dim + 0])
                box[dim + 2] = (self.pilImg.size[dim + 0] +
                                logoLocation[dim + 0])

        self.pilImg = self.pilImg.convert('RGB')
        self.pilImg.paste(logoImg, tuple(box))

    def add_latlon_grids(self, **kwargs):
        '''Add lat/lon grid lines into the PIL image

        Compute step of the grid
        Make matrices with binarized lat/lon
        Find edge (make line)
        Convert to maks
        Add mask to PIL

        Parameters
        ----------
        Any of Figure__init__() parameters:
        latGrid : numpy array
            array with values of latitudes
        lonGrid : numpy array
            array with values of longitudes
        nGridLines : int
            number of lines to draw

        Modifies
        ---------
        self.pilImg

        '''
        # modify default values
        self._set_defaults(kwargs)
        # test availability of grids
        if (self.latGrid is None or
            self.lonGrid is None or
            self.nGridLines is None or
            self.nGridLines == 0):
            return
        # get number of grid lines
        llSpacing = self.nGridLines
        # get vectors for grid lines
        latVec = np.linspace(self.latGrid.min(),
                             self.latGrid.max(), llSpacing)
        lonVec = np.linspace(self.lonGrid.min(),
                             self.lonGrid.max(), llSpacing)
        latI = np.zeros(self.latGrid.shape, 'int8')
        lonI = np.zeros(self.latGrid.shape, 'int8')
        # convert lat/lon to indeces
        for i in range(len(latVec)):
            latI[self.latGrid > latVec[i]] = i
            lonI[self.lonGrid > lonVec[i]] = i
        # find pixels on the rgid lines (binarize)
        latI = np.diff(latI)
        lonI = np.diff(lonI)
        # make grid from both lat and lon
        latI += lonI
        latI[latI != 0] = 1
        # add mask to the image
        self.apply_mask(mask_array=latI, mask_lut={1: [255, 255, 255]})

    def add_latlon_labels(self, **kwargs):
        '''Add lat/lon labels along upper and left side

        Compute step of lables
        Get lat/lon for these labels from latGrid, lonGrid
        Print lables to PIL

        Parameters
        ----------
        Figure__init__() parameters:
        latGrid : numpy array
        lonGrid : numpy array
        latlonLabels : int

        Modifies
        ---------
        self.pilImg

        '''
        # modify default values
        self._set_defaults(kwargs)
        # test availability of grids
        if (self.latGrid is None or
            self.lonGrid is None or
            self.latlonLabels == 0):
            return

        draw = ImageDraw.Draw(self.pilImg)
        font = ImageFont.truetype(self.fontFileName, self.fontSize)

        # get number of labels; step of lables
        llLabels = self.latlonLabels
        llShape = self.latGrid.shape
        latI = range(0, llShape[0], (llShape[0] / llLabels) - 1)
        lonI = range(0, llShape[1], (llShape[1] / llLabels) - 1)
        # get lons/lats from first row/column
        #lats = self.latGrid[latI, 0]
        #lons = self.lonGrid[0, lonI]
        for i in range(len(latI)):
            lat = self.latGrid[latI[i], 0]
            lon = self.lonGrid[0, lonI[i]]
            draw.text((0, 10 + latI[i]), '%4.2f' % lat, fill=255, font=font)
            draw.text((50 + lonI[i], 0), '%4.2f' % lon, fill=255, font=font)

    def clim_from_histogram(self, **kwargs):
        '''Estimate min and max pixel values from histogram

        if ratio=1.0, simply the minimum and maximum values are returned.
        if 0 < ratio < 1.0, get the histogram of the pixel values.
        Then get rid of (1.0-ratio)/2 from the both sides and
        return the minimum and maximum values.

        Parameters
        -----------
        Any of Figure.__init__() parameters

        Returns
        --------
        clim : numpy array 2D ((3x2) or (1x2))
            minimum and maximum pixel values for each band

        '''
        # modify default values
        self._set_defaults(kwargs)
        ratio = self.ratio

        # create a ratio list for each band
        if isinstance(ratio, float) or isinstance(ratio, int):
            ratioList = np.ones(self.array.shape[0]) * float(ratio)
        else:
            ratioList = []
            for iRatio in range(self.array.shape[0]):
                try:
                    ratioList.append(ratio[iRatio])
                except:
                    ratioList.append(ratio[0])

        # create a 2D array and set min and max values
        clim = [[0] * self.array.shape[0], [0] * self.array.shape[0]]
        for iBand in range(self.array.shape[0]):
            clim[0][iBand] = self.array[iBand, :, :].min()
            clim[1][iBand] = self.array[iBand, :, :].max()
            # if 0<ratio<1 try to compute histogram
            if (ratioList[iBand] > 0 or ratioList[iBand] < 1):
                try:
                    hist, bins = self._get_histogram(iBand)
                except:
                    self.logger.warning('Unable to compute histogram')
                else:
                    cumhist = hist.cumsum()
                    cumhist /= cumhist[-1]
                    clim[0][iBand] = bins[len(cumhist[cumhist <
                                              (1 - ratioList[iBand]) / 2])]
                    clim[1][iBand] = bins[len(cumhist[cumhist <
                                              1 - ((1 - ratioList[iBand]) /
                                                   2)])]
        self.color_limits = clim
        return clim

    def clip(self, **kwargs):
        '''Convert self.array to values between cmin and cmax

        if pixel value < cmin, replaced to cmin.
        if pixel value > cmax, replaced to cmax.

        Parameters
        -----------
        Any of Figure.__init__() parameters

        Modifies
        ---------
        self.array : numpy array
        self.cmin, self.cmax : allowed min/max values

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        for iBand in range(self.array.shape[0]):
            # if clipping integer matrix, make clipping ranges valid
            if self.array.dtype in ['int8', 'uint8', 'int16', 'uint16']:
                self.cmin[iBand] = np.ceil(self.cmin[iBand])
                self.cmin[iBand] = np.floor(self.cmin[iBand])

            # Clipping, allowing for reversed colorscale (cmin > cmax)
            clipMin = np.min([self.cmin[iBand], self.cmax[iBand]])
            clipMax = np.max([self.cmin[iBand], self.cmax[iBand]])
            self.array[iBand, :, :] = np.clip(self.array[iBand, :, :],
                                              clipMin, clipMax)

    def convert_palettesize(self, **kwargs):
        '''Convert self.array to palette color size in uint8

        Parameters
        -----------

        Any of Figure.__init__() parameters

        Modifies
        ---------
        self.array : numpy array (=>uint8)

        '''
        # modify default values
        self._set_defaults(kwargs)

        for iBand in range(self.array.shape[0]):
            self.array[iBand, :, :] = (
                (self.array[iBand, :, :].astype('float32') -
                 self.cmin[iBand]) *
                (self.numOfColor - 1) /
                (self.cmax[iBand] - self.cmin[iBand]))

        self.array = self.array.astype(np.uint8)

    def create_legend(self, **kwargs):
        ''' self.legend is replaced from None to PIL image

        PIL image includes colorbar, caption, and titleString.

        Parameters
        -----------
        Any of Figure.__init__() parameters

        Modifies
        ---------
        self.legend : PIL image

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        # set fonts size for colorbar
        font = ImageFont.truetype(self.fontFileName, self.fontSize)

        # create a pilImage for the legend
        self.pilImgLegend = Image.new('P', (self.width,
                                      int(self.height *
                                      self.LEGEND_HEIGHT)), 255)
        draw = ImageDraw.Draw(self.pilImgLegend)

        # set black color
        if self.array.shape[0] == 1:
            black = 254
        else:
            black = (0, 0, 0)

        # if 1 band, draw the color bar
        if self.array.shape[0] == 1:
            # make an array for color bar
            bar = np.outer(np.ones(max(int(self.pilImgLegend.size[1] *
                           self.CBAR_HEIGHT), self.CBAR_HEIGHTMIN)),
                           np.linspace(0, self.numOfColor,
                           int(self.pilImgLegend.size[0] *
                           self.CBAR_WIDTH)))
            # create a colorbar pil Image
            pilImgCbar = Image.fromarray(np.uint8(bar))
            # paste the colorbar pilImage on Legend pilImage
            self.pilImgLegend.paste(pilImgCbar,
                                    (int(self.pilImgLegend.size[0] *
                                     self.CBAR_LOCATION_X),
                                     int(self.pilImgLegend.size[1] *
                                     self.CBAR_LOCATION_Y)))
            # create a scale for the colorbar
            scaleLocation = np.linspace(0, 1, self.numOfTicks)
            scaleArray = scaleLocation
            if self.logarithm:
                scaleArray = (np.power(scaleArray, (1.0 / self.gamma)))
            scaleArray = (scaleArray * (self.cmax[0] -
                          self.cmin[0]) + self.cmin[0])
            scaleArray = map(self._round_number, scaleArray)
            # draw scales and lines on the legend pilImage
            for iTick in range(self.numOfTicks):
                coordX = int(scaleLocation[iTick] *
                             self.pilImgLegend.size[0] *
                             self.CBAR_WIDTH +
                             int(self.pilImgLegend.size[0] *
                             self.CBAR_LOCATION_X))

                box = (coordX, int(self.pilImgLegend.size[1] *
                                   self.CBAR_LOCATION_Y),
                       coordX, int(self.pilImgLegend.size[1] *
                                  (self.CBAR_LOCATION_Y +
                                   self.CBAR_HEIGHT)) - 1)
                draw.line(box, fill=black)
                box = (coordX - self.CBAR_LOCATION_ADJUST_X,
                       int(self.pilImgLegend.size[1] *
                           (self.CBAR_LOCATION_Y +
                            self.CBAR_HEIGHT)) +
                       self.CBAR_LOCATION_ADJUST_Y)
                draw.text(box, scaleArray[iTick], fill=black, font=font)

        # draw longname and units
        box = (int(self.pilImgLegend.size[0] * self.NAME_LOCATION_X),
               int(self.pilImgLegend.size[1] * self.NAME_LOCATION_Y))
        draw.text(box, str(self.caption), fill=black, font=font)

        # if titleString is given, draw it
        if self.titleString != '':
            # write text each line onto pilImgCanvas
            textHeight = int(self.pilImgLegend.size[1] *
                             self.TEXT_LOCATION_Y)
            for line in self.titleString.splitlines():
                draw.text((int(self.pilImgLegend.size[0] *
                               self.TEXT_LOCATION_X),
                           textHeight), line, fill=black, font=font)
                text = draw.textsize(line, font=font)
                textHeight += text[1]

    def create_pilImage(self, **kwargs):
        ''' self.create_pilImage is replaced from None to PIL image

        If three images are given, create a image with RGB mode.
            if self.pilImgLegend is not None, it is pasted.
        If one image is given, create a image with P(palette) mode.
            if self.pilImgLegend is not None,
            self.array is extended before create the pilImag and
            then paste pilImgLegend onto it.

        Parameters
        -----------
        Any of Figure.__init__() parameters

        Modifies
        ---------
        self.pilImg : PIL image
            PIL image with / without the legend
        self.array : replace to None

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        # if legend is created, expand array with empty space below the data
        if self.pilImgLegend is not None:
            appendArray = 255 * np.ones((self.array.shape[0],
                                         self.pilImgLegend.size[1],
                                         self.width), 'uint8')
            self.array = np.append(self.array, appendArray, 1)

        # create a new PIL image from three bands (RGB) or from one (palette)
        if self.array.shape[0] == 3:
            self.pilImg = Image.merge('RGB',
                                      (Image.fromarray(self.array[0, :, :]),
                                       Image.fromarray(self.array[1, :, :]),
                                       Image.fromarray(self.array[2, :, :])))
        else:
            self.pilImg = Image.fromarray(self.array[0, :, :])
            self.pilImg.putpalette(self.palette)

        # append legend
        if self.pilImgLegend is not None:
            self.pilImg.paste(self.pilImgLegend, (0, self.height))

        # remove array from memory
        #self.array = None

    def process(self, **kwargs):
        '''Do all common operations for preparation of a figure for saving

        #. Modify default values of parameters by the provided ones (if any)
        #. Clip to min/max
        #. Apply logarithm if required
        #. Convert data to uint8
        #. Create palette
        #. Apply mask for colouring land, clouds, etc if required
        #. Create legend if required
        #. Create PIL image
        #. Add logo if required

        Parameters
        -----------
        Any of Figure.__init__() parameters

        Modifies
        --------
        self.d
        self.array
        self.palette
        self.pilImgLegend
        self.pilImg

        '''
        # modify default parameters
        self._set_defaults(kwargs)

        # if the image is reprojected it has 0 values
        # we replace them with mask before creating PIL Image
        self.reprojMask = self.array[0, :, :] == 0

        # clip values to min/max
        self.clip()

        # apply logarithm
        if self.logarithm:
            self.apply_logarithm()

        # convert to uint8
        self.convert_palettesize()

        # create the paletter
        self._create_palette()

        # apply colored mask (land mask, cloud mask and something else)
        if self.mask_array is not None and self.mask_lut is not None:
            self.apply_mask()

        # add lat/lon grids lines if latGrid and lonGrid are given
        if self.latGrid is not None and self.lonGrid is not None:
            self.add_latlon_grids()

        # append legend
        if self.legend:
            self.create_legend()

        # create PIL image ready for saving
        self.create_pilImage(**kwargs)

        # add labels with lats/lons
        if (self.latGrid is not None and
            self.lonGrid is not None and
            self.latlonLabels > 0):
            self.add_latlon_labels()

        # add logo
        if self.logoFileName is not None:
            self.add_logo()

    def _make_transparent_color(self):
        ''' makes colors specified by self.transparency
        and self.reprojMask (if the image is reprojected) transparent

        Modifies
        --------
        self.pilImg : PIL image
            Adds transparency to PIL image

        '''
        self.pilImg = self.pilImg.convert('RGBA')
        datas = self.pilImg.getdata()
        newData = list()

        for item in datas:
            if (item[0] == self.transparency[0] and
                item[1] == self.transparency[1] and
                item[2] == self.transparency[2]):
                newData.append((255, 255, 255, 0))
            else:
                newData.append(item)

        self.pilImg.putdata(newData)

        # The alphaMask is set in process() before clip() the Image
        img = np.array(self.pilImg)
        img[:, :, 3][self.reprojMask] = 0
        self.pilImg = Image.fromarray(np.uint8(img))

    def save(self, fileName, **kwargs):
        ''' Save self.pilImg to a physical file

        If given extension is JPG, convert the image mode from Palette to RGB

        Parameters
        ----------
        fileName : string
            name of outputfile
        Any of Figure.__init__() parameters

        Modifies
        --------
        self.pilImg : None

        '''
        # modify default values
        self._set_defaults(kwargs)

        if not((fileName.split('.')[-1] in self.extensionList)):
            fileName = fileName + self.DEFAULT_EXTENSION

        fileExtension = fileName.split('.')[-1]
        if fileExtension in ['jpg', 'JPG', 'jpeg', 'JPEG']:
            self.pilImg = self.pilImg.convert('RGB')

        if self.transparency is not None:
            self._make_transparent_color()
        self.pilImg.save(fileName)

    def _create_palette(self):
        '''Create a palette based on Matplotlib colormap name

        default number of color palette is 250.
        it means 6 colors are possible to use for other purposes.
        the last palette (255) is white and the second last (254) is black.

        Modifies
        --------
        self.palette : numpy array (uint8)

        '''
        # create a colorList from matplotlib colormap name
        try:
            colorDic = cm.datad[self.cmapName]
        except:
            colorDic = cm.datad[self._cmapName]
            self.cmapName = self._cmapName

        colorList = [colorDic['red'], colorDic['green'], colorDic['blue']]

        # create a numpyarray for a palette based on the color list
        # default is all values are black (=0)
        lut = np.zeros([3, 256])
        # place colors to each number (palette)
        for iColor in range(3):
            iPalette = 0
            for i in range(len(colorList[iColor]) - 1):
                spaceNum = int(self.numOfColor *
                               (colorList[iColor][i + 1][0] -
                               colorList[iColor][i][0]))
                lut[iColor][iPalette:iPalette + spaceNum] = np.array(
                                                np.linspace(
                                                colorList[iColor][i][2],
                                                colorList[iColor][i + 1][1],
                                                num=spaceNum) * 255,
                                                dtype=np.uint8)
                iPalette += (spaceNum)
            # adjust the number of colors on the palette
            while iPalette < self.numOfColor:
                lut[iColor][iPalette] = lut[iColor][iPalette - 1]
                iPalette += 1
            while iPalette > self.numOfColor:
                lut[iColor][iPalette] = 0
                iPalette -= 1
        # the last palette color is replaced to white
        for iColor in range(3):
            lut[iColor][-1] = 255

        # set palette
        self.palette = lut.T.flatten().astype(np.uint8)

    def _get_histogram(self, iBand):
        '''Create a subset array and return the histogram.

        Parameters
        -----------
        iBand : int

        Returns
        --------
        hist : numpy array
        bins : numpy array

        '''
        array = self.array[iBand, :, :].flatten()
        array = array[array > array.min()]
        array = array[array < array.max()]
        step = max(int(round(float(len(array)) /
                       float(self.subsetArraySize))), 1.0)
        arraySubset = array[::step]
        hist, bins, patches = plt.hist(arraySubset, bins=100)
        plt.close()
        return hist.astype(float), bins

    def _round_number(self, val):
        '''Return writing format for scale on the colorbar

        Parameters
        ----------
        val : int / float / exponential

        Returns
        --------
        string

        '''
        frmts = {-2: '%.2f', -1: '%.1f', 0: '%.2f',
                 1: '%.1f', 2: '%d', 3: '%d'}
        if val == 0:
            frmt = '%d'
        else:
            digit = floor(log10(abs(val)))
            if digit in frmts:
                frmt = frmts[digit]
            else:
                frmt = '%4.2e'

        return str(frmt % val)

    def _set_defaults(self, idict):
        '''Check input params and set defaut values

        Look throught default parameters (self.d) and given parameters (dict)
        and paste value from input if the key matches

        Parameters
        ----------
        idict : dictionary
            parameter names and values

        Modifies
        ---------
            default self attributes

        '''
        for key in idict:
            if hasattr(self, key):
                if key in ['cmin', 'cmax'] and type(idict[key]) != list:
                    setattr(self, key, [idict[key]])
                else:
                    setattr(self, key, idict[key])

