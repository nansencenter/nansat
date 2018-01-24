/******************************************************************************
 *
 * Project:  GDAL
 * Purpose:  Implementation of a set of GDALDerivedPixelFunc(s) to be used
 *           with source raster band of virtual GDAL datasets.
 * Author:   Antonio Valentino <a_valentino@users.sf.net>
 *
 ******************************************************************************
 * Copyright (c) 2008-2011 Antonio Valentino
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *****************************************************************************/

#include <math.h>
#include <gdal.h>
#include <stdio.h>
#include <stdlib.h>

void GenericPixelFunction(double f(double*), void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace);

CPLErr RealPixelFunc(void **papoSources, int nSources, void *pData,
                     int nXSize, int nYSize,
                     GDALDataType eSrcType, GDALDataType eBufType,
                     int nPixelSpace, int nLineSpace)
{
    int iLine, nPixelSpaceSrc, nLineSpaceSrc;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;

    nPixelSpaceSrc = GDALGetDataTypeSize( eSrcType ) / 8;
    nLineSpaceSrc = nPixelSpaceSrc * nXSize;

    /* ---- Set pixels ---- */
    for( iLine = 0; iLine < nYSize; ++iLine ) {
        GDALCopyWords(((GByte *)papoSources[0]) + nLineSpaceSrc * iLine,
                      eSrcType, nPixelSpaceSrc,
                      ((GByte *)pData) + nLineSpace * iLine,
                      eBufType, nPixelSpace, nXSize);
    }

    /* ---- Return success ---- */
    return CE_None;
} /* RealPixelFunc */


CPLErr ImagPixelFunc(void **papoSources, int nSources, void *pData,
                     int nXSize, int nYSize,
                     GDALDataType eSrcType, GDALDataType eBufType,
                     int nPixelSpace, int nLineSpace)
{
    int iLine;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;

    if (GDALDataTypeIsComplex( eSrcType ))
    {
        int nPixelSpaceSrc = GDALGetDataTypeSize( eSrcType ) / 8;
        int nLineSpaceSrc = nPixelSpaceSrc * nXSize;

        void* pImag = ((GByte *)papoSources[0])
                    + GDALGetDataTypeSize( eSrcType ) / 8 / 2;

        /* ---- Set pixels ---- */
        for( iLine = 0; iLine < nYSize; ++iLine ) {
            GDALCopyWords(((GByte *)pImag) + nLineSpaceSrc * iLine,
                          eSrcType, nPixelSpaceSrc,
                          ((GByte *)pData) + nLineSpace * iLine,
                          eBufType, nPixelSpace, nXSize);
        }
    } else {
        double dfImag = 0;

        /* ---- Set pixels ---- */
        for( iLine = 0; iLine < nYSize; ++iLine ) {
            /* always copy from the same location */
            GDALCopyWords(&dfImag, eSrcType, 0,
                          ((GByte *)pData) + nLineSpace * iLine,
                          eBufType, nPixelSpace, nXSize);
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* ImagPixelFunc */


CPLErr ModulePixelFunc(void **papoSources, int nSources, void *pData,
                       int nXSize, int nYSize,
                       GDALDataType eSrcType, GDALDataType eBufType,
                       int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;
    double dfPixVal;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;

    if (GDALDataTypeIsComplex( eSrcType ))
    {
        double dfReal, dfImag;
        void *pReal = papoSources[0];
        void *pImag = ((GByte *)papoSources[0])
                    + GDALGetDataTypeSize( eSrcType ) / 8 / 2;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii = 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfReal = SRCVAL(pReal, eSrcType, ii);
                dfImag = SRCVAL(pImag, eSrcType, ii);

                dfPixVal = sqrt( dfReal * dfReal + dfImag * dfImag );

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        /* ---- Set pixels ---- */
        for( iLine = 0, ii = 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfPixVal = abs(SRCVAL(papoSources[0], eSrcType, ii));

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* ModulePixelFunc */


CPLErr PhasePixelFunc(void **papoSources, int nSources, void *pData,
                      int nXSize, int nYSize,
                      GDALDataType eSrcType, GDALDataType eBufType,
                      int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;
    double dfPixVal, dfReal;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;

    if (GDALDataTypeIsComplex( eSrcType ))
    {
        double dfImag;
        void *pReal = papoSources[0];
        void *pImag = ((GByte *)papoSources[0])
                    + GDALGetDataTypeSize( eSrcType ) / 8 / 2;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfReal = SRCVAL(pReal, eSrcType, ii);
                dfImag = SRCVAL(pImag, eSrcType, ii);

                dfPixVal = atan2(dfImag, dfReal);

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine +
                          iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        /* ---- Set pixels ---- */
        /*
        for( iLine = 0; iLine < nYSize; ++iLine ) {
            / * always copy from the same location * /
            GDALCopyWords(&dfImag, eSrcType, 0,
                          ((GByte *)pData) + nLineSpace * iLine,
                          eBufType, nPixelSpace, nXSize);
        }
        */
        /* ---- Set pixels ---- */
        double pi = atan2(0, -1);
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                void *pReal = papoSources[0];

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfReal = SRCVAL(pReal, eSrcType, ii);
                dfPixVal = (dfReal < 0) ? pi : 0;

                GDALCopyWords(&dfPixVal, GDT_Float64, dfPixVal,
                          ((GByte *)pData) + nLineSpace * iLine +
                          iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* PhasePixelFunc */


CPLErr ConjPixelFunc(void **papoSources, int nSources, void *pData,
                     int nXSize, int nYSize,
                     GDALDataType eSrcType, GDALDataType eBufType,
                     int nPixelSpace, int nLineSpace)
{
    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;

    if (GDALDataTypeIsComplex( eSrcType ) && GDALDataTypeIsComplex( eBufType ))
    {
        int iLine, iCol, ii;
        double adfPixVal[2];
        int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;
        void *pReal = papoSources[0];
        void *pImag = ((GByte *)papoSources[0]) + nOffset;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                adfPixVal[0] = +SRCVAL(pReal, eSrcType, ii); /* re */
                adfPixVal[1] = -SRCVAL(pImag, eSrcType, ii); /* im */

                GDALCopyWords(adfPixVal, GDT_CFloat64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        /* no complex data type */
        return RealPixelFunc(papoSources, nSources, pData, nXSize, nYSize,
                             eSrcType, eBufType, nPixelSpace, nLineSpace);
    }

    /* ---- Return success ---- */
    return CE_None;
} /* ConjPixelFunc */


CPLErr SumPixelFunc(void **papoSources, int nSources, void *pData,
                    int nXSize, int nYSize,
                    GDALDataType eSrcType, GDALDataType eBufType,
                    int nPixelSpace, int nLineSpace)
{
    int iLine, iCol, ii, iSrc;

    /* ---- Init ---- */
    if (nSources < 2) return CE_Failure;

    /* ---- Set pixels ---- */
    if (GDALDataTypeIsComplex( eSrcType ))
    {
        double adfSum[2];
        void *pReal, *pImag;
        int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {
                adfSum[0] = 0;
                adfSum[1] = 0;

                for( iSrc = 0; iSrc < nSources; ++iSrc ) {
                    pReal = papoSources[iSrc];
                    pImag = ((GByte *)pReal) + nOffset;

                    /* Source raster pixels may be obtained with SRCVAL macro */
                    adfSum[0] += SRCVAL(pReal, eSrcType, ii);
                    adfSum[1] += SRCVAL(pImag, eSrcType, ii);
                }

                GDALCopyWords(adfSum, GDT_CFloat64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        /* non complex */
        double dfSum;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {
                dfSum = 0;

                for( iSrc = 0; iSrc < nSources; ++iSrc ) {
                    /* Source raster pixels may be obtained with SRCVAL macro */
                    dfSum += SRCVAL(papoSources[iSrc], eSrcType, ii);
                }

                GDALCopyWords(&dfSum, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* SumPixelFunc */


CPLErr DiffPixelFunc(void **papoSources, int nSources, void *pData,
                     int nXSize, int nYSize,
                     GDALDataType eSrcType, GDALDataType eBufType,
                     int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;

    /* ---- Init ---- */
    if (nSources != 2) return CE_Failure;

    if (GDALDataTypeIsComplex( eSrcType ))
    {

        double adfPixVal[2];
        int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;
        void *pReal0 = papoSources[0];
        void *pImag0 = ((GByte *)papoSources[0]) + nOffset;
        void *pReal1 = papoSources[1];
        void *pImag1 = ((GByte *)papoSources[1]) + nOffset;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                adfPixVal[0] = SRCVAL(pReal0, eSrcType, ii)
                             - SRCVAL(pReal1, eSrcType, ii);
                adfPixVal[1] = SRCVAL(pImag0, eSrcType, ii)
                             - SRCVAL(pImag1, eSrcType, ii);

                GDALCopyWords(adfPixVal, GDT_CFloat64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        /* non complex */
        double dfPixVal;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfPixVal = SRCVAL(papoSources[0], eSrcType, ii)
                         - SRCVAL(papoSources[1], eSrcType, ii);

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* DiffPixelFunc */


CPLErr MulPixelFunc(void **papoSources, int nSources, void *pData,
                    int nXSize, int nYSize,
                    GDALDataType eSrcType, GDALDataType eBufType,
                    int nPixelSpace, int nLineSpace)
{
    int iLine, iCol, ii, iSrc;

    /* ---- Init ---- */
    if (nSources < 2) return CE_Failure;

    /* ---- Set pixels ---- */
    if (GDALDataTypeIsComplex( eSrcType ))
    {
        double adfPixVal[2], dfOldR, dfOldI, dfNewR, dfNewI;
        void *pReal, *pImag;
        int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {
                adfPixVal[0] = 1.;
                adfPixVal[1] = 0.;

                for( iSrc = 0; iSrc < nSources; ++iSrc ) {
                    pReal = papoSources[iSrc];
                    pImag = ((GByte *)pReal) + nOffset;

                    dfOldR = adfPixVal[0];
                    dfOldI = adfPixVal[1];

                    /* Source raster pixels may be obtained with SRCVAL macro */
                    dfNewR = SRCVAL(pReal, eSrcType, ii);
                    dfNewI = SRCVAL(pImag, eSrcType, ii);

                    adfPixVal[0] = dfOldR * dfNewR - dfOldI * dfNewI;
                    adfPixVal[1] = dfOldR * dfNewI + dfOldI * dfNewR;
                }

                GDALCopyWords(adfPixVal, GDT_CFloat64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        /* non complex */
        double dfPixVal;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {
                dfPixVal = 1;

                for( iSrc = 0; iSrc < nSources; ++iSrc ) {
                    /* Source raster pixels may be obtained with SRCVAL macro */
                    dfPixVal *= SRCVAL(papoSources[iSrc], eSrcType, ii);
                }

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* MulPixelFunc */


CPLErr CMulPixelFunc(void **papoSources, int nSources, void *pData,
                     int nXSize, int nYSize,
                     GDALDataType eSrcType, GDALDataType eBufType,
                     int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;

    /* ---- Init ---- */
    if (nSources != 2) return CE_Failure;

    /* ---- Set pixels ---- */
    if (GDALDataTypeIsComplex( eSrcType ))
    {
        double adfPixVal[2], dfReal0, dfImag0, dfReal1, dfImag1;

        int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;
        void *pReal0 = papoSources[0];
        void *pImag0 = ((GByte *)papoSources[0]) + nOffset;
        void *pReal1 = papoSources[1];
        void *pImag1 = ((GByte *)papoSources[1]) + nOffset;

        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfReal0 = SRCVAL(pReal0, eSrcType, ii);
                dfReal1 = SRCVAL(pReal1, eSrcType, ii);
                dfImag0 = SRCVAL(pImag0, eSrcType, ii);
                dfImag1 = SRCVAL(pImag1, eSrcType, ii);
                adfPixVal[0]  = dfReal0 * dfReal1 + dfImag0 * dfImag1;
                adfPixVal[1]  = dfReal1 * dfImag0 - dfReal0 * dfImag1;

                GDALCopyWords(adfPixVal, GDT_CFloat64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        /* non complex */
        double adfPixVal[2] = {0, 0};

        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                adfPixVal[0] = SRCVAL(papoSources[0], eSrcType, ii)
                             * SRCVAL(papoSources[1], eSrcType, ii);

                GDALCopyWords(adfPixVal, GDT_CFloat64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* CMulPixelFunc */


CPLErr InvPixelFunc(void **papoSources, int nSources, void *pData,
                    int nXSize, int nYSize,
                    GDALDataType eSrcType, GDALDataType eBufType,
                    int nPixelSpace, int nLineSpace)
{
    int iLine, iCol, ii;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;

    /* ---- Set pixels ---- */
    if (GDALDataTypeIsComplex( eSrcType ))
    {
        double adfPixVal[2], dfReal, dfImag, dfAux;

        int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;
        void *pReal = papoSources[0];
        void *pImag = ((GByte *)papoSources[0]) + nOffset;

        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfReal = SRCVAL(pReal, eSrcType, ii);
                dfImag = SRCVAL(pImag, eSrcType, ii);
                dfAux = dfReal * dfReal + dfImag * dfImag;
                adfPixVal[0]  = +dfReal / dfAux;
                adfPixVal[1]  = -dfImag / dfAux;

                GDALCopyWords(adfPixVal, GDT_CFloat64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        /* non complex */
        double dfPixVal;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfPixVal = 1. / SRCVAL(papoSources[0], eSrcType, ii);

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* InvPixelFunc */


CPLErr IntensityPixelFunc(void **papoSources, int nSources, void *pData,
                          int nXSize, int nYSize,
                          GDALDataType eSrcType, GDALDataType eBufType,
                          int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;
    double dfPixVal;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;

    if (GDALDataTypeIsComplex( eSrcType ))
    {
        double dfReal, dfImag;
        int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;
        void *pReal = papoSources[0];
        void *pImag = ((GByte *)papoSources[0]) + nOffset;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii = 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfReal = SRCVAL(pReal, eSrcType, ii);
                dfImag = SRCVAL(pImag, eSrcType, ii);

                dfPixVal = dfReal * dfReal + dfImag * dfImag;

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        /* ---- Set pixels ---- */
        for( iLine = 0, ii = 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfPixVal = SRCVAL(papoSources[0], eSrcType, ii);
                dfPixVal *= dfPixVal;

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* IntensityPixelFunc */


CPLErr SqrtPixelFunc(void **papoSources, int nSources, void *pData,
                     int nXSize, int nYSize,
                     GDALDataType eSrcType, GDALDataType eBufType,
                     int nPixelSpace, int nLineSpace)
{
    int iLine, iCol, ii;
    double dfPixVal;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;
    if (GDALDataTypeIsComplex( eSrcType )) return CE_Failure;

    /* ---- Set pixels ---- */
    for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
        for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

            /* Source raster pixels may be obtained with SRCVAL macro */;
            dfPixVal = sqrt( SRCVAL(papoSources[0], eSrcType, ii) );

            GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine +
                          iCol * nPixelSpace, eBufType, nPixelSpace, 1);
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* SqrtPixelFunc */


CPLErr Log10PixelFunc(void **papoSources, int nSources, void *pData,
                      int nXSize, int nYSize,
                      GDALDataType eSrcType, GDALDataType eBufType,
                      int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;

    if (GDALDataTypeIsComplex( eSrcType ))
    {
        /* complex input datatype */
        double dfReal, dfImag, dfPixVal;
        int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;
        void *pReal = papoSources[0];
        void *pImag = ((GByte *)papoSources[0]) + nOffset;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii = 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfReal = SRCVAL(pReal, eSrcType, ii);
                dfImag = SRCVAL(pImag, eSrcType, ii);

                dfPixVal = log10( dfReal * dfReal + dfImag * dfImag );

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        double dfPixVal;

        /* ---- Set pixels ---- */
        for( iLine = 0, ii = 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                /* Source raster pixels may be obtained with SRCVAL macro */
                dfPixVal = SRCVAL(papoSources[0], eSrcType, ii);
                dfPixVal = log10( abs( dfPixVal ) );

                GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* Log10PixelFunc */


CPLErr PowPixelFuncHelper(void **papoSources, int nSources, void *pData,
                          int nXSize, int nYSize,
                          GDALDataType eSrcType, GDALDataType eBufType,
                          int nPixelSpace, int nLineSpace,
                          double base, double fact)
{
    int iLine, iCol, ii;
    double dfPixVal;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;
    if (GDALDataTypeIsComplex( eSrcType )) return CE_Failure;

    /* ---- Set pixels ---- */
    for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
        for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

            /* Source raster pixels may be obtained with SRCVAL macro */
            dfPixVal = SRCVAL(papoSources[0], eSrcType, ii);
            dfPixVal = pow(base, dfPixVal / fact);

            GDALCopyWords(&dfPixVal, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine +
                          iCol * nPixelSpace, eBufType, nPixelSpace, 1);
        }
    }

    /* ---- Return success ---- */
    return CE_None;
} /* PowPixelFuncHelper */

CPLErr dB2AmpPixelFunc(void **papoSources, int nSources, void *pData,
                       int nXSize, int nYSize,
                       GDALDataType eSrcType, GDALDataType eBufType,
                       int nPixelSpace, int nLineSpace)
{
    return PowPixelFuncHelper(papoSources, nSources, pData,
                              nXSize, nYSize, eSrcType, eBufType,
                              nPixelSpace, nLineSpace, 10., 20.);
} /* dB2AmpPixelFunc */


CPLErr dB2PowPixelFunc(void **papoSources, int nSources, void *pData,
                       int nXSize, int nYSize,
                       GDALDataType eSrcType, GDALDataType eBufType,
                       int nPixelSpace, int nLineSpace)
{
    return PowPixelFuncHelper(papoSources, nSources, pData,
                              nXSize, nYSize, eSrcType, eBufType,
                              nPixelSpace, nLineSpace, 10., 10.);
} /* dB2PowPixelFunc */

/************************************************************************/
/*                     Nansat pixelfunctions                            */
/************************************************************************/

CPLErr BetaSigmaToIncidence(void **papoSources, int nSources, void *pData,
        int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;
    double incidence;
    double beta0, sigma0;

    /* ---- Init ---- */
    if (nSources != 2) return CE_Failure;
    #define PI 3.14159265;

        /*printf("%d",eSrcType);*/

        if (GDALDataTypeIsComplex( eSrcType ))
        {
            double b0Real, b0Imag;
            double s0Real, s0Imag;
            void *b0pReal = papoSources[0];
            void *s0pReal = papoSources[1];
            void *b0pImag = ((GByte *)papoSources[0])
                        + GDALGetDataTypeSize( eSrcType ) / 8 / 2;
            void *s0pImag = ((GByte *)papoSources[1])
                        + GDALGetDataTypeSize( eSrcType ) / 8 / 2;

            /* ---- Set pixels ---- */
            for( iLine = 0, ii = 0; iLine < nYSize; ++iLine ) {
                for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                    /* Source raster pixels may be obtained with SRCVAL macro */
                    b0Real = SRCVAL(b0pReal, eSrcType, ii);
                    b0Imag = SRCVAL(b0pImag, eSrcType, ii);
                    s0Real = SRCVAL(s0pReal, eSrcType, ii);
                    s0Imag = SRCVAL(s0pImag, eSrcType, ii);

                    beta0 = b0Real*b0Real + b0Imag*b0Imag;
                    sigma0 = s0Real*s0Real + s0Imag*s0Imag;

            if (beta0 != 0) incidence = asin(sigma0/beta0)*180/PI
            else incidence = -10000; // NB: this is also hard-coded in
                                             //     mapper_radarsat2.py, and
                                             //     should be the same in other
                                             //     mappers where this function
                                             //     is needed...
            GDALCopyWords(&incidence, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                          eBufType, nPixelSpace, 1);
                }
            }
        } else {

        /* ---- Set pixels ---- */
        for( iLine = 0; iLine < nYSize; iLine++ )
        {
        for( iCol = 0; iCol < nXSize; iCol++ )
        {
            ii = iLine * nXSize + iCol;
            /* Source raster pixels may be obtained with SRCVAL macro */
            beta0 = SRCVAL(papoSources[0], eSrcType, ii);
            sigma0 = SRCVAL(papoSources[1], eSrcType, ii);

            if (beta0 != 0) incidence = asin(sigma0/beta0)*180/PI
            else incidence = -10000; // NB: this is also hard-coded in
                                                 //     mapper_radarsat2.py, and
                                                 //     should be the same in other
                                                 //     mappers where this function
                                     //     is needed...
                        GDALCopyWords(&incidence, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                          eBufType, nPixelSpace, 1);
        }
        }
        }

    /* ---- Return success ---- */
        return CE_None;
}


CPLErr UVToMagnitude(void **papoSources, int nSources, void *pData,
        int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;
    double magnitude;
    double u, v;

    /* ---- Init ---- */
    if (nSources != 2) return CE_Failure;

    /* ---- Set pixels ---- */
    for( iLine = 0; iLine < nYSize; iLine++ )
    {
        for( iCol = 0; iCol < nXSize; iCol++ )
        {
            ii = iLine * nXSize + iCol;
            /* Source raster pixels may be obtained with SRCVAL macro */
            u = SRCVAL(papoSources[0], eSrcType, ii);
            v = SRCVAL(papoSources[1], eSrcType, ii);

            magnitude = sqrt(u*u + v*v);

            GDALCopyWords(&magnitude, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                          eBufType, nPixelSpace, 1);
        }
    }

    /* ---- Return success ---- */
return CE_None;
}



CPLErr Sigma0HHBetaToSigma0VV(void **papoSources, int nSources, void *pData,
        int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;
    double sigma0HH, beta0, incidence, factor, sigma0VV;

    /* ---- Init ---- */
    if (nSources != 2) return CE_Failure;
        /*fprintf("nSources: %d\n", nSources);*/

    /* ---- Set pixels ---- */
    for( iLine = 0; iLine < nYSize; iLine++ )
    {
        for( iCol = 0; iCol < nXSize; iCol++ )
        {
            ii = iLine * nXSize + iCol;
            /* Source raster pixels may be obtained with SRCVAL macro */
            sigma0HH = SRCVAL(papoSources[0], eSrcType, ii);
            beta0 = SRCVAL(papoSources[1], eSrcType, ii);

                        /* get incidence angle first */
                    if (beta0 != 0){
                            incidence = asin(sigma0HH/beta0);
                        } else {
                            incidence = 0;
                        }

                    /* Polarisation ratio from Thompson et al. with alpha=1 */
                    factor = pow( (1 + 2 * pow(tan(incidence), 2)) / (1 + 1 * pow(tan(incidence), 2)), 2);
                    sigma0VV = sigma0HH * factor;

                    GDALCopyWords(&sigma0VV, GDT_Float64, 0,
                                ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                            eBufType, nPixelSpace, 1);
        }
    }

    /* ---- Return success ---- */
        return CE_None;
}


CPLErr RawcountsToSigma0_CosmoSkymed_SBI(void **papoSources, int nSources, void *pData,
        int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{

    int ii, iLine, iCol;
    /* int iReal, iImag; */
    double imPower, real, imag;

    /* ---- Init ---- */
    if (nSources != 2) return CE_Failure;

    /* ---- Set pixels ---- */
    for( iLine = 0; iLine < nYSize; iLine++ ){
        for( iCol = 0; iCol < nXSize; iCol++ ){
        ii = iLine * nXSize + iCol;
        /* Source raster pixels may be obtained with SRCVAL macro */
        real = SRCVAL(papoSources[0], eSrcType, ii);
        imag = SRCVAL(papoSources[1], eSrcType, ii);

            /*printf("%d",iReal); OK!*/

            /*real = (double) iReal;*/
            /*imag = (double) iImag;*/

            /*printf("%.1f",imag); OK!*/

            imPower = pow(real,2.0) + pow(imag,2.0);
            /*printf("%.1f",imPower); */

        GDALCopyWords(&imPower, GDT_Float64, 0,
                ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                eBufType, nPixelSpace, 1);

        }
    }

    /* ---- Return success ---- */
    return CE_None;

}

CPLErr RawcountsToSigma0_CosmoSkymed_QLK(void **papoSources, int nSources, void *pData,
        int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{

    int ii, iLine, iCol;
    double raw_counts, imPower;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;

    /* ---- Set pixels ---- */
    for( iLine = 0; iLine < nYSize; iLine++ )
    {
        for( iCol = 0; iCol < nXSize; iCol++ )
        {
            ii = iLine * nXSize + iCol;
            /* Source raster pixels may be obtained with SRCVAL macro */
            raw_counts = SRCVAL(papoSources[0], eSrcType, ii);
                        imPower = pow(raw_counts,2.);

            GDALCopyWords(&imPower, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                          eBufType, nPixelSpace, 1);
                }
        }

    /* ---- Return success ---- */
    return CE_None;

}


CPLErr ComplexData(void **papoSources, int nSources, void *pData,
        int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;
    double adfPixVal[2];
    void *pReal = papoSources[0];
    void *pImag = papoSources[1];

    for( iLine = 0, ii= 0; iLine < nYSize; ++iLine ) {
        for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

            /* Source raster pixels may be obtained with SRCVAL macro */
            adfPixVal[0] = SRCVAL(pReal, eSrcType, ii);
            adfPixVal[1] = SRCVAL(pImag, eSrcType, ii);

            GDALCopyWords(&adfPixVal, GDT_CFloat64, 0,
                          ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                          eBufType, nPixelSpace, 1);
        }
    }

    /* ---- Return success ---- */
return CE_None;
}

CPLErr IntensityInt(void **papoSources, int nSources, void *pData,
                    int nXSize, int nYSize,
                    GDALDataType eSrcType, GDALDataType eBufType,
                    int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol;
    int dfPixVal;
    int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;
    void *pReal = papoSources[0];
    void *pImag = ((GByte *)papoSources[0]) + nOffset;

    /* ---- Init ---- */
    if (nSources != 1) return CE_Failure;


        // ---- Set pixels ----

    if (GDALDataTypeIsComplex( eSrcType ))
    {
        int dfReal, dfImag;
        int nOffset = GDALGetDataTypeSize( eSrcType ) / 8 / 2;
        void *pReal = papoSources[0];
        void *pImag = ((GByte *)papoSources[0]) + nOffset;

        // ---- Set pixels ----

        for( iLine = 0, ii = 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                // Source raster pixels may be obtained with SRCVAL macro
                dfReal = SRCVAL(pReal, eSrcType, ii);
                dfImag = SRCVAL(pImag, eSrcType, ii);

                dfPixVal = dfReal * dfReal + dfImag * dfImag;

                GDALCopyWords(&dfPixVal, GDT_Int16, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    } else {
        // ---- Set pixels ----
        for( iLine = 0, ii = 0; iLine < nYSize; ++iLine ) {
            for( iCol = 0; iCol < nXSize; ++iCol, ++ii ) {

                // Source raster pixels may be obtained with SRCVAL macro
                dfPixVal = SRCVAL(papoSources[0], eSrcType, ii);
                dfPixVal *= dfPixVal;

                GDALCopyWords(&dfPixVal, GDT_Int32, 0,
                              ((GByte *)pData) + nLineSpace * iLine +
                              iCol * nPixelSpace, eBufType, nPixelSpace, 1);
            }
        }
    }
    /* ---- Return success ---- */
    return CE_None;
} /* IntensityInt */


CPLErr OnesPixelFunc(void **papoSources, int nSources, void *pData,
                    int nXSize, int nYSize,
                    GDALDataType eSrcType, GDALDataType eBufType,
                    int nPixelSpace, int nLineSpace)
{
    char one=1;
    int iLine, iCol;

    /* ---- Set all pixels to 1 ---- */
    for( iLine = 0; iLine < nYSize; iLine++ ){
        for( iCol = 0; iCol < nXSize; iCol++ ){

        GDALCopyWords(&one, GDT_Byte, 0,
                ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                eBufType, nPixelSpace, 1);
        }
    }
    /* ---- Return success ---- */
    return CE_None;
}



/************************************************************************/
/*                       Convert Rrs to Rrsw                            */
/************************************************************************/
/* scientifc function */
double NormReflectanceToRemSensReflectanceFunction(double *b){
    return b[0] / (0.52 + 1.7 * b[0]);
}

double RawcountsIncidenceToSigma0Function(double *b){
    double pi = 3.14159265;
    return (pow(b[0], 2.0) * sin(b[1] *  pi / 180.0));
}

double Sentinel1CalibrationFunction(double *b){

    // With noise removal -- I am not sure if the noise (b[2]) should be
    // squared or not but have sent an email to esa..
    //return ( pow(b[1],2.0) - pow(b[2],2.0) ) / pow(b[0], 2.0);
    // Without noise removal
    return pow(b[0],2.0) / pow(b[1], 2.0);

}

double Sigma0HHToSigma0VVFunction(double *b){
    double pi = 3.14159265;
    double s0hh, factor;
    s0hh = (pow(b[0], 2.0) * sin(b[1] *  pi / 180.0));
    /* Polarisation ratio from Thompson et al. with alpha=1 */
    factor = pow( (1 + 2 * pow(tan(b[1]*pi/180.0), 2)) / (1 + 1 * pow(tan(b[1]*pi/180.0), 2)), 2);
    return s0hh * factor;
}

double Sentinel1Sigma0HHToSigma0VVFunction( double *b ){

    double s0hh, s0vv;
    double bcal[3];
    double s0hh2s0vv[2];

    bcal[0] = b[0]; // sigmaNought LUT
    bcal[1] = b[2]; // DN
    //bcal[2] = b[3]; // noise
    s0hh = Sentinel1CalibrationFunction(bcal);

    s0hh2s0vv[0] = s0hh;
    s0hh2s0vv[1] = b[1];

    s0vv = Sigma0HHToSigma0VVFunction(s0hh2s0vv);

    return s0vv;

}

double Sigma0NormalizedIceFunction(double *b){
    double pi = 3.14159265;
    double sigma0 = (pow(b[0], 2.0) * sin(b[1] *  pi / 180.0));
    return sigma0 * pow((tan(b[1] * pi / 180.0) / tan(31.0 * pi / 180.0)), 1.5);
}

double Sigma0VVNormalizedWaterFunction(double *b){
    double pi = 3.14159265;
    double sigma0 = (pow(b[0], 2.0) * sin(b[1] *  pi / 180.0));
    return sigma0 * pow((sin(b[1] * pi / 180.0) / sin(31.0 * pi / 180.0)), 4.0);
}

double Sigma0HHNormalizedWaterFunction(double *b){
    double pi = 3.14159265;
    double sigma0 = (pow(b[0], 2.0) * sin(b[1] *  pi / 180.0));
    return sigma0 * pow((tan(b[1] * pi / 180.0) / tan(31.0 * pi / 180.0)), 4.0);
}

double UVToDirectionFromFunction(double *b){
        /* Convention 0-360 degrees positive clockwise from north*/
    double pi = 3.14159265;
    //return (b[0]==9999 || b[1]==9999) ? 9999 : 180.0 - atan2(-b[0],b[1])*180./pi;
    return 180.0 - atan2(-b[0],b[1])*180./pi;
}

double UVToDirectionToFunction(double *b){
        /* Convention 0-360 degrees positive clockwise from north*/
    double pi = 3.14159265;
    return 360.0 - atan2(-b[0],b[1])*180./pi;
        /*
           Below code is hirlam specific - we don't know if the invalid data is
           actually 9999. One option is to make mapper specific pixelfunctions
           but for now only return the direction as if all data was good.
        */
    //return (b[0]==9999 || b[1]==9999) ? 9999 : 360.0 - atan2(-b[0],b[1])*180./pi;
}

/* pixel function */
CPLErr UVToDirectionTo(void **papoSources, int nSources, void *pData,
        int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace){

    GenericPixelFunction(UVToDirectionToFunction,
        papoSources, nSources,  pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}

CPLErr UVToDirectionFrom(void **papoSources, int nSources, void *pData,
        int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace){

    GenericPixelFunction(UVToDirectionFromFunction,
        papoSources, nSources,  pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}


CPLErr NormReflectanceToRemSensReflectance(void **papoSources, int nSources, void *pData,
        int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace){

    GenericPixelFunction(NormReflectanceToRemSensReflectanceFunction,
        papoSources, nSources,  pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}

CPLErr Sentinel1Calibration(void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
                GDALDataType eSrcType, GDALDataType eBufType,
                int nPixelSpace, int nLineSpace){

    GenericPixelFunction(Sentinel1CalibrationFunction,
        papoSources, nSources, pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}

CPLErr Sentinel1Sigma0HHToSigma0VV(void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
                GDALDataType eSrcType, GDALDataType eBufType,
                int nPixelSpace, int nLineSpace){

    GenericPixelFunction(Sentinel1Sigma0HHToSigma0VVFunction,
        papoSources, nSources, pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}

CPLErr RawcountsIncidenceToSigma0(void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace){

    GenericPixelFunction(RawcountsIncidenceToSigma0Function,
        papoSources, nSources, pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}

CPLErr Sigma0HHToSigma0VV(void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
                GDALDataType eSrcType, GDALDataType eBufType,
                int nPixelSpace, int nLineSpace){
    // Works for ASAR!
    GenericPixelFunction(Sigma0HHToSigma0VVFunction,
        papoSources, nSources, pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}

CPLErr Sigma0NormalizedIce(void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace){

    GenericPixelFunction(Sigma0NormalizedIceFunction,
        papoSources, nSources, pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}

CPLErr Sigma0VVNormalizedWater(void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace){

    GenericPixelFunction(Sigma0VVNormalizedWaterFunction,
        papoSources, nSources, pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}

CPLErr Sigma0HHNormalizedWater(void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace){

    GenericPixelFunction(Sigma0HHNormalizedWaterFunction,
        papoSources, nSources, pData,
        nXSize, nYSize, eSrcType, eBufType,
        nPixelSpace, nLineSpace);

    return CE_None;
}



/************************************************************************/
/* Generic Pixel Function is called from a pixel function and calls
 * corresponding scientific function */
/************************************************************************/

// all data (band) size must be same and full size of bands (XSize x YSize).
void GenericPixelFunction(double f(double*), void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{
    int ii, iLine, iCol, iSrc;
    double *bVal, result;
    bVal = malloc(nSources * sizeof (double));

    /* ---- Set pixels ---- */
    for( iLine = 0; iLine < nYSize; iLine++ ){
        for( iCol = 0; iCol < nXSize; iCol++ ){
        ii = iLine * nXSize + iCol;
        /* Source raster pixels may be obtained with SRCVAL macro */
            for (iSrc = 0; iSrc < nSources; iSrc ++){
                bVal[iSrc] = SRCVAL(papoSources[iSrc], eSrcType, ii);
                //if (iLine==0 && iCol==0){
                //    printf("%d ",iSrc);
                //    printf("%.4f\n",bVal[iSrc]);
                //}
            }

        result = f(bVal);

        GDALCopyWords(&result, GDT_Float64, 0,
                ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                eBufType, nPixelSpace, 1);
    }
    }
}

// From the 1st to (N-1)th bands are full size (XSize x YSize),
// and the last band is a one-pixel band (1 x 1).
void GenericPixelFunctionPixel(double f(double*), void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{
    int iLine, iCol, iSrc;
    double *bVal, result;
    bVal = malloc(nSources * sizeof (double));

    /* ---- Set pixels ---- */
    /* Set the first value form one-pixel band */
    bVal[0] = SRCVAL(papoSources[nSources-1], eSrcType, 0);
    for( iLine = 0; iLine < nYSize; iLine++ )
    {
        for( iCol = 0; iCol < nXSize; iCol++ )
        {
            for (iSrc = 1; iSrc < nSources; iSrc ++)
                /* Source raster pixels may be obtained with SRCVAL macro */
                bVal[iSrc] = SRCVAL(papoSources[iSrc-1], eSrcType, iLine * nXSize + iCol);

            result = f(bVal);

            GDALCopyWords(&result, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                          eBufType, nPixelSpace, 1);
        }
    }
}

// From the 1st to (N-1)th bands are full size (XSize x YSize),
// and the last band is a line band (XSize x 1).
void GenericPixelFunctionLine(double f(double*), void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{
    int iLine, iCol, iSrc;
    double *bVal, result;
    bVal = malloc(nSources * sizeof (double));

    /* ---- Set pixels ---- */
    for( iLine = 0; iLine < nYSize; iLine++ )
    {
        for( iCol = 0; iCol < nXSize; iCol++ )
        {
            /* Source raster pixels may be obtained with SRCVAL macro */
            bVal[0] = SRCVAL(papoSources[nSources-1], eSrcType, iCol);

            for (iSrc = 1; iSrc < nSources; iSrc ++)
                /* Source raster pixels may be obtained with SRCVAL macro */
                bVal[iSrc] = SRCVAL(papoSources[iSrc-1], eSrcType, iLine * nXSize + iCol);

            result = f(bVal);

            GDALCopyWords(&result, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                          eBufType, nPixelSpace, 1);
        }
    }
}

// From the 1st to (N-2)th bands are full size (XSize x YSize),
// the last 2nd band is a line band (XSize x 1) and the last is one pixel band.
void GenericPixelFunctionPixelLine(double f(double*), void **papoSources,
        int nSources, void *pData, int nXSize, int nYSize,
        GDALDataType eSrcType, GDALDataType eBufType,
        int nPixelSpace, int nLineSpace)
{
    int iLine, iCol, iSrc;
    double  *bVal, result;
    bVal = malloc(nSources * sizeof (double));

    /* ---- Set pixels ---- */
    bVal[0] = SRCVAL(papoSources[nSources-1], eSrcType, 0);
    for( iLine = 0; iLine < nYSize; iLine++ )
    {
        for( iCol = 0; iCol < nXSize; iCol++ )
        {
            bVal[1] = SRCVAL(papoSources[nSources-2], eSrcType, iCol);

            for(iSrc = 2; iSrc < nSources; iSrc++ )
                bVal[iSrc] = SRCVAL(papoSources[iSrc-2], eSrcType, iLine * nXSize + iCol);

            result = f(bVal);

            GDALCopyWords(&result, GDT_Float64, 0,
                          ((GByte *)pData) + nLineSpace * iLine + iCol * nPixelSpace,
                          eBufType, nPixelSpace, 1);
        }
    }
}



/************************************************************************/
/*                     GDALRegisterDefaultPixelFunc()                   */
/************************************************************************/

/**
 * This adds a default set of pixel functions to the global list of
 * available pixel functions for derived bands:
 *
 * - "real": extract real part from a single raster band (just a copy if the
 *           input is non-complex)
 * - "imag": extract imaginary part from a single raster band (0 for
 *           non-complex)
 * - "mod": extract module from a single raster band (real or complex)
 * - "phase": extract phase from a single raster band (0 for non-complex)
 * - "conj": computes the complex conjugate of a single raster band (just a
 *           copy if the input is non-complex)
 * - "sum": sum 2 or more raster bands
 * - "diff": computes the difference between 2 raster bands (b1 - b2)
 * - "mul": multilpy 2 or more raster bands
 * - "cmul": multiply the first band for the complex comjugate of the second
 * - "inv": inverse (1./x). Note: no check is performed on zero division
 * - "intensity": computes the intensity Re(x*conj(x)) of a single raster band
 *                (real or complex)
 * - "sqrt": perform the square root of a single raster band (real only)
 * - "log10": compute the logarithm (base 10) of the abs of a single raster
 *            band (real or complex): log10( abs( x ) )
 * - "dB2amp": perform scale conversion from logarithmic to linear
 *             (amplitude) (i.e. 10 ^ ( x / 20 ) ) of a single raster
 *                 band (real only)
 * - "dB2pow": perform scale conversion from logarithmic to linear
 *             (power) (i.e. 10 ^ ( x / 10 ) ) of a single raster
 *             band (real only)
 *
 * @see GDALAddDerivedBandPixelFunc
 *
 * @return CE_None, invalid (NULL) parameters are currently ignored.
 */
CPLErr CPL_STDCALL GDALRegisterDefaultPixelFunc()
{
    GDALAddDerivedBandPixelFunc("real", RealPixelFunc);
    GDALAddDerivedBandPixelFunc("imag", ImagPixelFunc);
    GDALAddDerivedBandPixelFunc("mod", ModulePixelFunc);
    GDALAddDerivedBandPixelFunc("phase", PhasePixelFunc);
    GDALAddDerivedBandPixelFunc("conj", ConjPixelFunc);
    GDALAddDerivedBandPixelFunc("sum", SumPixelFunc);
    GDALAddDerivedBandPixelFunc("diff", DiffPixelFunc);
    GDALAddDerivedBandPixelFunc("mul", MulPixelFunc);
    GDALAddDerivedBandPixelFunc("cmul", CMulPixelFunc);
    GDALAddDerivedBandPixelFunc("inv", InvPixelFunc);
    GDALAddDerivedBandPixelFunc("intensity", IntensityPixelFunc);
    GDALAddDerivedBandPixelFunc("sqrt", SqrtPixelFunc);
    GDALAddDerivedBandPixelFunc("log10", Log10PixelFunc);
    GDALAddDerivedBandPixelFunc("dB2amp", dB2AmpPixelFunc);
    GDALAddDerivedBandPixelFunc("dB2pow", dB2PowPixelFunc);

    GDALAddDerivedBandPixelFunc("BetaSigmaToIncidence", BetaSigmaToIncidence);
    GDALAddDerivedBandPixelFunc("UVToMagnitude", UVToMagnitude);
    GDALAddDerivedBandPixelFunc("UVToDirectionTo", UVToDirectionTo);
    GDALAddDerivedBandPixelFunc("UVToDirectionFrom", UVToDirectionFrom);
    GDALAddDerivedBandPixelFunc("Sigma0HHBetaToSigma0VV", Sigma0HHBetaToSigma0VV); //Radarsat-2
    GDALAddDerivedBandPixelFunc("Sigma0HHToSigma0VV", Sigma0HHToSigma0VV); // ASAR
    GDALAddDerivedBandPixelFunc("RawcountsIncidenceToSigma0", RawcountsIncidenceToSigma0);
    GDALAddDerivedBandPixelFunc("RawcountsToSigma0_CosmoSkymed_QLK", RawcountsToSigma0_CosmoSkymed_QLK);
    GDALAddDerivedBandPixelFunc("RawcountsToSigma0_CosmoSkymed_SBI", RawcountsToSigma0_CosmoSkymed_SBI);
    GDALAddDerivedBandPixelFunc("ComplexData", ComplexData);
    GDALAddDerivedBandPixelFunc("NormReflectanceToRemSensReflectance", NormReflectanceToRemSensReflectance);
    GDALAddDerivedBandPixelFunc("Sigma0NormalizedIce", Sigma0NormalizedIce);
    GDALAddDerivedBandPixelFunc("Sigma0HHNormalizedWater", Sigma0HHNormalizedWater);
    GDALAddDerivedBandPixelFunc("Sigma0VVNormalizedWater", Sigma0VVNormalizedWater);
    GDALAddDerivedBandPixelFunc("Sentinel1Calibration", Sentinel1Calibration);
    GDALAddDerivedBandPixelFunc("Sentinel1Sigma0HHToSigma0VV", Sentinel1Sigma0HHToSigma0VV);
    GDALAddDerivedBandPixelFunc("IntensityInt", IntensityInt);
    GDALAddDerivedBandPixelFunc("OnesPixelFunc", OnesPixelFunc);
    return CE_None;
}

