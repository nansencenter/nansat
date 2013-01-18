/******************************************************************************
 * $Id: pixfunplugin.c 2958 2011-02-11 10:45:29Z valentino $
 *
 * Project:  GDAL
 * Purpose:  Provide a fake GDAL driver to register a small set of pixel
 *           functions to be used with the virtual driver.
 *           It is a hack aimed to enable python users to use a small set
 *           of custom pixel functions without C++ coding.
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

#include <gdal.h>

/* CPL_CVSID("$Id: pixfunplugin.c 2958 2011-02-11 10:45:29Z valentino $"); */

CPL_C_START
CPL_DLL void GDALRegister_PIXFUN(void);
CPL_C_END

extern CPLErr CPL_STDCALL GDALRegisterDefaultPixelFunc();

/************************************************************************/
/*                         GDALRegister_PIXFUN()                        */
/************************************************************************/
void GDALRegister_PIXFUN(void)
{
    const char PLUGINNAME[]="Pixfun";

    if (! GDAL_CHECK_VERSION(PLUGINNAME))
        return;

    GDALRegisterDefaultPixelFunc();
    CPLDebug("PIXFUN", "Plugin %s %s", PLUGINNAME, "$Revision: 2958 $");

} /* GDALRegister_PIXFUN */
