#include <Python.h>
#include <gdal.h>

extern CPLErr CPL_STDCALL GDALRegisterDefaultPixelFunc();

/* Docstrings */
static char module_docstring[] =
	"";
static char pixfun_docstring[] =
	"";

/* The only available function */
static PyObject *registerPixelFunctions(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
	{"registerPixelFunctions", registerPixelFunctions, METH_VARARGS, pixfun_docstring},
	{NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC init_pixfun_py2(void)
{
	PyObject *m = Py_InitModule3("_pixfun_py2", module_methods, module_docstring);
	if (m == NULL)
		return;
}

static PyObject *registerPixelFunctions(PyObject *self, PyObject *args)
{
	GDALRegisterDefaultPixelFunc();
	Py_INCREF(Py_None);
	return Py_None;
}
