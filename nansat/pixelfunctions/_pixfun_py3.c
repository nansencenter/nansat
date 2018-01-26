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
/* deprecated in Py3
static PyMethodDef module_methods[] = {
	{"registerPixelFunctions", registerPixelFunctions, METH_VARARGS, pixfun_docstring},
	{NULL, NULL, 0, NULL}
};
*/

static PyMethodDef module_methods[] = {
    {"registerPixelFunctions", (PyCFunction) registerPixelFunctions, METH_NOARGS, pixfun_docstring},
    {NULL, NULL, 0, NULL}
};


/* Initialize the module */
/* DEPRECATED
PyMODINIT_FUNC init_pixfun(void)
{
	PyObject *m = Py_InitModule3("_pixfun", module_methods, module_docstring);
	if (m == NULL)
		return;
}
*/

/* deprecated :
PyMODINIT_FUNC init_uniqueCombinations(void)
{
    Py_InitModule3("uniqueCombinations", uniqueCombinations_funcs,
                   "Extension module uniqueCombinations v. 0.01");
}
*/

static struct PyModuleDef _pixfun_py3 =
{
    PyModuleDef_HEAD_INIT,
    "_pixfun_py3", /* name of module */
    "usage: _pixfun_py3.registerPixelFunctions\n", /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit__pixfun_py3(void)
{
    return PyModule_Create(&_pixfun_py3);
}



static PyObject *registerPixelFunctions(PyObject *self, PyObject *args)
{
	GDALRegisterDefaultPixelFunc();
	Py_INCREF(Py_None);
	return Py_None;
}

/***********************************/

/* deprecated:
static PyMethodDef uniqueCombinations_funcs[] = {
    {"uniqueCombinations", (PyCFunction)uniqueCombinations,
     METH_NOARGS, uniqueCombinations_docs},
    {NULL}
};
use instead of the above: */

/* NEW
static PyMethodDef module_methods[] = {
    {"uniqueCombinations", (PyCFunction) uniqueCombinations, METH_NOARGS, uniqueCombinations_docs},
    {NULL}
};
*/


/* deprecated :
PyMODINIT_FUNC init_uniqueCombinations(void)
{
    Py_InitModule3("uniqueCombinations", uniqueCombinations_funcs,
                   "Extension module uniqueCombinations v. 0.01");
}
*/

/* NEW
static struct PyModuleDef Combinations =
{
    PyModuleDef_HEAD_INIT,
    "Combinations", // name of module
    "usage: Combinations.uniqueCombinations(lstSortableItems, comboSize)\n", // module documentation, may be NULL
    -1,   // size of per-interpreter state of the module, or -1 if the module keeps state in global variables.
    module_methods
};

PyMODINIT_FUNC PyInit_Combinations(void)
{
    return PyModule_Create(&Combinations);
}
*/

