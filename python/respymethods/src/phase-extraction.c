#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>


static PyObject *times2array(PyObject *self, PyObject *args) {

	PyArrayObject* array;
    
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &array)) {
        return NULL;
    }

    // Get input array properties
    int ndim = PyArray_NDIM(array);
    npy_intp* dims = PyArray_DIMS(array);
    int dtype = PyArray_TYPE(array);
	npy_intp size = PyArray_SIZE(array);

    double *datain = (double*)PyArray_DATA(array);


    // Create new array: PyArray_SimpleNew(ndim, dims, typenum)
    PyObject *result = PyArray_SimpleNew(ndim, dims, dtype);
    double *dataout = (double*)PyArray_DATA((PyArrayObject*)result);


	for (npy_intp i=0; i<size; i++){
		
		dataout[i] = datain[i]*2.0;		
	}

	return result;
	
}



// Method definitions
static PyMethodDef cpytest_methods[] = {
	{"times2array", times2array, METH_VARARGS, "multiply times 2 each element of an array"},
	{NULL, NULL, 0, NULL}
};

// Module defintion
static struct PyModuleDef PhaseExtraction = {
	PyModuleDef_HEAD_INIT,
	"PhaseExtraction", // __name__
	"Hello", // __doc__
	32, // m_size, set to -1 meaning that the module does not support sub-interpreters. https://docs.python.org/3/c-api/module.html#c.PyModuleDef
	cpytest_methods
};

// Module inititalization
PyMODINIT_FUNC PyInit_PhaseExtraction(void) {
	import_array();
	return PyModuleDef_Init(&PhaseExtraction);
}


