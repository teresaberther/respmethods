#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>


static PyObject *two_point_interp(PyObject *self, PyObject *args) {

	PyArrayObject* array;
    
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &array)) {
        return NULL;
    }

    // Get input array properties
    int ndim = PyArray_NDIM(array);
    npy_intp* dims = PyArray_DIMS(array);
    int dtype = PyArray_TYPE(array);
	npy_intp arr_len = PyArray_SIZE(array);
	int typenum_int = NPY_INT;

    double *datain = (double*)PyArray_DATA(array);

	// allocate an array for events for loc min and max (null derivative)
    PyObject *events_arr = PyArray_SimpleNew(ndim, dims, typenum_int);
	int *null_der = (int*)PyArray_DATA((PyArrayObject*)events_arr);

    // Allocate the resulting array of phase values
    PyObject *phase_arr = PyArray_SimpleNew(ndim, dims, dtype);
    double *phase_tmp = (double*)PyArray_DATA((PyArrayObject*)phase_arr);

	// first pass: assign 1 and -1 to peaks and troughs, respectively
	npy_int acc_ = 0;
	npy_int peakhere = 0;
	npy_int troughhere = 0;
	null_der[0] = 0;
	
	for (npy_int j=1; j<arr_len; j++){
		
		peakhere = ((datain[j]>datain[j+1]) & (datain[j]>datain[j-1]));	
		troughhere = -1*((datain[j]<datain[j+1]) & (datain[j]<datain[j-1]));
		null_der[j] = peakhere + troughhere;
		acc_ = acc_ + peakhere - troughhere;

	}

	// allocate an array to store indexes of events and calculate in this way the interval
	npy_intp dims_ev[1] = {acc_};
	 
    PyObject *idxs_arr = PyArray_SimpleNew(ndim, dims_ev, typenum_int);
	int *idx_der = (int*)PyArray_DATA((PyArrayObject*)idxs_arr);

	// second pass: store the indexes of events (peaks and troughs)
	acc_ = 0;
	for (npy_int j=0; j<arr_len; j++){
		if (null_der[j]){
			idx_der[acc_] = j;
			acc_ ++;				
		}		
	}


	// last pass: interpolate phase
	npy_double this_entry = NAN;
	npy_int acc_event = 0;
	npy_int npoints = 0;
	npy_double acc_phase = NAN;

	phase_tmp[0] = NAN;
	for (npy_int j=1; j<arr_len; j++){
		if (null_der[j]){
			this_entry = M_PI*(1-null_der[j])/(2);			
			if (acc_event<acc_){
				npoints = idx_der[acc_event+1]-idx_der[acc_event];
				acc_phase = (M_PI/(double)npoints); 							
				acc_event ++;				
			}
			else {
				acc_phase = NAN;
			}			
		}
		if (null_der[j-1]==-1){
			this_entry = -M_PI;			
		}
		else {
			this_entry = this_entry+acc_phase;			
		}
		phase_tmp[j] = this_entry;		
	}

	return phase_arr;	

}


// Method definitions
static PyMethodDef phase_extraction_methods[] = {
	{"two_point_interp", two_point_interp, METH_VARARGS, "extract instantaneous phase via 2 point interpolation"},
	{NULL, NULL, 0, NULL}
};

// Module defintion
static struct PyModuleDef PhaseExtraction = {
	PyModuleDef_HEAD_INIT,
	"PhaseExtraction", // __name__
	"", // __doc__
	32, // m_size  https://docs.python.org/3/c-api/module.html#c.PyModuleDef
	phase_extraction_methods
};

// Module inititalization
PyMODINIT_FUNC PyInit_PhaseExtraction(void) {
	import_array();
	return PyModuleDef_Init(&PhaseExtraction);
}


