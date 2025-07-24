#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "raweep.h"


static PyObject *
raweep_read_block ( PyObject *self, PyObject *args ) {

    const uint8_t * bytes;
    const size_t bsize;
    const uint64_t nsamp, nchan, offset;
    int64_t off;
    int32_t * data;
    PyObject * pydata;
    
    if ( !PyArg_ParseTuple ( args, "y#KKK", &bytes, &bsize, &nsamp, &nchan, &offset ) )
        return NULL;
    
    
    /* Reserves memory for the block data. */
    data = PyMem_Malloc ( nsamp * nchan * sizeof ( *data ) );
    
    /* Reads the block. */
    off  = read_block ( data, bytes, offset, nsamp, nchan );
    
    if ( off < 1 )
        PyErr_SetString ( PyExc_RuntimeError, "The compression method for the first channel cannot be inter-channel residuals (method 3/11)." );
    

    /* Builds a Python tuple from the data stram and the offset. */
    pydata = Py_BuildValue ( "y#K", data, nchan * nsamp * sizeof ( uint32_t ), off );

    /* Destroys the archived data. */
    PyMem_Free ( data );
    
    
    /* Returns the tuple. */
    return pydata;
    
}



/* Defines the method table (available methods). */
static PyMethodDef
raweep_methods [] = {
    { "read_block",  raweep_read_block, METH_VARARGS, "Decodes a block of data from a raw3 stream." },
    { NULL, NULL, 0, NULL }        /* Sentinel */
};

/* Describes the module */
static struct
PyModuleDef raweepmodule = {
    PyModuleDef_HEAD_INIT,
    "raweep", 
    PyDoc_STR ( "Module to read raw data from EEP (CNT) files." ),
    -1,
    raweep_methods
};

/* Module's initialization function. */
PyMODINIT_FUNC
PyInit_raweep ( void ) {
    return PyModule_Create ( &raweepmodule );
}


int
main ( int argc, char *argv [] ) {
    wchar_t *program = Py_DecodeLocale ( argv [0], NULL );
    if ( program == NULL ) {
        fprintf ( stderr, "Fatal error: cannot decode argv[0]\n" );
        exit (1);
    }

    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName ( program );

    /* Add a built-in module, before Py_Initialize */
    if ( PyImport_AppendInittab ( "raweep", PyInit_raweep ) == -1 ) {
        fprintf ( stderr, "Error: could not extend in-built modules table\n" );
        exit (1);
    }

    /* Initializes the Python interpreter. */
    Py_Initialize ();

    /* Imports the module. */
    PyObject *pmodule = PyImport_ImportModule ( "raweep" );
    
    if ( !pmodule ) {
        PyErr_Print ();
        fprintf ( stderr, "Error: could not import module 'raweep'\n" );
    }

    PyMem_RawFree ( program );
    return 0;
}
