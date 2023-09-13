# include <stdint.h>
#include <Python.h>
#include "levenshtein.c"



PyObject *Py_levenshtein(PyObject *self, PyObject *args){
    const char* a;
    const char* b;
    int aN;
    int bN;
    double ins;
    double rm;
    double sub;
    
    if (!PyArg_ParseTuple(args, "yyiiddd", &a, &b, &aN, &bN, &ins, &rm, &sub))
        return NULL;
    
    double res = weighted_levenshtein(a, b, aN, bN, ins, rm, sub);
        
    return PyFloat_FromDouble(res);
}


static PyMethodDef methods[] = {
    {
        "c_levenshtein", 
        Py_levenshtein, 
        METH_VARARGS, 
        "Computes levenshtein distance with specified weights"
     },
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef _levenshtein = {
    PyModuleDef_HEAD_INIT,
    "_levenshtein", 
    "C implementation of levenshtein distance", 
    -1, 
    methods
};


PyMODINIT_FUNC PyInit__levenshtein(){
    return PyModule_Create(&_levenshtein);
};




