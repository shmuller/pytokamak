#include <Python.h>
#include <numpy/arrayobject.h>

//#include "../fitfun.h"
#include <math.h>

typedef struct {
    int n;
    double *P;
    double *x;
    double *y;
} data;


void parse_args(PyObject *args, data *D)
{
	PyObject *p;
	
    p = PyTuple_GET_ITEM(args, 0);
    D->P = PyArray_DATA(p);

    p = PyTuple_GET_ITEM(args, 1);
    D->x = PyArray_DATA(p);
    D->n = PyArray_DIM(p, 0);

    p = PyTuple_GET_ITEM(args, 2);
    D->y = PyArray_DATA(p);
}


void _exp3(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y;
    double a = P[0], b = P[1], c = 1./P[2];

    for (i=D->n; i--; ) {
        *y++ = a*(1. - exp((*x++ - b)*c));
    }
}


double _exp3_diff(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y;
    double a = P[0], b = P[1], c = 1./P[2], cum = 0., tmp;

    for (i=D->n; i--; ) {
        tmp = a*(1. - exp((*x++ - b)*c)) - *y++;
        cum += tmp*tmp;
    }
    return sqrt(cum) / D->n;
}


static PyObject* exp3(PyObject *self, PyObject *args)
{
    data D;
    parse_args(args, &D);
    _exp3(&D);
    Py_RETURN_NONE;
}

static PyObject* exp3_diff(PyObject *self, PyObject *args)
{
    data D;
    parse_args(args, &D);
    double d = _exp3_diff(&D);
    return Py_BuildValue("d", d);
}


static PyMethodDef methods[] = {
    {"exp3", exp3, METH_VARARGS, "Exponential with 3 parameters"},
    {"exp3_diff", exp3_diff, METH_VARARGS, "Difference to exponential with 3 parameters"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initfitfun(void)
{
    Py_InitModule("fitfun", methods);
}
