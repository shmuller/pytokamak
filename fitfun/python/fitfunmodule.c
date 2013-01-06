#include <Python.h>
#include <numpy/arrayobject.h>

//#include "../fitfun.h"
#include <math.h>

typedef struct {
    int n;
    double *P;
    double *x;
    double *y;
    double *w;
    double *a;
} data;


void IV3(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y;
    double P0 = P[0], P1 = P[1], iP2 = 1./P[2];

    for (i=D->n; i--; ) {
        *y++ = P0*(1. - exp((*x++ - P1)*iP2));
    }
}

double IV3_diff(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y;
    double P0 = P[0], P1 = P[1], iP2 = 1./P[2];
    double dy2 = 0., dy;

    for (i=D->n; i--; ) {
        dy = P0*(1. - exp((*x++ - P1)*iP2)) - *y++;
        dy2 += dy*dy;
    }
    return sqrt(dy2) / D->n;
}


void IV6(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *a = D->a;
    double P0 = P[3], P1 = P[4], P2 = P[5];
    double dP0 = P[0]-P0, dP1 = P[1]-P1, dP2 = P[2]-P2;
    double ai, P0i, P1i, P2i;

    for (i=D->n; i--; ) {
        ai = *a++;
        P0i = P0 + ai*dP0;
        P1i = P1 + ai*dP1;
        P2i = P2 + ai*dP2;
        *y++ = P0i*(1. - exp((*x++ - P1i)/P2i));
    }
}

double IV6_diff(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *a = D->a;
    double P0 = P[3], P1 = P[4], P2 = P[5];
    double dP0 = P[0]-P0, dP1 = P[1]-P1, dP2 = P[2]-P2;
    double ai, P0i, P1i, P2i;
    double dy2 = 0., dy;

    for (i=D->n; i--; ) {
        ai = *a++;
        P0i = P0 + ai*dP0;
        P1i = P1 + ai*dP1;
        P2i = P2 + ai*dP2;
        dy = P0i*(1. - exp((*x++ - P1i)/P2i)) - *y++;
        dy2 += dy*dy;
    }
    return sqrt(dy2) / D->n;
}



void parse_args(PyObject *args, data *D)
{
    PyObject *obj;
	
    obj = PyTuple_GET_ITEM(args, 0);
    D->P = PyArray_DATA(obj);

    obj = PyTuple_GET_ITEM(args, 1);
    D->x = PyArray_DATA(obj);
    D->n = PyArray_DIM(obj, 0);

    obj = PyTuple_GET_ITEM(args, 2);
    D->y = PyArray_DATA(obj);
}

void parse_args_a(PyObject *args, data *D)
{
    PyObject *obj;

    parse_args(args, D);

    obj = PyTuple_GET_ITEM(args, 3);
    D->a = PyArray_DATA(obj);
}


static PyObject* meth_IV3(PyObject *self, PyObject *args)
{
    data D;
    parse_args(args, &D);
    IV3(&D);
    Py_RETURN_NONE;
}

static PyObject* meth_IV3_diff(PyObject *self, PyObject *args)
{
    data D;
    parse_args(args, &D);
    double d = IV3_diff(&D);
    return Py_BuildValue("d", d);
}

static PyObject* meth_IV6(PyObject *self, PyObject *args)
{
    data D;
    parse_args_a(args, &D);
    IV6(&D);
    Py_RETURN_NONE;
}

static PyObject* meth_IV6_diff(PyObject *self, PyObject *args)
{
    data D;
    parse_args_a(args, &D);
    double d = IV6_diff(&D);
    return Py_BuildValue("d", d);
}



static PyMethodDef methods[] = {
    {"IV3", meth_IV3, METH_VARARGS, "IV curve with 3 parameters"},
    {"IV3_diff", meth_IV3_diff, METH_VARARGS, "Difference to IV curve with 3 parameters"},
    {"IV6", meth_IV6, METH_VARARGS, "IV curve with 6 parameters"},
    {"IV6_diff", meth_IV6_diff, METH_VARARGS, "Difference to IV curve with 6 parameters"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initfitfun(void)
{
    Py_InitModule("fitfun", methods);
}

