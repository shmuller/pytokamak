#include <Python.h>
#include <numpy/arrayobject.h>

//#include "../fitfun.h"
#include <math.h>

#include <minpack.h>


void IV3(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y;
    double P0 = P[0], P1 = P[1], iP2 = 1./P[2];

    for (i=D->m; i--; ) {
        *y++ = P0*(1. - exp((*x++ - P1)*iP2));
    }
}

void IV3_diff(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *ydata = D->ydata;
    double P0 = P[0], P1 = P[1], iP2 = 1./P[2];

    for (i=D->m; i--; ) {
        *y++ = P0*(1. - exp((*x++ - P1)*iP2)) - *ydata++;
    }
}

double IV3_rms(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y;
    double P0 = P[0], P1 = P[1], iP2 = 1./P[2];
    double dy2 = 0., dy;

    for (i=D->m; i--; ) {
        dy = P0*(1. - exp((*x++ - P1)*iP2)) - *y++;
        dy2 += dy*dy;
    }
    return sqrt(dy2) / D->m;
}

void IV4(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *a = D->a;
    double P0 = P[0], P1 = P[1], iP2 = 1./P[2];
    double dP0 = P[3]-P0;
    double ai, P0i;

    for (i=D->m; i--; ) {
        ai = *a++;
        P0i = P0 + ai*dP0;
        *y++ = P0i*(1. - exp((*x++ - P1)*iP2));
    }
}

void IV4_diff(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *a = D->a, *ydata = D->ydata;
    double P0 = P[0], P1 = P[1], iP2 = 1./P[2];
    double dP0 = P[3]-P0;
    double ai, P0i;

    for (i=D->m; i--; ) {
        ai = *a++;
        P0i = P0 + ai*dP0;
        *y++ = P0i*(1. - exp((*x++ - P1)*iP2)) - *ydata++;
    }
}

double IV4_rms(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *a = D->a;
    double P0 = P[0], P1 = P[1], iP2 = 1./P[2];
    double dP0 = P[3]-P0;
    double ai, P0i;
    double dy2 = 0., dy;

    for (i=D->m; i--; ) {
        ai = *a++;
        P0i = P0 + ai*dP0;
        dy = P0i*(1. - exp((*x++ - P1)*iP2)) - *y++;
        dy2 += dy*dy;
    }
    return sqrt(dy2) / D->m;
}

void IV6(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *a = D->a;
    double P0 = P[0], P1 = P[1], P2 = P[2];
    double dP0 = P[3]-P0, dP1 = P[4]-P1, dP2 = P[5]-P2;
    double ai, P0i, P1i, P2i;

    for (i=D->m; i--; ) {
        ai = *a++;
        P0i = P0 + ai*dP0;
        P1i = P1 + ai*dP1;
        P2i = P2 + ai*dP2;
        *y++ = P0i*(1. - exp((*x++ - P1i)/P2i));
    }
}

void IV6_diff(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *a = D->a, *ydata = D->ydata;
    double P0 = P[0], P1 = P[1], P2 = P[2];
    double dP0 = P[3]-P0, dP1 = P[4]-P1, dP2 = P[5]-P2;
    double ai, P0i, P1i, P2i;

    for (i=D->m; i--; ) {
        ai = *a++;
        P0i = P0 + ai*dP0;
        P1i = P1 + ai*dP1;
        P2i = P2 + ai*dP2;
        *y++ = P0i*(1. - exp((*x++ - P1i)/P2i)) - *ydata++;
    }
}

double IV6_rms(data *D)
{
    int i;
    double *P = D->P, *x = D->x, *y = D->y, *a = D->a;
    double P0 = P[0], P1 = P[1], P2 = P[2];
    double dP0 = P[3]-P0, dP1 = P[4]-P1, dP2 = P[5]-P2;
    double ai, P0i, P1i, P2i;
    double dy2 = 0., dy;

    for (i=D->m; i--; ) {
        ai = *a++;
        P0i = P0 + ai*dP0;
        P1i = P1 + ai*dP1;
        P2i = P2 + ai*dP2;
        dy = P0i*(1. - exp((*x++ - P1i)/P2i)) - *y++;
        dy2 += dy*dy;
    }
    return sqrt(dy2) / D->m;
}



static PyObject* parse_args(PyObject *args, data *D)
{
    PyObject *obj;
	
    obj = PyTuple_GET_ITEM(args, 0);
    D->P = PyArray_DATA(obj);
    D->n = PyArray_DIM(obj, 0);

    obj = PyTuple_GET_ITEM(args, 1);
    D->x = PyArray_DATA(obj);
    D->m = PyArray_DIM(obj, 0);

    obj = PyTuple_GET_ITEM(args, 2);
    D->y = PyArray_DATA(obj);
    return obj;
}

static PyObject* parse_args_a(PyObject *args, data *D)
{
    PyObject *obj;

    obj = PyTuple_GET_ITEM(args, 3);
    D->a = PyArray_DATA(obj);

    return parse_args(args, D);
}

static PyObject* parse_args_ydata(PyObject *args, data *D)
{
    PyObject *obj;

    obj = PyTuple_GET_ITEM(args, 3);
    D->ydata = PyArray_DATA(obj);

    /*
    int i, m = PyArray_DIM(obj, 0);
    double *ydata = D->ydata;

    printf("m = %d\n", m);
    for (i=0; i < m; ++i)
        printf("%f\n", *ydata++);
    */

    return parse_args(args, D);
}

static PyObject* parse_args_ydata_a(PyObject *args, data *D)
{
    PyObject *obj;

    obj = PyTuple_GET_ITEM(args, 4);
    D->a = PyArray_DATA(obj);

    return parse_args_ydata(args, D);
}


#define meth_template_passthru(fun, parser)                   \
static PyObject* meth_##fun(PyObject *self, PyObject *args) { \
    data D;                                                   \
    PyObject *obj = parser(args, &D);                         \
    fun(&D);                                                  \
    Py_INCREF(obj);                                           \
    return obj;                                               \
}

#define meth_template_double(fun, parser)                     \
static PyObject* meth_##fun(PyObject *self, PyObject *args) { \
    data D;                                                   \
    parser(args, &D);                                         \
    double d = fun(&D);                                       \
    return Py_BuildValue("d", d);                             \
}

meth_template_passthru(IV3, parse_args)
meth_template_passthru(IV3_diff, parse_args_ydata)
meth_template_double(IV3_rms, parse_args)

meth_template_passthru(IV4, parse_args_a)
meth_template_passthru(IV4_diff, parse_args_ydata_a)
meth_template_double(IV4_rms, parse_args_a)

meth_template_passthru(IV6, parse_args_a)
meth_template_passthru(IV6_diff, parse_args_ydata_a)
meth_template_double(IV6_rms, parse_args_a)


static PyMethodDef methods[] = {
    {"IV3", meth_IV3, METH_VARARGS, "IV curve with 3 parameters"},
    {"IV3_diff", meth_IV3_diff, METH_VARARGS, "Difference to IV curve with 3 parameters"},
    {"IV3_rms", meth_IV3_rms, METH_VARARGS, "rms for IV curve with 3 parameters"},
    {"IV4", meth_IV4, METH_VARARGS, "IV curve with 4 parameters"},
    {"IV4_diff", meth_IV4_diff, METH_VARARGS, "Difference to IV curve with 4 parameters"},
    {"IV4_rms", meth_IV4_rms, METH_VARARGS, "rms for IV curve with 4 parameters"},
    {"IV6", meth_IV6, METH_VARARGS, "IV curve with 6 parameters"},
    {"IV6_diff", meth_IV6_diff, METH_VARARGS, "Difference to IV curve with 6 parameters"},
    {"IV6_rms", meth_IV6_rms, METH_VARARGS, "rms for IV curve with 6 parameters"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initfitfun(void)
{
    Py_InitModule("fitfun", methods);
}

