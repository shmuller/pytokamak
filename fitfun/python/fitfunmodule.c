#include <Python.h>
#include <numpy/arrayobject.h>

#include <math.h>

#include <minpack.h>

#define template(fun, defs, body, expr)                        \
void fun(data *D) {                                            \
    int i;                                                     \
    double *P = D->P, *x = D->x, *y = D->y;                    \
    double defs;                                               \
    for (i=D->m; i--; ) {                                      \
        body                                                   \
        *y++ = expr;                                           \
    }                                                          \
}

#define template_diff(fun, defs, body, expr)                   \
void fun##_diff(data *D) {                                     \
    int i;                                                     \
    double *P = D->P, *x = D->x, *y = D->y, *ydata = D->ydata; \
    double defs;                                               \
    for (i=D->m; i--; ) {                                      \
        body                                                   \
        *y++ = expr - *ydata++;                                \
    }                                                          \
}

#define template_rms(fun, defs, body, expr)                    \
double fun##_rms(data *D) {                                    \
    int i;                                                     \
    double dy2 = 0., dy;                                       \
    double *P = D->P, *x = D->x, *y = D->y;                    \
    double defs;                                               \
    for (i=D->m; i--; ) {                                      \
        body                                                   \
        dy = expr - *y++;                                      \
        dy2 += dy*dy;                                          \
    }                                                          \
    return sqrt(dy2) / D->m;                                   \
}


#define EMPTY

#define e2_defs P0 = P[0], P1 = P[1]
#define e2_expr P0*exp(-P1*(*x++))

template(e2, e2_defs, EMPTY, e2_expr)
template_diff(e2, e2_defs, EMPTY, e2_expr)
template_rms(e2, e2_defs, EMPTY, e2_expr)


#define IV3_defs P0 = P[0], P1 = P[1], iP2 = 1./P[2]
#define IV3_expr P0*(1. - exp((*x++ - P1)*iP2))

template(IV3, IV3_defs, EMPTY, IV3_expr)
template_diff(IV3, IV3_defs, EMPTY, IV3_expr)
template_rms(IV3, IV3_defs, EMPTY, IV3_expr)


#define IV4_defs *a = D->a, P0 = P[0], P1 = P[1], iP2 = 1./P[2], \
                 dP0 = P[3]-P0, ai, P0i
#define IV4_body ai = *a++; \
                 P0i = P0 + ai*dP0;
#define IV4_expr P0i*(1. - exp((*x++ - P1)*iP2))

template(IV4, IV4_defs, IV4_body, IV4_expr)
template_diff(IV4, IV4_defs, IV4_body, IV4_expr)
template_rms(IV4, IV4_defs, IV4_body, IV4_expr)


#define IV5_defs *a = D->a, P0 = P[0], P1 = P[1], iP2 = 1./P[2], \
                 dP0 = P[3]-P0, dP1 = P[4]-P1, ai, P0i, P1i
#define IV5_body ai = *a++; \
                 P0i = P0 + ai*dP0; \
                 P1i = P1 + ai*dP1;
#define IV5_expr P0i*(1. - exp((*x++ - P1i)*iP2))

template(IV5, IV5_defs, IV5_body, IV5_expr)
template_diff(IV5, IV5_defs, IV5_body, IV5_expr)
template_rms(IV5, IV5_defs, IV5_body, IV5_expr)


#define IV6_defs *a = D->a, P0 = P[0], P1 = P[1], P2 = P[2], \
                 dP0 = P[3]-P0, dP1 = P[4]-P1, dP2 = P[5]-P2, \
                 ai, P0i, P1i, P2i
#define IV6_body ai = *a++; \
                 P0i = P0 + ai*dP0; \
                 P1i = P1 + ai*dP1; \
                 P2i = P2 + ai*dP2;
#define IV6_expr P0i*(1. - exp((*x++ - P1i)/P2i))

template(IV6, IV6_defs, IV6_body, IV6_expr)
template_diff(IV6, IV6_defs, IV6_body, IV6_expr)
template_rms(IV6, IV6_defs, IV6_body, IV6_expr)


double pow_075(double x) 
{
    double tmp;
    if (x >= 0.) return 0.0;

    tmp = sqrt(1. - x - x);
    return (tmp + 2.)*sqrt(tmp - 1.);
}


#define IVdbl_defs P0 = P[0], P1 = P[1], iP2 = 1./P[2], P3 = P[3], P4 = P[4], \
                   A = P3, B = P4, arg, exp_arg
#define IVdbl_body arg = (P1 - *x++)*iP2; exp_arg = exp(arg);
#define IVdbl_expr P0*(exp_arg - 1. - A*pow_075(arg)) / (exp_arg + B)

template(IVdbl, IVdbl_defs, IVdbl_body, IVdbl_expr)
template_diff(IVdbl, IVdbl_defs, IVdbl_body, IVdbl_expr)
template_rms(IVdbl, IVdbl_defs, IVdbl_body, IVdbl_expr)


#define fit_template(name)   \
void name##_fit(data *D) {   \
    leastsq(name##_diff, D); \
}

fit_template(e2)
fit_template(IV3)
fit_template(IV4)
fit_template(IV5)
fit_template(IV6)
fit_template(IVdbl)


#define get_arr(args, i) PyArray_DATA(PyTuple_GET_ITEM(args, i))

void parse_args(PyObject *args, data *D)
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
}

void parse_args_a(PyObject *args, data *D)
{
    parse_args(args, D);
    D->a = get_arr(args, 3);
}

void parse_args_ydata(PyObject *args, data *D)
{
    parse_args(args, D);
    D->ydata = get_arr(args, 3);
}

void parse_args_ydata_a(PyObject *args, data *D)
{
    parse_args_ydata(args, D);
    D->a = get_arr(args, 4);
}


#define meth_template_passthru(fun, parser, i)                \
static PyObject* meth_##fun(PyObject *self, PyObject *args) { \
    data D = {0};                                             \
    parser(args, &D);                                         \
    fun(&D);                                                  \
    PyObject *obj = PyTuple_GET_ITEM(args, i);                \
    Py_INCREF(obj);                                           \
    return obj;                                               \
}

#define meth_template_double(fun, parser)                     \
static PyObject* meth_##fun(PyObject *self, PyObject *args) { \
    data D = {0};                                             \
    parser(args, &D);                                         \
    double d = fun(&D);                                       \
    return Py_BuildValue("d", d);                             \
}

#define meth_template_passthru_fit(fun, parser, i)                          \
static PyObject* meth_##fun(PyObject *self, PyObject *args, PyObject *kw) { \
    PyObject *obj;                                                          \
    data D = {0};                                                           \
    parser(args, &D);                                                       \
    if ((obj = PyDict_GetItemString(kw, "do_var")))                         \
        D.do_var = PyArray_DATA(obj);                                       \
    fun(&D);                                                                \
    obj = PyTuple_GET_ITEM(args, i);                                        \
    Py_INCREF(obj);                                                         \
    return obj;                                                             \
}

meth_template_passthru(e2, parse_args, 2)
meth_template_passthru(e2_diff, parse_args_ydata, 2)
meth_template_passthru_fit(e2_fit, parse_args_ydata, 0)
meth_template_double(e2_rms, parse_args)

meth_template_passthru(IV3, parse_args, 2)
meth_template_passthru(IV3_diff, parse_args_ydata, 2)
meth_template_passthru_fit(IV3_fit, parse_args_ydata, 0)
meth_template_double(IV3_rms, parse_args)

meth_template_passthru(IV4, parse_args_a, 2)
meth_template_passthru(IV4_diff, parse_args_ydata_a, 2)
meth_template_passthru_fit(IV4_fit, parse_args_ydata_a, 0)
meth_template_double(IV4_rms, parse_args_a)

meth_template_passthru(IV5, parse_args_a, 2)
meth_template_passthru(IV5_diff, parse_args_ydata_a, 2)
meth_template_passthru_fit(IV5_fit, parse_args_ydata_a, 0)
meth_template_double(IV5_rms, parse_args_a)

meth_template_passthru(IV6, parse_args_a, 2)
meth_template_passthru(IV6_diff, parse_args_ydata_a, 2)
meth_template_passthru_fit(IV6_fit, parse_args_ydata_a, 0)
meth_template_double(IV6_rms, parse_args_a)

meth_template_passthru(IVdbl, parse_args, 2)
meth_template_passthru(IVdbl_diff, parse_args_ydata, 2)
meth_template_passthru_fit(IVdbl_fit, parse_args_ydata, 0)
meth_template_double(IVdbl_rms, parse_args)


static PyMethodDef methods[] = {
    {"e2", meth_e2, METH_VARARGS, "Exp with 2 parameters"},
    {"e2_diff", meth_e2_diff, METH_VARARGS, "Difference to Exp with 2 parameters"},
    {"e2_rms", meth_e2_rms, METH_VARARGS, "rms for Exp with 2 parameters"},
    {"e2_fit", (PyCFunction) meth_e2_fit, METH_VARARGS|METH_KEYWORDS, 
        "Fit Exp with 2 parameters"},
    {"IV3", meth_IV3, METH_VARARGS, "IV curve with 3 parameters"},
    {"IV3_diff", meth_IV3_diff, METH_VARARGS, "Difference to IV curve with 3 parameters"},
    {"IV3_rms", meth_IV3_rms, METH_VARARGS, "rms for IV curve with 3 parameters"},
    {"IV3_fit", (PyCFunction) meth_IV3_fit, METH_VARARGS|METH_KEYWORDS, 
        "Fit IV curve with 3 parameters"},
    {"IV4", meth_IV4, METH_VARARGS, "IV curve with 4 parameters"},
    {"IV4_diff", meth_IV4_diff, METH_VARARGS, "Difference to IV curve with 4 parameters"},
    {"IV4_rms", meth_IV4_rms, METH_VARARGS, "rms for IV curve with 4 parameters"},
    {"IV4_fit", (PyCFunction) meth_IV4_fit, METH_VARARGS|METH_KEYWORDS, 
        "Fit IV curve with 4 parameters"},
    {"IV5", meth_IV5, METH_VARARGS, "IV curve with 5 parameters"},
    {"IV5_diff", meth_IV5_diff, METH_VARARGS, "Difference to IV curve with 5 parameters"},
    {"IV5_rms", meth_IV5_rms, METH_VARARGS, "rms for IV curve with 5 parameters"},
    {"IV5_fit", (PyCFunction) meth_IV5_fit, METH_VARARGS|METH_KEYWORDS, 
        "Fit IV curve with 5 parameters"},
    {"IV6", meth_IV6, METH_VARARGS, "IV curve with 6 parameters"},
    {"IV6_diff", meth_IV6_diff, METH_VARARGS, "Difference to IV curve with 6 parameters"},
    {"IV6_rms", meth_IV6_rms, METH_VARARGS, "rms for IV curve with 6 parameters"},
    {"IV6_fit", (PyCFunction) meth_IV6_fit, METH_VARARGS|METH_KEYWORDS, 
        "Fit IV curve with 6 parameters"},
    {"IVdbl", meth_IVdbl, METH_VARARGS, "Double probe IV curve"},
    {"IVdbl_diff", meth_IVdbl_diff, METH_VARARGS, "Difference to double probe IV curve"},
    {"IVdbl_rms", meth_IVdbl_rms, METH_VARARGS, "rms for double probe IV curve"},
    {"IVdbl_fit", (PyCFunction) meth_IVdbl_fit, METH_VARARGS|METH_KEYWORDS, 
        "Fit double probe IV curve"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initfitfun(void)
{
    Py_InitModule("fitfun", methods);
}

