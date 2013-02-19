#include <Python.h>
#include <numpy/arrayobject.h>

#include "../mag_fit.h"


static PyObject* meth_quick_075(PyObject *self, PyObject *args)
{
    PyObject *p[2];
	
	if (!PyArg_ParseTuple(args, "OO", p, p+1)) {
        return NULL;
    }

    int n_x   = PyArray_DIM(p[0], 0);
    double *x = PyArray_DATA(p[0]);
    double *y = PyArray_DATA(p[1]);

    int i;

    init_quick_075();

    for (i=n_x; i--; ) {
        *y++ = quick_075(*x++);
    }

    Py_RETURN_NONE;
}


static PyObject* meth_mag_doppel(PyObject *self, PyObject *args)
{
    PyObject *p[3];
	
	if (!PyArg_ParseTuple(args, "OOO", p, p+1, p+2)) {
        return NULL;
    }

    int n_x      = PyArray_DIM(p[0], 0);

    double *x    = PyArray_DATA(p[0]);
    double *yfit = PyArray_DATA(p[1]);
    double *pars = PyArray_DATA(p[2]);

    mag_doppel(x, yfit, n_x, pars);

    Py_RETURN_NONE;
}


static PyObject* meth_magfit(PyObject *self, PyObject *args)
{
	PyObject *p[5];
	
	if (!PyArg_ParseTuple(args, "OOOOO", p, p+1, p+2, p+3, p+4)) {
        return NULL;
    }
    
    int n_x        = PyArray_DIM(p[0], 0);
    int n_pars     = PyArray_DIM(p[3], 0);

    double *x      = PyArray_DATA(p[0]);
    double *y      = PyArray_DATA(p[1]);
    double *yfit   = PyArray_DATA(p[2]);
    double *pars   = PyArray_DATA(p[3]);
    int    *do_var = PyArray_DATA(p[4]);
    
    int i, rc;
    static int verbose = 0;

    double chi_sqr = 0., eps_abs = 0., eps_rel = 1e-4;
    int iter_max = 200;

    void *mem = calloc((n_x + 2*n_pars)*sizeof(double), 1);

    double *sig = (double*) mem;
    double *nb_values = sig + n_x;
    double *nb_const = nb_values + n_pars;

    for (i=0; i<n_x; ++i) sig[i] = 1.;

    rc = magfit(mag_doppel, x, y, sig, yfit, &n_x, pars, do_var, 
                &n_pars, nb_values, nb_const, 
                &chi_sqr, &iter_max, &eps_abs, &eps_rel, &verbose);

    free(mem);

    return Py_BuildValue("iidd", rc, iter_max, eps_abs, eps_rel);
}


static PyObject* meth_call_magfit(PyObject *self, PyObject *args)
{
	PyObject *p[5];
	
	if (!PyArg_ParseTuple(args, "OOOOO", p, p+1, p+2, p+3, p+4)) {
        return NULL;
    }
    
    int n_x        = PyArray_DIM(p[0], 0);
    int n_params   = PyArray_DIM(p[3], 0);

    double *x      = PyArray_DATA(p[0]);
    double *y      = PyArray_DATA(p[1]);
    double *yfit   = PyArray_DATA(p[2]);
    double *params = PyArray_DATA(p[3]);
    int    *fixed  = PyArray_DATA(p[4]);
    
    int i;

    int rc;
    static int mag_infos = 0;

    double chi_sqr = 0., eps_abs = 0., eps_rel = 1e-4;
    int iter_max = 200;

    int n_pars = 20;
    void *mem = calloc((n_x+3*n_pars)*sizeof(double) + n_pars*sizeof(int), 1);

    double *w = (double*) mem;
    double *c_params = w + n_x;
    double *c_nb_val = c_params + n_pars;
    double *c_nb_cst = c_nb_val + n_pars;
    int *do_var = (int*) c_nb_cst + n_pars;

    for (i=0; i<n_x; ++i) w[i] = 1.;
    
    c_params[0] = params[0];   // n_e
    c_params[1] = params[1];   // T_e
    c_params[2] = params[3];   // alpha1
    c_params[3] = params[4];   // beta
    c_params[4] = params[5];   // d_V
    c_params[5] = params[7];   // alpha2
    c_params[6] = params[2];   // psi1
    c_params[7] = params[6];   // psi2

    // optional parameters
    if (n_params >= 11 && fabs(params[10]) > 1e-8) c_params[17] = params[10];  // p_height
    if (n_params >= 12 && fabs(params[11]) > 1e-8) c_params[10] = params[11];  // a_width
    if (n_params >= 13 && fabs(params[12]) > 1e-8) c_params[11] = params[12];  // a_height
    if (n_params >= 14 && fabs(params[13]) > 1e-8) c_params[12] = params[13];  // a_mass
    if (n_params >= 15 && fabs(params[14]) > 1e-8) c_params[13] = params[14];  // a_z
    if (n_params >= 16 && fabs(params[15]) > 1e-8) c_params[14] = params[15];  // B_t

    for (i=0; i<n_pars; ++i) do_var[i] = 1;

    do_var[0] = fixed[0];
    do_var[1] = fixed[1];
    do_var[2] = fixed[3];
    do_var[3] = fixed[4];
    do_var[4] = fixed[5];
    do_var[5] = fixed[7];
    do_var[6] = fixed[2];
    do_var[7] = fixed[6];

    for (i=0; i<n_pars; ++i) do_var[i] = !do_var[i];

    rc = magfit(mag_doppel, x, y, w, yfit, &n_x,
			    c_params, do_var, &n_pars, c_nb_val, c_nb_cst,
	            &chi_sqr, &iter_max, &eps_abs, &eps_rel, &mag_infos);

    // write back parameters
    params[0] = c_params[0];
    params[1] = c_params[1];
    params[3] = c_params[2];
    params[4] = c_params[3];
    params[5] = c_params[4];
    params[7] = c_params[5];
    params[2] = c_params[6];
    params[6] = c_params[7];

    // if there is enough space, also write back V_wall - V_plasma
    if (n_params > 8) params[8] = c_params[15];

    free(mem);

    Py_RETURN_NONE;
}



static PyMethodDef methods[] = {
    {"quick_075", meth_quick_075, METH_VARARGS, "Tabulated values of x^0.75"},
    {"mag_doppel", meth_mag_doppel, METH_VARARGS, "Fit function"},
    {"magfit", meth_magfit, METH_VARARGS, "Interface to mag_fit library"},
    {"call_magfit", meth_call_magfit, METH_VARARGS, "IDL interface replication"},
    {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initmag_fit(void)
{
    Py_InitModule("mag_fit", methods);
}

