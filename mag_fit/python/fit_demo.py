import numpy as np
import matplotlib.pyplot as plt

from LP import mag_fit

def get_case1():
    VI = np.genfromtxt('../idl/demo.dat', skip_header=1)
    V, I = VI[:VI.size/2], VI[VI.size/2:]

    l_sonde = 0.040
    b_sonde = 0.005
    is_lsf = True
    return V, I, l_sonde, b_sonde, is_lsf


def get_case2():
    VI = np.genfromtxt('../idl/lsm_sample.dat', skip_header=0)
    V, I = VI.T

    cnd = V < 85.
    V, I = V[cnd], I[cnd]

    l_sonde = 0.002
    b_sonde = 0.0009
    is_lsf = False
    return V, I, l_sonde, b_sonde, is_lsf


def prepare(get_case):
    x, y, l_sonde, b_sonde, is_lsf = get_case()

    dtor = 0.0174533

    # startwerte fuer Parameter in IDL
    idl_params = np.zeros(17)

    idl_params[0] = 1.e18
    idl_params[1] = 10.
    idl_params[2] = 89.4 * dtor
    idl_params[3] = 1.
    idl_params[4] = 1.
    idl_params[5] = 0.
    idl_params[6] = 89.6 * dtor
    idl_params[7] = 1.

    if is_lsf:
        idl_params[10] = -l_sonde
        idl_params[11] = b_sonde
        idl_params[12] = l_sonde
    else:
        idl_params[2] = 10. * dtor
        idl_params[3] = 0.05
        idl_params[10] = -l_sonde
        idl_params[11] = b_sonde
        idl_params[12] = l_sonde / np.cos(idl_params[2])


    idl_fixed = np.ones(17, 'i')

    idl_fixed[0] = 0	# 0 --> Fit in entsprechendem Parameter
    idl_fixed[1] = 0
    # kein Fit im Winkel: idl_fixed[2] = 0
    idl_fixed[3] = 0
    idl_fixed[4] = 0
    idl_fixed[5] = 0
    idl_fixed[7] = 0

    return x, y, idl_params, idl_fixed


def fit(get_case):
    x, y, idl_params, idl_fixed = prepare(get_case)
    yfit = np.zeros_like(y)

    rc = mag_fit.call_magfit(x, y, yfit, idl_params, idl_fixed)

    print idl_params

    plt.plot(x, y, '-+', x, yfit)


fit(get_case1)
fit(get_case2)

plt.show()

