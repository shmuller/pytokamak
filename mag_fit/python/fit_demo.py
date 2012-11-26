import numpy as np
import matplotlib.pyplot as plt

from LP import mag_fit

fname = '../idl/lsm_sample.dat'

VI = np.genfromtxt('../idl/lsm_sample.dat', skip_header=0)

V, I = VI.T

cnd = V < 85
V, I = V[cnd], I[cnd]

l_sonde = 0.002
b_sonde = 0.0009

is_lsf = False

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

x, y = V, I
yfit = np.zeros_like(y)

rc = mag_fit.call_magfit( x, y, yfit, idl_params, idl_fixed)

print idl_params

plt.plot(V, I, '-+', x, yfit)

plt.show()


