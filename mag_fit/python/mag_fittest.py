import numpy as np

from LP import probe_xpr, mag_fit

from sm_pyplot.tight_figure import get_fig, show

XPR = probe_xpr.ProbeXPR(shn=20326)

S = XPR['tip2']

#i0 = 126
i0 = 201

S0 = S[S.V.iE[i0]:S.V.iE[i0+1] + 1]

x, y = S0.V.x, -S0.x

c_params0 = np.array([5e18, 18., 1e-8, 12., 18., 1e-8, 
    180., 1e-6, 0., 0., 0.0009, -1., 2., 1., 2., 0., 1., 0.0016/2, 0., 0.])

do_var = np.array([1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], 'i')

yfit = np.zeros_like(y)

c_params = c_params0.copy()

info = mag_fit.magfit(x, y, yfit, c_params, do_var)

print info

print c_params

fig = get_fig()
fig.axes[0].plot(x, y, x, yfit)
show()

