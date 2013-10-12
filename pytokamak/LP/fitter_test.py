import numpy as np
import probe_xpr, probe, LP.fitfun

import scipy.optimize as opt

XPR = probe_xpr.ProbeXPR(shn=28989)

self = XPR.IV_series_simple[3]

iE = self.S.V.iE

def setup(j):
    S = self.S[iE[j]:iE[j+1]]

    fitter_IV = probe.FitterIV(S.V.x, S.x)

    P0 = fitter_IV.get_guess()
    X, Y = fitter_IV.get_norm()
    out = np.empty_like(X)
    return P0, (X, out, Y)


def fit_lsq(j):
    P0, args = setup(j)
    
    def wrapper(P):
        return LP.fitfun.IV3_diff(P, *args).copy()

    P, info = opt.leastsq(wrapper, P0)
    return P

def fit_cus(j):
    P0, args = setup(j)
    P = LP.fitfun.IV3_fit(P0.copy(), *args)
    return P


def fit_compare(j):
    P0, args = setup(j)

    P = LP.fitfun.IV3_fit(P0.copy(), *args)

    def wrapper(P):
        return LP.fitfun.IV3_diff(P, *args).copy()

    P2, info = opt.leastsq(wrapper, P0)
    return P, P2


#j = 2501

for j in xrange(iE.size-1):
    try:
        P, P2 = fit_compare(j)
        print P - P2
    except:
        pass

