import numpy as np
from sig import Signal

from sm_pyplot.tight_figure import show

from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline

import dierckx

class Spline(InterpolatedUnivariateSpline):
    def __init__(self, *args, **kw):
        data = kw.pop('data', None)
        if data is None:
            InterpolatedUnivariateSpline.__init__(self, *args, **kw)
        else:
            self._data = tuple(data)
            self._reset_class()

    def eval(self, x, nu=0):
        tck = self._eval_args
        y = np.zeros_like(x)
        self.wrk = np.zeros_like(tck[0])
        ier = 0

        dierckx.splder(*(tck + (nu, x, y, self.wrk, ier)))
        return y

    def deriv(self, nu=1):
        x = self._eval_args[0][:1]
        self.eval(x, nu)
        
        data = list(self._data)
        data[5] -= nu
        data[7] -= 2*nu
        n = data[7]
        data[8] = data[8][nu:n+nu]
        data[9] = self.wrk[:n]
        return Spline(data=data)


class Spline2D(RectBivariateSpline):
    def eval(self, x, y):
        tx, ty, c = self.tck[:3]
        kx, ky = self.degrees
        mx, my = x.size, y.size
        
        lwrk = mx*(kx+1)+my*(ky+1)
        kwrk = mx+my
        self.wrk = np.zeros(lwrk)
        self.iwrk = np.zeros(kwrk, 'i')

        z = np.zeros(mx*my)
        ier = 0
        dierckx.bispev(tx, ty, c, kx, ky, x, y, z, self.wrk, self.iwrk, ier)
        return z.reshape(mx, my)


def spline_test():
    x = np.linspace(0., 10., 100)
    y = np.sin(x)
    #y = np.exp(-x)
    S = Signal(y, x)

    S2 = S[::9]

    spl = Spline(S2.t, S2.x)
    y2 = spl(x)

    dy2 = spl(x, 2)

    #dy2 = np.empty_like(x)
    #for i in xrange(x.size):
    #    dy2[i] = -spl.derivatives(x[i])[1]

    dy3 = spl.eval(x, 2)

    splprime = spl.deriv(nu=2)
    dy4 = splprime.eval(x)


    ax = S.plot()
    #S2.plot(ax=ax)
    #Signal(y2, x).plot(ax=ax)
    Signal(dy2, x).plot(ax=ax)
    Signal(dy3, x).plot(ax=ax)
    Signal(dy4, x).plot(ax=ax)


if __name__ == "__main__":
    #spline_test()
    
    x = np.linspace(0., 10., 100)
    y = np.linspace(1., 2., 100)
    X, Y = np.ix_(x, y)
    Z = np.sin(X)*Y

        
    spl = Spline2D(x[::9], y[::9], Z[::9,::9])

    
    S = Signal(Z[:,::33], x)

    Z2 = spl.eval(x, y)

    ax = S.plot()
    Signal(Z2[:,::33], x).plot(ax=ax)

    show()
