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
    def __init__(self, *args, **kw):
        data = kw.pop('data', None)
        if data is None:
            RectBivariateSpline.__init__(self, *args, **kw)
        else:
            self.tck = data['tck']
            self.degrees = data['degrees']

    def eval(self, x, y):
        tx, ty, c = self.tck[:3]
        kx, ky = self.degrees
        mx, my = x.size, y.size
        
        lwrk = mx*(kx+1)+my*(ky+1)
        kwrk = mx+my
        wrk = np.zeros(lwrk)
        iwrk = np.zeros(kwrk, 'i')

        z = np.zeros(mx*my)
        ier = 0
        dierckx.bispev(tx, ty, c, kx, ky, x, y, z, wrk, iwrk, ier)
        return z.reshape(mx, my)

    def deriv(self, nux=0, nuy=0):
        tx, ty, c = self.tck[:3]
        nx, ny, n = tx.size, ty.size, c.size
        kx, ky = self.degrees
        x, y, z = tx[:1], ty[:1], np.zeros(1)
        
        lwrk = n + (kx+1-nux) + (ky+1-nuy)
        self.wrk = np.zeros(lwrk)
        self.iwrk = np.zeros(2)
        ier = 0

        dierckx.parder(tx, ty, c, kx, ky, nux, nuy, x, y, z, self.wrk, self.iwrk, ier)

        nx, ny = nx-2*nux, ny-2*nuy
        kx, ky = kx-nux, ky-nuy
                
        data = dict(
            degrees = (kx, ky),
            tck = (tx[nux:nx+nux], ty[nuy:ny+nuy], self.wrk[:(nx-kx-1)*(ny-ky-1)].copy()))
        return Spline2D(data=data)


def spline_test(nu=2):
    x = np.linspace(0., 10., 100)
    y = np.sin(x)
    #y = np.exp(-x)
    S = Signal(y, x)

    S2 = S[::9]

    spl = Spline(S2.t, S2.x)
    y2 = spl(x)

    dy2 = spl(x, nu=nu)

    #dy2 = np.empty_like(x)
    #for i in xrange(x.size):
    #    dy2[i] = -spl.derivatives(x[i])[1]

    dy3 = spl.eval(x, nu=nu)

    splprime = spl.deriv(nu=nu)
    dy4 = splprime.eval(x)


    ax = S.plot()
    #S2.plot(ax=ax)
    #Signal(y2, x).plot(ax=ax)
    Signal(dy2, x).plot(ax=ax)
    Signal(dy3, x).plot(ax=ax)
    Signal(dy4, x).plot(ax=ax)


def spline2d_test(nu=2):
    x = np.linspace(0., 10., 100)
    y = np.linspace(0., 10., 100)
    X, Y = np.ix_(x, y)
    Z = np.sin(X)*Y*0.1
    #Z = np.sin(Y)*X*0.1
        
    spl = Spline2D(x[::9], y[::9], Z[::9,::9])
    splx = spl.deriv(nux=nu)
    sply = spl.deriv(nuy=nu)
    
    S = Signal(Z[:,::33], x)
    #S = Signal(Z[::33,:].T, y)

    Z2 = spl.eval(x, y)
    Z2x = splx.eval(x, y)
    Z2y = sply.eval(x, y)

    ax = S.plot()
    
    ax = Signal(Z2[:,::33], x).plot(ax=ax)
    ax = Signal(Z2x[:,::33], x).plot(ax=ax)
    #ax = Signal(Z2[::33,:].T, y).plot(ax=ax)
    #ax = Signal(Z2y[::33,:].T, y).plot(ax=ax)


if __name__ == "__main__":
    nu = 2
    spline_test(nu=nu)
    spline2d_test(nu=nu)
    show()
