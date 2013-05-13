import numpy as np
from sig import memoized_property, Signal

from sm_pyplot.tight_figure import get_axes, show

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

    def _eval(self, x, nu=0, y=None):
        '''
        Unsafe evaluation: Assume that x is sorted and within spline bbox
        '''
        tck = self._eval_args
        if y is None:
            y = np.zeros_like(x)
        self.wrk = np.zeros_like(tck[0])
        ier = 0

        dierckx.splder(*(tck + (nu, x, y, self.wrk, ier)))
        return y

    def eval(self, x, nu=0):
        '''
        Safe evaluation: Check for x values outside bbox and make sure x is
        sorted before passed to _eval().
        '''
        x = np.atleast_1d(x)
        shape = x.shape
        x = x.ravel()
        perm = x.argsort()
        iperm = np.zeros_like(perm)
        iperm[perm] = np.arange(perm.size)
        x = x[perm]
        
        t = self._data[8]
        i = x.searchsorted(t[[0, -1]])
        if x[i[1]-1] == t[-1]:
            i[1] += 1

        y = np.zeros_like(x)
        y.fill(np.nan)

        self._eval(x[i[0]:i[1]], nu, y[i[0]:i[1]])
        return y[iperm].reshape(shape)

    def deriv(self, nu=1):
        x = self._eval_args[0][:1]
        y = self._eval(x, nu)
        
        data = list(self._data)
        data[5] -= nu
        data[7] -= 2*nu
        n = data[7]
        data[8] = data[8][nu:n+nu]
        data[9] = self.wrk[:n]
        return Spline(data=data)

    def plot(self, ax=None):
        ax = get_axes(ax)
        x = self.get_knots()
        y = self._eval(x)
        return Signal(y, x).plot(ax=ax)


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
        wrk = np.zeros(lwrk)
        iwrk = np.zeros(2, 'i')
        ier = 0

        dierckx.parder(tx, ty, c, kx, ky, nux, nuy, x, y, z, wrk, iwrk, ier)

        nx, ny = nx-2*nux, ny-2*nuy
        kx, ky = kx-nux, ky-nuy
                
        data = dict(
            degrees = (kx, ky),
            tck = (tx[nux:nx+nux], ty[nuy:ny+nuy], wrk[:(nx-kx-1)*(ny-ky-1)].copy()))
        return Spline2D(data=data)

    def derivx(self, nu=1):
        return self.deriv(nux=nu)

    def derivy(self, nu=1):
        return self.deriv(nuy=nu)

    def plot(self, ax=None):
        ax = get_axes(ax)
        tx, ty = self.tck[:2]
        kx, ky = self.degrees
        x, y = tx[kx:-kx], ty[ky:-ky]
        z = self.eval(x, y)
        ax.contourf(y, x, z, 30)
        return ax


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
    splprime.plot(ax=ax)
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
    nu = 1
    spline_test(nu=nu)
    spline2d_test(nu=nu)
    show()
