import numpy as np
from utils import memoized_property, BoundingBox
from sig import Signal

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

    def get_bbox(self):
        x = self.get_knots()
        return BoundingBox(x[0], x[-1])

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

    def eval(self, x, nu=0, fill=np.nan):
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
        y.fill(fill)

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

    def assignal(self):
        x = self.get_knots()
        y = self._eval(x)
        return Signal(y, x)
        
    def plot(self, ax=None):
        ax = get_axes(ax)
        return self.assignal().plot(ax=ax)


class Spline2D(RectBivariateSpline):
    def __init__(self, *args, **kw):
        data = kw.pop('data', None)
        if data is None:
            RectBivariateSpline.__init__(self, *args, **kw)
        else:
            self.tck = data['tck']
            self.degrees = data['degrees']

        # predefine workspaces for fast 1-point evaluation
        kx, ky = self.degrees
        self.wrk1 = np.zeros(kx + ky + 2)
        self.iwrk1 = np.zeros(2, 'i')
        self.z1 = np.zeros(1)

    def astuple(self):
        return self.tck[:3] + self.degrees

    def get_bbox(self):
        tx, ty = self.tck[:2]
        return BoundingBox(np.array(((tx[0], ty[0]), (tx[-1], ty[-1]))))

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

    def eval_test(self, x, y):
        tx, ty, c = self.tck[:3]
        kx, ky = self.degrees
        mx, my = x.size, y.size
        
        lwrk = mx*(kx+1)+my*(ky+1)
        kwrk = mx+my
        wrk = np.zeros(lwrk)
        iwrk = np.zeros(kwrk, 'i')

        wx = np.zeros((mx, kx+1), order='F')
        wy = np.zeros((my, ky+1), order='F')
        lx = np.zeros(mx, 'i')
        ly = np.zeros(my, 'i')

        z = np.zeros(mx*my)

        dierckx.fpbisp(tx, ty, c, x, y, z, wx, wy, lx, ly)
        return z.reshape(mx, my), wx, wy

    def eval_tensor(self, x, y):
        tx, ty, c = self.tck[:3]
        kx, ky = self.degrees
        nx, ny = tx.size, ty.size
        mx, my = x.size, y.size
        C = np.zeros((nx-kx-1, my))

        c = c.reshape((nx-kx-1, ny-ky-1))
        ier = 0
        wrk = np.zeros(ny)

        for i in xrange(nx-kx-1):
            dierckx.splev(ty, c[i], ky, y, C[i], ier)

        z = np.zeros((my, mx))
        C = C.T.copy()
        for j in xrange(my):
            dierckx.splev(tx, C[j], kx, x, z[j], ier)
        return z.T, C

    def eval1(self, x, y):
        # shortcut for fast evaluation at 1 point
        tx, ty, c = self.tck[:3]
        kx, ky = self.degrees
        z = self.z1
        dierckx.bispev(tx, ty, c, kx, ky, x, y, z, self.wrk1, self.iwrk1, 0)
        return z

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
    #nu = 1
    #spline_test(nu=nu)
    #spline2d_test(nu=nu)
    #show()
    x = y = np.linspace(0., 10., 100)
    X, Y = np.ix_(x, y)
    Z = np.sin(X)*Y*0.1
    self = Spline2D(x[::9], y[::9], Z[::9,::9])

    spl = Spline(x[::9], Z[::9,9])
    t, c, k = spl._eval_args

    kx, ky = self.degrees
    tx, ty, C = self.tck
    #xi, yi = tx[kx:-kx], ty[kx+1:kx+2]
    xi, yi = np.arange(10.) + 0.5, np.arange(5.) + 0.5
    zi = self.eval(xi, yi)
    zi2, wx, wy = self.eval_test(xi, yi)
    
    zi3, CC = self.eval_tensor(xi, yi)
