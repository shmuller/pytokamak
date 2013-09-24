import numpy as np
import operator
from utils import memoized_property, BoundingBox
from sig import Signal

from sm_pyplot.tight_figure import get_axes, show

from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline

import dierckx

def atleast_1d_dtype(x, dtype):
    x = np.asanyarray(x, dtype)
    if len(x.shape) == 0:
        x.reshape(1)
    return x


class Spline(object):
    def __init__(self, x, y, k=3, s=0.):
        iopt = 0
        m = x.size
        w = np.ones(m)
        nest = m+k+1
        lwrk = m*(k+1)+nest*(7+3*k)

        t = np.zeros(nest)
        c = np.zeros(nest)
        wrk = np.zeros(lwrk)
        iwrk = np.zeros(nest, np.int32)
        n, fp, ier = dierckx.curfit(iopt, x, y, w, x[0], x[-1], k, s, t, c, wrk, iwrk)
        t.resize(n)
        c.resize(n-k-1)
        self.tck = t, c, k

    @classmethod
    def _from_tck(cls, tck):
        self = cls.__new__(cls)
        self.tck = tck
        return self

    def _op_factory(op):
        def apply(self, other):
            t, c, k = self.tck
            return self._from_tck((t, op(c, other), k))
        return apply

    __add__ = _op_factory(operator.add)
    __sub__ = _op_factory(operator.sub)
    __mul__ = _op_factory(operator.mul)
    __div__ = _op_factory(operator.div)

    def _eval(self, x, nu=0, y=None):
        '''
        Unsafe evaluation: Assume that x is double, sorted and within spline bbox
        '''
        tck = self.tck
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
        x = atleast_1d_dtype(x, np.float64)
        shape = x.shape
        x = x.ravel()
        perm = x.argsort()
        iperm = np.zeros_like(perm)
        iperm[perm] = np.arange(perm.size)
        x = x[perm]
        
        t = self.tck[0]
        i = x.searchsorted(t[[0, -1]])
        if x[-1] == t[-1]:
            i[1] = x.size

        y = np.zeros_like(x)
        y.fill(fill)

        s = slice(i[0], i[1])
        self._eval(x[s], nu, y[s])
        return y[iperm].reshape(shape)

    __call__ = eval

    def deriv(self, nu=1):
        t, c, k = self.tck
        self._eval(t[:1], nu)

        n = t.size - 2*nu
        return self.from_tck((t[nu:n+nu], self.wrk[:n], k - nu))

    def roots(self, out=None, maxroots=100):
        t, c, k = self.tck
        if k == 3:
            ier = 0
            if out is None:
                out = np.zeros(maxroots)
                m = dierckx.sproot(t, c, out, ier)
                out.resize(m)
                return out
            else:
                return dierckx.sproot(t, c, out, ier)
        else:
            raise NotImplementedError("finding roots unsupported for "
                                      "non-cubic splines")

    def as_pp(self):
        t, c, k = self.tck
        x = t[k:-k]
        k1 = k+1
        d = np.zeros((x.size, k1))
        ier = 0
        for xx, dd in zip(x, d):
            dierckx.spalde(t, c, xx, dd, ier)
        return x, d

    def as_signal(self):
        t, c, k = self.tck
        x = t[k:-k]
        y = self._eval(x)
        return Signal(y, x)
        
    def plot(self, ax=None):
        ax = get_axes(ax)
        return self.as_signal().plot(ax=ax)


class SplineOld(InterpolatedUnivariateSpline):
    def __init__(self, *args, **kw):
        data = kw.pop('data', None)
        if data is None:
            InterpolatedUnivariateSpline.__init__(self, *args, **kw)
        else:
            self._data = tuple(data)
            self._reset_class()

    def _op_factory(op):
        def apply(self, other):
            data = list(self._data)
            data[9] = op(data[9], other)
            return self.__class__(data=data)
        return apply

    __add__ = _op_factory(operator.add)
    __sub__ = _op_factory(operator.sub)
    __mul__ = _op_factory(operator.mul)
    __div__ = _op_factory(operator.div)

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
        if x[-1] == t[-1]:
            i[1] = x.size

        y = np.zeros_like(x)
        y.fill(fill)

        s = slice(i[0], i[1])
        self._eval(x[s], nu, y[s])
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

    def roots(self):
        t, c, k = self._eval_args
        if k == 3:
            ier = 0
            z = np.zeros(100)
            m = dierckx.sproot(t, c, z, ier)
            return z[:m]
        raise NotImplementedError("finding roots unsupported for "
                                  "non-cubic splines")

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

    def _op_factory(op):
        def apply(self, other):
            tx, ty, c = self.tck
            data = dict(tck=(tx, ty, op(c, other)), degrees=self.degrees)
            return self.__class__(data=data)
        return apply

    __add__ = _op_factory(operator.add)
    __sub__ = _op_factory(operator.sub)
    __mul__ = _op_factory(operator.mul)
    __div__ = _op_factory(operator.div)

    def astuple(self):
        return self.tck[:3] + self.degrees

    def get_bbox(self):
        tx, ty = self.tck[:2]
        return BoundingBox(np.array(((tx[0], ty[0]), (tx[-1], ty[-1]))))

    def _args(self):
        tx, ty, c = self.tck
        nx, ny = tx.size, ty.size
        kx, ky = self.degrees
        wx, wy = nx-kx-1, ny-ky-1
        ier = 0
        return tx, ty, c, nx, ny, kx, ky, wx, wy, ier

    def eval(self, x, y):
        tx, ty, c, nx, ny, kx, ky, wx, wy, ier = self._args()
        mx, my = x.size, y.size
        lwrk = mx*(kx+1)+my*(ky+1)
        kwrk = mx+my
        wrk = np.zeros(lwrk)
        iwrk = np.zeros(kwrk, 'i')

        z = np.zeros(mx*my)
        dierckx.bispev(tx, ty, c, kx, ky, x, y, z, wrk, iwrk, ier)
        return z.reshape(mx, my)

    def eval_y(self, y):
        tx, ty, c, nx, ny, kx, ky, wx, wy, ier = self._args()
        my = y.size
        C = np.zeros(my*wx)

        dierckx.splevv(ty, c, ky, y, C, wx, ier)
        C = C.reshape((my, wx))
        return [Spline._from_tck((tx, c, kx)) for c in C]

    def eval_x(self, x):
        tx, ty, c, nx, ny, kx, ky, wx, wy, ier = self._args()
        mx = x.size
        C = np.zeros(mx*wy)

        c = c.reshape((wx, wy)).T.ravel().copy()
        dierckx.splevv(tx, c, kx, x, C, wy, ier)
        C = C.reshape((mx, wy))
        return [Spline._from_tck((ty, c, ky)) for c in C]

    def eval_yx(self, x, y):
        tx, ty, c, nx, ny, kx, ky, wx, wy, ier = self._args()
        mx, my = x.size, y.size
        C = np.zeros(my*wx)
        z = np.zeros(mx*my)

        dierckx.splevv(ty, c, ky, y, C, wx, ier)
        dierckx.splevv(tx, C, kx, x, z, my, ier)
        z = z.reshape((mx, my))
        return z, C

    def eval_xy(self, x, y):
        tx, ty, c, nx, ny, kx, ky, wx, wy, ier = self._args()
        mx, my = x.size, y.size
        C = np.zeros(mx*wy)
        z = np.zeros(mx*my)

        c = c.reshape((wx, wy)).T.ravel().copy()
        dierckx.splevv(tx, c, kx, x, C, wy, ier)
        dierckx.splevv(ty, C, ky, y, z, mx, ier)
        z = z.reshape((my, mx)).T
        return z, C

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

    kx, ky = self.degrees
    tx, ty, c = self.tck
    #xi, yi = tx[kx:-kx], ty[kx+1:kx+2]
    xi, yi = np.arange(10.) + 0.5, np.arange(5.) + 0.5
    zi = self.eval(xi, yi)
    
    zi_yx, c_yx = self.eval_yx(xi, yi)
    zi_xy, c_xy = self.eval_xy(xi, yi)
    print np.abs(zi_yx - zi_xy).max()

    spl_y = self.eval_y(y[9:10])[0]
    spl = Spline(x[::9], Z[::9,9])
    print np.abs(spl_y(x) - spl(x)).max()

