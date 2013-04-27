import numpy as np
import numpy.ma as ma

from pprint import pformat

import scipy.optimize as opt
import scipy.fftpack

fft = scipy.fftpack.fft
#fft = np.fft.fft

try:
    import bottleneck as bn
    median = bn.median
    movefuns = dict(
        median  = bn.move_median,
        nanmean = bn.move_nanmean,
        nanstd  = bn.move_nanstd,
        nanmin  = bn.move_nanmin,
        nanmax  = bn.move_nanmax)
except ImportError:
    median = np.median
    movefuns = dict(
        median  = None,
        nanmean = None,
        nanstd  = None,
        nanmin  = None,
        nanmax  = None)

try:
    import matplotlib.mlab as mlab
except ImportError:
    pass

from sm_pyplot.tight_figure import get_axes, MathSelector

usetex = False

math_sel = MathSelector(usetex=usetex)


from mediansmooth import *
from cookb_signalsmooth import smooth

from collections import MutableMapping, Iterable, OrderedDict

def ensure_tuple(d, *keys):
    for k in keys:
        if not isinstance(d[k], tuple):
            d[k] = (d[k],)


class memoized_property(object):
    """A read-only @property that is only evaluated once."""
    def __init__(self, fget, doc=None):
        self.fget = fget
        self.__doc__ = doc or fget.__doc__
        self.__name__ = fget.__name__

    def __get__(self, obj, cls):
        if obj is None:
            return self
        obj.__dict__[self.__name__] = result = self.fget(obj)
        return result


class ArrayView(np.ndarray):
    def __new__(subtype, x, fields):
        dtype = {f: x.dtype.fields[f] for f in fields}
        return np.ndarray.__new__(subtype, x.shape, dtype,
                                  buffer=x, strides=x.strides)


class DictView(MutableMapping):
    def __init__(self, base, valid_keys):
        self.base, self.valid_keys = base, valid_keys

    def __getitem__(self, key):
        if key in self.valid_keys:
            return self.base[key]
        else:
            raise KeyError(key)

    def __len__(self):
        return len(self.valid_keys)

    def __iter__(self):
        for key in self.valid_keys:
            yield key

    def __setitem__(self, key, value):
        if key in self.valid_keys:
            self.base[key] = value
        else:
            raise KeyError(key)

    def __delitem__(self, key):
        self.valid_keys.remove(key)

    def __repr__(self):
        return self.__class__.__name__ + " with: " + self.valid_keys.__repr__()


class GeneratorDict(dict):
    def __init__(self, *args, **kw):
        self.generator = kw.pop('generator')
        dict.__init__(self, *args, **kw)

    def __getitem__(self, indx):
        try:
            return dict.__getitem__(self, indx)
        except KeyError:
            self[indx] = result = self.generator(indx)
            return result


class Container(Iterable):
    def __init__(self):
        self.x = OrderedDict()

    def __getitem__(self, indx):
        return self.x[indx]

    def __iter__(self):
        return self.x.itervalues()

    def __add__(self, other):
        s = self.__class__()
        s.x = self.x.copy()
        s.x.update(other.x)
        return s

    def __iadd__(self, other):
        self.x.update(other.x)
        return self

    def __repr__(self):
        return self.__class__.__name__ + " with: " + self.x.keys().__repr__()

    @staticmethod
    def _item(v, attr, cnd):
        x = np.array([getattr(v, attr)])
        if cnd(v):
            return x
        else:
            return np.empty((0,), x.dtype)

    def collect_as_list(self, attr, cnd=lambda v: True): 
        return np.concatenate([v.collect_as_list(attr, cnd) if isinstance(v, Container)
            else self._item(v, attr, cnd) for v in self.x.itervalues()])

    def collect_as_dict(self, attr):
        return {k: getattr(v, attr) for k, v in self.x.iteritems()}


def dict_pop(d, **kw):
    return {k: d.pop(k, kw[k]) for k in kw.keys()}


class NodeInterpolator:
    def __init__(self, n):
        self.n = n

    @memoized_property
    def pascal(self):
        '''Pascal's triangle'''
        n = self.n
        p = np.zeros((n,n), np.int)
        p[:,0] = 1
        oldrow = p[0]
        for i in xrange(1, n):
            row = p[i]
            for j in xrange(1, i+1):
                row[j] = oldrow[j-1] + oldrow[j]
            oldrow = row
        return p

    @memoized_property
    def turned_pascal(self):
        return self.pascal[::-1,::-1].T.copy()

    def calc_c_vect(self, dX, c, out):
        '''Vectorized polynomial evaluation'''
        for k in xrange(c.shape[0]):
            pascal = self.turned_pascal[k]
            o = out[k]
            for m in xrange(k+1):
                o *= dX
                o += pascal[m]*c[m]

    def add_nodes(self, c, x, xi, xmap=None):
        X = np.r_[x, xi]
        x, perm = np.unique(X, return_index=True)

        d, n = c.shape[:2]
        nans = np.zeros((d, perm.size - n))
        nans.fill(np.nan)
        c = np.concatenate((c, nans), axis=1)

        iins = np.flatnonzero(perm >= n+1)
        imap = perm[iins]        
        iins -= np.arange(1, iins.size+1)

        # check
        assert all(iins == np.searchsorted(X[:n+1], X[imap]) - 1)

        if xmap is None:
            dX = X[imap] - X[iins]
        else:
            dX = xmap[X[imap]] - xmap[X[iins]]

        c2 = c[:,iins]
        out = np.zeros_like(c2)

        self.calc_c_vect(dX, c2, out)

        c[:,imap] = out
        
        c = c[:,perm][:,:-1]
        return c, x

    def add_indices(self, c, x, i0, ii):
        return self.add_nodes(c, i0, ii, xmap=x)


class PiecewisePolynomial:
    def __init__(self, c, x, **kw):
        self.c, self.x, self.kw = c, x, kw

        if isinstance(self.c, (tuple, list)):
            self.c = np.concatenate([c[None] for c in self.c], axis=0)

        self.N = self.c.shape[1]
        self.shape = self.c.shape[2:]
        self.bcast = (-1,) + (1,)*len(self.shape)
        self.NI = NodeInterpolator(self.c.shape[0])

        self.fill = kw.get('fill', None)
        self.i0 = kw.get('i0', np.arange(x.size))
        self.i1 = self.i0[1:]

    @memoized_property
    def xi(self):
        return self.x[self.i0]

    def __array__(self):
        return self.c

    def __array_wrap__(self, c):
        return self.__class__(c, self.x, **self.kw)

    def copy(self):
        return self.__array_wrap__(self.c.copy())

    def __add__(self, other):
        c = self.c.copy()
        c[-1] += other
        return self.__array_wrap__(c)

    def __sub__(self, other):
        c = self.c.copy()
        c[-1] -= other
        return self.__array_wrap__(c)

    def __mul__(self, other):
        return self.__array_wrap__(self.c*other)

    def __div__(self, other):
        return self.__array_wrap__(self.c/other)

    def _cmp_factory(op):
        def cmp(self, other):
            i, y = self.eval()
            return np.all(op(y.reshape((-1, 2) + self.shape), other), axis=1)
        return cmp

    __lt__ = _cmp_factory(np.less         )
    __le__ = _cmp_factory(np.less_equal   )
    __eq__ = _cmp_factory(np.equal        )
    __ne__ = _cmp_factory(np.not_equal    )
    __ge__ = _cmp_factory(np.greater_equal)
    __gt__ = _cmp_factory(np.greater      )

    def __getitem__(self, index):
        if not isinstance(index, tuple): 
            index = (index,)
        index = (slice(None), slice(None)) + index
        return self.__class__(self.c[index], self.x, **self.kw)

    def __call__(self, X, **kw):
        X = np.atleast_1d(X)
        ind, outl, outr = self._findind(X, **kw)
        Y = self._polyval(X, ind)

        if self.fill is not None:
            Y[outl | outr] = self.fill
        return Y

    def _findind(self, X, side='right'):
        ind = np.searchsorted(self.xi, X, side) - 1

        indm, indM = 0, self.N - 1
        outl, outr = ind < indm, ind > indM
        ind[outl], ind[outr] = indm, indM
        return ind, outl, outr

    def _polyval(self, X, ind):
        dX = (X - self.xi[ind]).reshape(self.bcast)
        
        c = self.c[:, ind]
        Y = c[0].copy()
        for a in c[1:]:
            Y *= dX
            Y += a
        return Y

    def eval_at_event(self, x_event):
        ind = self._findind(np.array([x_event]))[0]
        s = slice(self.i0[ind], self.i1[ind])
        x = self.x[s]
        y = self._polyval(x, ind.repeat(x.size))
        return s, y

    @memoized_property
    def T(self):
        return self.__class__(self.c.swapaxes(2,3), self.x, **self.kw)

    def add_nodes(self, xi):
        c, x = self.NI.add_nodes(self.c, self.x, xi)
        return self.__class__(c, x, **self.kw)

    def add_indices(self, ii):
        c, self.kw['i0'] = self.NI.add_indices(self.c, self.x, self.i0, ii)
        return self.__class__(c, self.x, **self.kw)

    def _mask(self, w):
        ind0, ind1 = np.searchsorted(self.i0, w)
        return np.concatenate(map(np.arange, ind0, ind1))

    @staticmethod
    def cat(a, axis=0):
        return a[0].__array_wrap__(ma.concatenate(a, axis))

    def _eval_prepare_output(self, il, yl, ir, yr, pad):
        dim = (0, 1)[pad]
        ipad = np.zeros((ir.shape[0], dim), np.int_)
        ypad = np.empty((yr.shape[0], dim) + self.shape, yr.dtype)
        ypad.view('d').fill(np.nan)

        shape = (il.size + ir.size + ipad.size,)
        i = self.cat((il[:,None], ir[:,None], ipad), 1).reshape(shape)
        y = self.cat((yl[:,None], yr[:,None], ypad), 1).reshape(shape + self.shape)
        return i, y

    def eval(self, w=None, pad=False):
        try:
            ind = self._mask(w)[:-1]
            indl, indr = ind, ind + 1
        except:
            indl, indr = slice(0, self.N), slice(1, self.N + 1)

        il = self.i0[indl]
        yl = self.c[-1, indl]
        
        ir = self.i0[indr]
        xr = self.xi[indr]
        yr = self._polyval(xr, indl)
        return self._eval_prepare_output(il, yl, ir, yr, pad)

    @staticmethod
    def _pop_eval_params(kw):
        return dict_pop(kw, w=None, pad=False)
    
    def plot(self, ax=None, x=None, **kw):
        i, y = self.eval(**self._pop_eval_params(kw))

        if x is None:
            x = self.x
        x = x[i]

        ax = get_axes(ax)
        ax.plot(x, y, **kw)
        return ax


class PiecewisePolynomialEndpoints(PiecewisePolynomial):
    def __init__(self, *args, **kw):
        PiecewisePolynomial.__init__(self, *args, **kw)

        self.shift = kw.get('shift', 0)
        self.i1 = kw.get('i1', None)
        if self.i1 is not None:
            if self.i1.size != self.i0.size:
                raise Exception("i1 and i0 must have same size")
            self.i0 = np.concatenate((self.i0, self.i1[-1:]))

    def _shiftind(self, ind, shift):
        if isinstance(ind, slice):
            args = ind.indices(ind.stop)
            ind = np.arange(args[0] + shift, args[1] + shift, args[2])
        else:
            # making a copy is necessary!
            ind = ind + shift
        indm = max(0, shift)
        indM = min(0, shift) + self.N - 1

        outl, outr = ind < indm, ind > indM
        ind[outl], ind[outr] = indm, indM
        return ind, outl, outr

    def _findind(self, X, side='right', shift=None):
        if shift is None:
            shift = self.shift
        ind = np.searchsorted(self.xi, X, side) - 1
        return self._shiftind(ind, shift)

    def eval(self, w=None, pad=False, ext=False, shift=None):
        if shift is None:
            if ext:
                shift = 0
            else:
                shift = self.shift

        try:
            ind = self._mask(w)[:-1]
            indl, indr = ind, ind + 1
        except:
            indl, indr = slice(0, self.N), slice(1, self.N + 1)

        il = self.i0[indl]
        if shift == 0:
            indl_shift = indl
            yl = self.c[-1, indl_shift]
        else:
            indl_shift, outl, outr = self._shiftind(indl, shift)
            xl = self.xi[indl]
            yl = self._polyval(xl, indl_shift)

        if ext:
            ir = self.i1[indl]
            xr = self.x[ir]
        else:
            ir = self.i0[indr]
            xr = self.xi[indr]

        yr = self._polyval(xr, indl_shift)
        return self._eval_prepare_output(il, yl, ir, yr, pad)

    @staticmethod
    def _pop_eval_params(kw):
        return dict_pop(kw, w=None, pad=False, ext=False, shift=None)
        
    def plot_ext(self, *args, **kw):
        return self.plot(*args, pad=True, ext=True, **kw)


class PiecewisePolynomialEndpoints2(PiecewisePolynomial):
    def __init__(self, *args, **kw):
        PiecewisePolynomial.__init__(self, *args, **kw)

        self.shift = kw.get('shift', 0)
        self.i1 = kw.get('i1')
        
        self._add_tail()
        self._shift()
        
    def _add_tail(self):
        ind = np.flatnonzero(self.i1 > self.i0[-1])
        rep = [-1]*(ind.size - 1)
        i0_new = self.i1[ind]
        dX = self.x[i0_new[:-1]] - self.x[self.i0[-1]]

        c = self.c[:,rep]
        c_new = np.zeros_like(c)

        self.NI.calc_c_vect(dX, c, c_new)

        self.i0 = np.concatenate((self.i0, i0_new))
        self.i1 = np.concatenate((self.i1, self.i1[rep]))
        self.c  = np.concatenate((self.c, c_new), axis=1)
        self.N = self.c.shape[1]

    def _shift(self):
        shift = -self.shift
        i0 = self.i0[:-1]
        base = np.concatenate((np.zeros(shift, np.int_), np.arange(self.N-shift)))

        dX = self.x[i0] - self.x[i0[base]]
        c = self.c[:,base]
        self.c = np.zeros_like(c)
        self.NI.calc_c_vect(dX, c, self.c)

    def plot_ext(self, *args, **kw):
        return self.plot(*args, **kw)
        

def PP_test():
    def test_plot(PP, w, ax=None):
        ax = PP.plot(ax=ax)

        PP.plot(ax=ax, w=w, linewidth=2)

        xi = np.linspace(-1., 12., 200)
        yi = PP(xi)
        ax.plot(xi, yi)
        return ax

    x = np.arange(12.)
    c = (np.arange(1., 8.), np.arange(0., 7.))

    i0 = np.arange(8)

    w = np.array([[5], [8]])

    PP = PiecewisePolynomial(c, x, i0=i0)

    ax = test_plot(PP, w)

    i0 = np.arange(7)
    i1 = i0 + 5

    PPE = PiecewisePolynomialEndpoints(c, x, i0=i0, i1=i1, shift=-2)

    ax = None
    ax = PPE.plot_ext(ax=ax)
    ax = PPE.plot_ext(ax=ax, w=w, linewidth=1.5)

    #ax = PPE.plot(ax=ax, shift=0)

    test_plot(PPE, w, ax=ax)

    return PP, PPE


class Amp:
    def __init__(self, fact=1., offs=0., fixpoints=None):
        if fixpoints is not None:
            (x0, y0), (x1, y1) = fixpoints
            self.fact = (y1-y0)/(x1-x0)
            self.offs = y0 - self.fact*x0
        else:
            self.fact, self.offs = fact, offs

    def __call__(self, x):
        return x*self.fact + self.offs   # try x.__mul__ if x is object

    def apply(self, x):
        x *= self.fact
        x += self.offs
        return x

    def __mul__(self, other):
        if isinstance(other, Amp):
            fact = self.fact*other.fact
            offs = self.fact*other.offs + self.offs
            return Amp(fact, offs)
        else:
            return self(other)

    def __imul__(self, other):
        if isinstance(other, Amp):
            self.offs *= other.fact
            self.offs += other.offs
            self.fact *= other.fact
            #self.offs += self.fact*other.offs
            #self.fact *= other.fact
        else:
            self.offs *= other
            self.fact *= other
        return self

    def inv(self):
        ifact = 1./self.fact
        ioffs = -self.offs*ifact
        return Amp(ifact, ioffs)

    def copy(self):
        return Amp(self.fact, self.offs)

    def __repr__(self):
        fmtstr = "%s (fact={fact:.3e}, offs={offs:.3e})"
        return (fmtstr % self.__class__.__name__).format(**self.__dict__)


class Signal:
    def __init__(self, x, t, **kw):
        self.x, self.t, self.kw = x, t, kw
        self.size, self.shape = x.size, x.shape
        self.number = kw.get('number', -1)
        self.name = kw.get('name', "")
        self.type = kw.get('type', None)
        self.units = kw.get('units', "")
        self.tunits = kw.get('tunits', "s")

    def __repr__(self):
        fmtstr = "%s {name} with shape {shape}"
        return (fmtstr % self.__class__.__name__).format(**self.__dict__)

    def __str__(self):
        return self.__repr__() + ", with:\n%s" % pformat(self.__dict__)

    def __call__(self, t):
        t = np.atleast_1d(t)
        return self.__class__(self.PP(t), t, **self.kw)

    def __getitem__(self, indx):
        if not isinstance(indx, tuple):
            indx = (indx,)
        return self.__class__(self.x[indx], self.t[indx[0]], **self.kw)

    def __array__(self):
        return self.x

    def __array_wrap__(self, x):
        return self.__class__(x, self.t, **self.kw)

    def shift_t(self, dt):
        return self.__class__(self.x, self.t + dt, **self.kw)

    def _op_factory(op, calcfun):
        def apply(self, other, out=None):
            try:
                if other.units != self.units:
                    raise Exception("Unit mismatch")
                other = other.x
            except AttributeError:
                pass
            
            return op(self.x, other, out)

        def iapply(self, other):
            apply(self, other, self.x)
            return self

        def wapply(self, other):
            return self.__array_wrap__(apply(self, other))
        
        return locals()[calcfun]

    __lt__   = _op_factory(np.less         , 'apply')
    __le__   = _op_factory(np.less_equal   , 'apply')
    __eq__   = _op_factory(np.equal        , 'apply')
    __ne__   = _op_factory(np.not_equal    , 'apply')
    __ge__   = _op_factory(np.greater_equal, 'apply')
    __gt__   = _op_factory(np.greater      , 'apply')

    __add__  = _op_factory(np.add     , 'wapply')
    __iadd__ = _op_factory(np.add     , 'iapply')
    __sub__  = _op_factory(np.subtract, 'wapply')
    __isub__ = _op_factory(np.subtract, 'iapply')
    __div__  = _op_factory(np.divide  , 'wapply')
    __idiv__ = _op_factory(np.divide  , 'iapply')

    __pow__  = _op_factory(np.power   , 'wapply')

    def __neg__(self):
        return self.__array_wrap__(-self.x)

    def __mul__(self, other):
        return self.__array_wrap__(other*self.x)
    
    def __imul__(self, other):
        try:
            other.apply(self.x)
        except AttributeError:
            np.multiply(self.x, other, self.x)
        return self

    def copy(self):
        return self.__array_wrap__(self.x.copy())

    def astype(self, dtype):
        return self.__array_wrap__(self.x.astype(dtype))

    def masked(self, mask):
        return self.__array_wrap__(ma.masked_array(self.x, mask))

    def unmasked(self):
        if not isinstance(self.x, ma.masked_array):
            return self
        else:
            return self.__array_wrap__(self.x.data)

    def normed(self):
        return self.__array_wrap__(self.x / np.abs(self.range).max())

    def normed_zero_one(self):
        xm, xM = self.range
        return self.__array_wrap__((self.x - xm) / (xM - xm))

    def standardized(self):
        return self.__array_wrap__((self.x - self.x.mean()) / self.x.std())

    def trim(self, s):
        self.x, self.t = self.x[s], self.t[s]
        return self

    @memoized_property
    def range(self):
        return np.nanmin(self.x), np.nanmax(self.x)

    @staticmethod
    def _percentiles(x, p):
        x = np.sort(x)
        i = np.round(np.asarray(p)*(x.size-1)).astype('i')
        return x[i]

    def plot_range(self, r=0.001):
        return tuple(self._percentiles(self.x, [r, 1. - r]))

    def norm_to_region(self, cnd):
        self.x[:] -= self.x[cnd].mean()
        return self

    @memoized_property
    def bcast(self):
        return (-1,) + (1,)*(len(self.shape) - 1)

    @memoized_property
    def PP(self):
        x, t = self.x.astype('d'), self.t.astype('d')
        x0, x1 = x[:-1], x[1:]
        t0, t1 = t[:-1], t[1:]
        dx_dt = (x1 - x0) / (t1 - t0).reshape(self.bcast)
        return PiecewisePolynomial((dx_dt, x0), t)

    def interp(self, ti):
        xi = self.PP(ti)
        return self.__class__(xi, ti, **self.kw)

    def ppolyfit(self, i0, i1, deg):
        N = len(i1)
        c = np.empty((deg + 1, N))
        for i in xrange(N):
            s = slice(i0[i], i1[i])
            t = self.t[s]
            c[:,i] = np.polyfit(t - t[0], self.x[s], deg)
        return c

    def as_PP(self, PP):
        c = self.ppolyfit(PP.i0, PP.i1, PP.c.shape[0] - 1)
        return PP.__class__(c, PP.x, **PP.kw)
        
    def smooth(self, w=100):
        self.x[:] = smooth(self.x, window_len=2*w+1)
        return self

    def mediansmooth(self, w=100):
        mediansmooth(self.x, w)
        return self

    def despike(self, w=2):
        return self.mediansmooth(w)

    def _move_factory(name):
        movefun = movefuns[name]
        def move(self, w=100):
            return self.__array_wrap__(np.roll(movefun(self.x, 2*w+1), -w))
        return move

    move_median = _move_factory('median')
    move_mean = _move_factory('nanmean')
    move_std  = _move_factory('nanstd')
    move_min  = _move_factory('nanmin')
    move_max  = _move_factory('nanmax')

    def remove_median(self, w=100):
        return self - self.move_median(w)

    def rms(self, w=100):
        m = self.move_median(w)
        return np.sqrt(((self - m)**2).move_mean(w))
    
    def rel_rms(self, w=100):
        m = self.move_median(w)
        return np.sqrt(((self - m)**2).move_mean(w)) / m

    def deriv(self, name=""):
        delta = lambda x: np.r_[x[1]-x[0], x[2:]-x[:-2], x[-1]-x[-2]]
        dx_dt = delta(self.x)/delta(self.t)
        return Signal(dx_dt, self.t, name=name, 
                      type="Derivative of " + self.type)

    @staticmethod
    def cross(lvl, x0, x1):
        return (x0 < lvl) & (lvl < x1)
        
    @staticmethod
    def group(ind, threshold):
        if len(ind) == 0:
            ind0 = ind1 = ind
        else:
            di = np.diff(np.r_[0,ind])
            j = np.flatnonzero(di > threshold)
            ind0, ind1 = ind[j], ind[np.r_[j[1:]-1, -1]] + 1
        return ind0, ind1

    def crossings(self, lvl, threshold):
        x0, x1 = self.x[:-1], self.x[1:]
        cnd = self.cross(lvl, x0, x1) | self.cross(lvl, x1, x0)
        ind0, ind1 = self.group(np.flatnonzero(cnd), threshold)

        x0, x1 = self.x[ind0], self.x[ind1]
        cnd = self.cross(lvl, x0, x1) | self.cross(lvl, x1, x0)
        ind0, ind1 = ind0[cnd], ind1[cnd]

        is_rising = self.x[ind0] < self.x[ind1]
        return ind0, ind1, is_rising

    def apply_fun(self, fun, ind0=0, ind1=None):
        slices = map(slice, ind0, ind1)
        res = [fun(self.x[s]) for s in slices]
        return np.array(res)
    
    def apply_argfun(self, argfun, ind0=0, ind1=None):
        if len(ind0) == 0:
            return ind0
        else:
            return ind0 + self.apply_fun(argfun, ind0, ind1)
    
    def local_argmin(self, *args):
        return self.apply_argfun(np.argmin, *args)

    def local_argmax(self, *args):
        return self.apply_argfun(np.argmax, *args)

    @memoized_property
    def xlab(self):
        xlab = "t"
        if len(self.tunits) > 0:
            xlab += " (%s)" % self.tunits
        return xlab

    @memoized_property
    def ylab(self):
        ylab = self.name
        if len(self.units) > 0:
            ylab += " (%s)" % self.units
        return ylab

    def plot(self, ax=None, **kw):
        ax = get_axes(ax, xlab=self.xlab, ylab=self.ylab)
        kw.setdefault('label', self.name)
        ax.plot(self.t, self.x, **kw)
        return ax

    @memoized_property
    def fs(self):
        return (self.t.size - 1) / (self.t[-1] - self.t[0])

    def specgram(self, ax=None, NFFT=2048, step=512, 
            cmap='spectral', **kw):
        Pxx, freqs, bins = mlab.specgram(self.x, NFFT, 1e-3*self.fs,
                mlab.detrend_linear, mlab.window_hanning, NFFT - step)
        
        df = freqs[1] - freqs[0]
        extent = self.t[0], self.t[-1], freqs[0] - df/2, freqs[-1] - df/2
        Z = 10. * np.log10(Pxx)

        ax = get_axes(ax, xlab='t (s)', ylab='f (kHz)')

        im = ax.imshow(Z, aspect='auto', cmap=cmap, origin='lower', 
                extent=extent, interpolation='none', **kw)
        
        ax.figure.tight_layout(pad=0.2)
        return ax

    def xcorr(self, other):
        dt = self.t - self.t[0]
        Dt = np.concatenate((-dt[:0:-1], dt))
        x = y = self.standardized().x
        if other is not self:
            y = other.standardized().x
        C = np.correlate(x / x.size, y, 'full')
        return Signal(C, Dt, tunits=self.tunits, type='Correlation',
                             name='xcorr(%s, %s)' % (self.name, other.name))

    def autocorr(self):
        return self.xcorr(self)

    def lag(self, other, sign=1):
        """
        Returns C, where C.t is the amount of time by which 'other' is lagging 
        behind 'self'
        """
        C = self.xcorr(other)
        return C[np.argmax(sign*C.x)]


class PeriodPhaseFinder:
    def __init__(self, x):
        self.x = x

    def nextpow2(self):
        return 2**np.ceil(np.log2(self.x.size)).astype('i')

    def find_f(self):
        Nfft = self.nextpow2()
        X = fft(self.x-self.x.mean(), Nfft)
        iM = np.abs(X[1:Nfft/2+1]).argmax()+1
        f = np.float(iM)/Nfft
        return f

    def calc_iE(self, p):
        d, i = p
        M = self.x.size
        return np.round(np.arange(i, M-1, d/2)).astype('i')

    def cumdist(self, p):
        iE = self.calc_iE(p)
        dx = np.diff(self.x[iE])
        D = -np.sqrt(dx.dot(dx)/dx.size)
        return D

    @memoized_property
    def guess(self):
        f0 = self.find_f()
        d0 = 1./f0
        x0 = self.x[:np.round(d0)]
        i0 = min(x0.argmin(), x0.argmax())
        p0 = d0, i0
        D0 = self.cumdist(p0)
        return p0, -D0

    @memoized_property
    def guess2(self):
        d0, i0 = self.guess[0]
        N = 200
        dd = np.linspace(0.99*d0, 1.01*d0, N)
        DD = np.empty(N)
        for j in xrange(N):
            DD[j] = self.cumdist((dd[j], i0))
        j = DD.argmin()
        p = dd[j], i0
        D = DD[j]
        return p, -D

    @memoized_property
    def res(self):
        p1 = self.guess2[0]
        p, D = opt.fmin(self.cumdist, p1, full_output=True)[:2]
        return p, -D

    @memoized_property
    def p(self):
        return self.res[0]

    @memoized_property
    def D(self):
        return self.res[1]

    @memoized_property
    def iE(self):
        return self.calc_iE(self.p)

    def print_res(self):
        print self.guess
        print self.guess2
        print self.res


class PositionSignal(Signal):
    def __init__(self, x, t, **kw):
        kw.setdefault('type', 'Position')
        kw.setdefault('units', 'm')

        Signal.__init__(self, x, t, **kw)
    
        self.baseline_slice = kw.get('baseline_slice', slice(None, 1000))
        self.lvl_fact = kw.get('lvl_fact', 0.1)
        self.dist_threshold = kw.get('dist_threshold', 1000)

    def get_baseline(self):
        return self.x[self.baseline_slice]

    def get_crossings(self):
        x = self.get_baseline()
        x0, xM = x.mean(), self.x.max()
        lvl = x0 + self.lvl_fact*(xM - x0)
        return self.crossings(lvl, self.dist_threshold)

    @memoized_property
    def t_ind(self):
        ind0, ind1, is_rising = self.get_crossings()
        i0, i1 = ind0[is_rising], ind1[~is_rising]
        iM = self.local_argmax(i0, i1)
        return i0, iM, i1

    @memoized_property
    def Nplunges(self):
        return len(self.t_ind[1])

    def tM(self, plunge=None):
        i0, iM, i1 = self.t_ind
        if plunge is not None:
            if not isinstance(plunge, Iterable):
                plunge = [plunge]
            iM = iM[plunge]
        return self.t[iM]

    def region_boundaries(self, inout=None):
        i0, iM, i1 = self.t_ind
        if inout == 'in':
            return np.array([i0, iM])
        elif inout == 'out':
            return np.array([iM, i1])
        else:
            return np.array([i0, i1])

    def plunges(self, plunge=None, inout=None):
        w = self.region_boundaries(inout)
        if plunge is not None:
            if not isinstance(plunge, Iterable):
                plunge = [plunge]
            w = w[:,plunge]
        return w

    def regions(self, fun=None, **kw):
        a, b = self.region_boundaries(**kw)
        return map(fun, a, b)

    def get_slices(self, **kw):
        return self.regions(fun=slice, **kw)

    def get_mask(self, **kw):
        return np.concatenate(self.regions(fun=np.arange, **kw))

    def plot_plunges(self, ax=None, **kw):
        ax = Signal.plot(self, ax, **kw)
        i0, iM, i1 = self.t_ind
        ax.plot(self.t[i0], self.x[i0], 'r*')
        ax.plot(self.t[i1], self.x[i1], 'g*')
        ax.plot(self.t[iM], self.x[iM], 'm*')
        return ax


class VoltageSignal(Signal):
    def __init__(self, x, t, **kw):
        kw.setdefault('type', 'Voltage')
        kw.setdefault('units', 'V')

        self.dt = kw.pop('dt', 0.1)
        self.min_ptp = kw.pop('min_ptp', 50)

        Signal.__init__(self, x, t, **kw)

    @memoized_property
    def PPF(self):
        return PeriodPhaseFinder(self.x)

    @memoized_property
    def iE(self): 
        return self.PPF.iE

    @memoized_property
    def is_swept(self):
        D = self.min_ptp
        cnd = self.t < self.t[0] + self.dt
        return self.x[cnd].ptp() > D and self.PPF.D > D

    def plot_sweeps(self, ax=None, **kw):
        ax = Signal.plot(self, ax, **kw)
        ax.plot(self.t[self.iE], self.x[self.iE], 'r+')
        return ax


class CurrentSignal(Signal):
    def __init__(self, x, t, **kw):
        kw.setdefault('type', 'Current')
        kw.setdefault('units', 'A')

        Signal.__init__(self, x, t, **kw)

        self.V = kw.get('V', None)
        self.C = kw.get('C', None)

    def __getitem__(self, index):
        s = Signal.__getitem__(self, index)
        if self.V:
            s.V = self.V[index]
        return s

    def __add__(self, other):
        return self.__class__(self.x + other.x, self.t, V=self.V,
                              name=self.name + '+' + other.name)

    def Isat(self, Vmax=-150):
        x = ma.masked_array(self.x, self.V.x > Vmax)
        return self.__class__(x, self.t, V=self.V, name=self.name)

    def capa_pickup(self):
        self.dV_dt = self.V.copy().smooth(10).deriv()

        cnd = self.t - self.t[0] < 0.05
        dV_dtc = self.dV_dt.x[cnd]
        N = (dV_dtc*dV_dtc).sum()

        dI = self.copy()
        dI.norm_to_region(cnd)
        self.C = (dI.x[cnd]*dV_dtc).sum()/N

    @property
    def I_capa(self):
        if self.C is None:
            self.capa_pickup()
        I_capa = self.C*self.dV_dt.x
        return self.__class__(I_capa, self.t, V=self.V, name=self.name + '_capa')

    @property
    def I_corr(self):
        return self - self.I_capa



