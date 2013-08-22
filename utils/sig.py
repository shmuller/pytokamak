import numpy as np
import numpy.ma as ma

import operator
from math import sqrt

from pprint import pformat

import scipy.optimize as opt
import scipy.fftpack
import scipy.signal

fft = scipy.fftpack.fft
fftconvolve = scipy.signal.fftconvolve

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
    import matplotlib.colors as colors
except ImportError:
    pass

from sm_pyplot.tight_figure import get_axes, MathSelector

usetex = False

math_sel = MathSelector(usetex=usetex)

from mediansmooth import *
from cookb_signalsmooth import smooth
from utils import memoized_property, dict_pop, DictView


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

    def __call__(self, X, masked=False, **kw):
        X = np.atleast_1d(X)
        ind, outl, outr = self._findind(X, **kw)
        Y = self._polyval(X, ind)

        if self.fill is not None:
            Y[outl | outr] = self.fill
        if masked:
            Y = ma.masked_array(Y, outl | outr)
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


def _op_factory(op, neutral):
    def _op(x, other):
        if np.isscalar(other):
            if other == neutral:
                return x
            elif np.isscalar(x):
                return op(x, other)
        raise TypeError("cannot apply operation")
    return _op

add = _op_factory(operator.add, 0)
sub = _op_factory(operator.sub, 0)
mul = _op_factory(operator.mul, 1)
div = _op_factory(operator.div, 1)

class Amp:
    def __init__(self, fact=1, offs=0, fixpoints=None):
        if fixpoints is not None:
            (x0, y0), (x1, y1) = fixpoints
            self.fact = (y1 - y0) / (x1 - x0)
            self.offs = y0 - self.fact * x0
        else:
            self.fact, self.offs = fact, offs

    def __call__(self, x):
        return x * self.fact + self.offs   # try x.__mul__ if x is object

    def __getitem__(self, indx):
        try:
            fact = self.fact[indx]
        except TypeError:
            fact = self.fact
        try:
            offs = self.offs[indx]
        except TypeError:
            offs = self.offs
        return self.__class__(fact, offs)

    def apply(self, x):
        x *= self.fact
        x += self.offs
        return x

    def _add_factory(op):
        def wapply(self, other):
            return self.__class__(self.fact, op(self.offs, other))

        def iapply(self, other):
            self.offs = op(self.offs, other)
            return self
        return wapply, iapply

    def _mul_factory(op):
        def wapply(self, other):
            return self.__class__(op(self.fact, other), op(self.offs, other))

        def iapply(self, other):
            self.fact = op(self.fact, other)
            self.offs = op(self.offs, other)
            return self
        return wapply, iapply

    __add__, __iadd__ = _add_factory(add)
    __sub__, __isub__ = _add_factory(sub)
    __mul  , __imul   = _mul_factory(mul)
    __div__, __idiv__ = _mul_factory(div)

    def __mul__(self, other):
        if isinstance(other, Amp):
            fact = self.fact * other.fact
            offs = self.fact * other.offs + self.offs
            return self.__class__(fact, offs)
        else:
            return self.__mul(other)
            
    def __imul__(self, other):
        if isinstance(other, Amp):
            self.offs *= other.fact
            self.offs += other.offs
            self.fact *= other.fact
            return self
        else:
            return self.__imul(other)

    def inv(self):
        ifact = 1. / self.fact
        ioffs = -self.offs * ifact
        return self.__class__(ifact, ioffs)

    def copy(self):
        return self.__class__(self.fact, self.offs)

    def __repr__(self):
        fmtstr = "%s (fact={fact}, offs={offs})"
        return (fmtstr % self.__class__.__name__).format(**self.__dict__)


class Detrender:
    def __init__(self, N, detrend='linear'):
        xm = (N - 1) / 2.
        dx = 1. / N
        fact = dx * sqrt(12.*dx/(1 - dx*dx))
        self.x = np.arange(-xm, xm + 1) * fact

        detrenders = dict(none=self.detrend_none,
                          mean=self.detrend_mean, 
                          linear=self.detrend_linear)
        self.detrend = detrenders[detrend]

    def detrend_none(self, y):
        return y

    def detrend_mean(self, y):
        return y - y.mean(axis=0)

    def detrend_linear(self, y):
        x = self.x
        y = y - y.mean(axis=0)
        k = np.dot(x, y)
        if k.ndim == 0:
            return y - k * x
        else:
            return y - k[None] * x[:, None]


class Window:
    def __init__(self, N, window='hanning'):
        windows = dict(none=np.ones,
                       hanning=np.hanning)

        self.window = windows[window](N)
        self.norm = 1. / np.sum(self.window**2)


class Signal2D:
    cdict = dict(blue =((0,1,1),(.5,0,0),(.75,0,0),(1,0,0)), 
                 green=((0,0,0),(.5,1,1),(.75,1,1),(1,0,0)), 
                 red  =((0,0,0),(.5,1,1),(.75,0,0),(1,1,1)))

    cmap = colors.LinearSegmentedColormap('spectral_colormap', cdict, 256)

    #clist = [(0,0,1)]*4 + [(1,1,0)]*3 + [(0,1,0)]*2 + [(1,0,0)]*1
    #cmap = colors.ListedColormap(clist)

    def __init__(self, x, y, Z, **kw):
        self.x, self.y, self.Z, self.kw = x, y, Z, kw
        self.type = kw.get('type', '')
        self.units = kw.get('units', '')
        self.xtype = kw.get('xtype', 't')
        self.ytype = kw.get('ytype', '')
        self.xunits = kw.get('xunits', 's')
        self.yunits = kw.get('yunits', '')

    def __array_wrap__(self, x, y, Z):
        return self.__class__(x, y, Z, **self.kw)

    def __getitem__(self, indx):
        if not isinstance(indx, tuple):
            indx = (indx,)
        if len(indx) == 1:
            indx = indx + (slice(None),)
        x, y, Z = self.x[indx[0]], self.y[indx[1]], self.Z[indx[::-1]]
        if y.ndim == 0:
            return Signal(Z, x, type=self.type, units=self.units, tunits=self.xunits)
        elif x.ndim == 0:
            return Signal(Z, y, type=self.type, units=self.units, tunits=self.yunits)
        else:
            return self.__array_wrap(x, y, Z)

    @memoized_property
    def xlab(self):
        xlab = self.xtype
        if len(self.xunits) > 0:
            xlab += " (%s)" % self.xunits
        return xlab

    @memoized_property
    def ylab(self):
        ylab = self.ytype
        if len(self.yunits) > 0:
            ylab += " (%s)" % self.yunits
        return ylab

    def plot(self, ax=None, ylim=None, cmap='spectral', **kw):
        if cmap == 'custom':
            cmap = self.cmap
        x, y, Z = self.x, self.y, self.Z
        
        if ylim is not None:
            cnd = (ylim[0] <= y) & (y <= ylim[1])
            y, Z = y[cnd], Z[cnd]
        
        dy = y[1] - y[0]
        y = np.r_[y[0], y[1:] - dy/2, y[-1]]
        
        kw.setdefault('xlab', self.xlab)
        kw.setdefault('ylab', self.ylab)
        ax = get_axes(ax, **kw)
        im = ax.pcolorfast(x, y, Z, cmap=cmap)
        ax._current_image = im
        return ax

    def plot_with_signal(self, S, ax=None, yloc=None, *args, **kw):
        ax = self.plot(ax, *args, **kw)
        ax2 = ax.get_2nd_axes(ylab=S.ylab, yloc=yloc)
        S.plot(ax2)
        ax.figure.tight_layout()
        return ax


class Windower:
    def __init__(self, w=2048, step=512, detrend='linear', window='hanning'):
        self.w, self.step = w, step
        self.detrender = Detrender(w, detrend=detrend)
        self.win = Window(w, window=window)
            
    def __call__(self, S):
        return self.windowed_calc(S)

    def calc(self, xi):
        raise NotImplementedError

    def prepare_output(self, t, S, Z):
        raise NotImplementedError

    def windowed_calc(self, S):
        w, step = self.w, self.step

        S = S.filled()
        x = S.x

        window = self.win.window.reshape((-1,) + (1,)*(x.ndim - 1))
        detrend = self.detrender.detrend
        calc = self.calc

        i0 = np.arange(0, x.shape[0] - w + 1, step)
        i1 = i0 + w

        Z = np.empty((self.M, i0.size))
        for i in xrange(i0.size):
            xi = x[i0[i]:i1[i]]
            Z[:,i] = calc(window * detrend(xi))
                
        n = i0.size - 1
        ind0 = (w + step) // 2
        ind = np.r_[0, np.arange(ind0, ind0 + n*step, step), n*step + w - 1]
        t = S.t[ind]
        
        return self.prepare_output(t, S, Z)


class Specgram(Windower):
    def __init__(self, *args, **kw):
        Windower.__init__(self, *args, **kw)
        self.M = self.w // 2 + 1

    def calc(self, xi):
        Xi = fft(xi)[:self.M]
        return Xi.real**2 + Xi.imag**2

    def prepare_output(self, t, S, Z):
        fs = S.fs
        f = np.arange(self.M) * (fs / self.w)

        Z *= self.win.norm
        Z[1:-1] *= 2.
        Z = 10. * np.log10(Z)
        funits = dict(s='Hz', ms='kHz')[S.tunits]
        return Signal2D(t, f, Z, type='Spectral intensity', units='dB',
                        xtype='t', xunits=S.tunits, ytype='f', yunits=funits)


class XCorr(Windower):
    def __init__(self, *args, **kw):
        kw.setdefault('window', 'none')
        Windower.__init__(self, *args, **kw)
        self.M = self.w * 2 - 1
        if self.w < 500:
            self._correlate = self._correlate_direct
        else:
            self._correlate = self._correlate_fft

    def _correlate_direct(self, x, y):
        return np.correlate(x / x.size, y, 'full')

    def _correlate_fft(self, x, y):
        return fftconvolve(x / x.size, y[::-1], 'full')

    def calc(self, xi):
        xi = (xi - xi.mean(axis=0)) / xi.std(axis=0)
        if xi.ndim == 1:
            return self._correlate(xi, xi)
        else:
            return self._correlate(xi[:, 0], xi[:, 1])

    def prepare_output(self, t, S, C): 
        dt = 1. / S.fs
        Dt = np.arange(-(self.w - 1), self.w) * dt
        return Signal2D(t, Dt, C, type='Correlation',
                        xtype='t', xunits=S.tunits, ytype='$\Delta$t', yunits=S.tunits)


class Filter:
    def __init__(self, fm=None, fM=None, f_Nyq=1., order=4):
        from scipy import signal

        if fm is None:
            b, a = signal.butter(order, fM/f_Nyq, 'lowpass')
        elif fM is None:
            b, a = signal.butter(order, fm/f_Nyq, 'highpass')
        else:
            b, a = signal.butter(order, [fm/f_Nyq, fM/f_Nyq], 'bandpass')

        def filter(x):
            return signal.filtfilt(b, a, x)
        self.filter = filter

    def __call__(self, x):
        mask = ma.getmask(x)
        if mask is ma.nomask:
            return self.filter(x)
        else:
            return ma.masked_array(self.filter(x.filled(0.)), mask)


class SignalBase:
    '''Signal base class - Makes no assumptions on x and t
    '''
    fmtstr = "%s with shape {shape}"

    def __init__(self, x, t=None, *args, **kw):
        self._x = np.atleast_1d(x)
        if t is None:
            t = np.arange(x.shape[0])
        self._t = np.atleast_1d(t)
        self.size, self.shape = self._x.size, self._x.shape
        self.args, self.kw = args, kw
        self._result_class = self.__class__

    @property
    def x(self):
        return self._x

    @property
    def t(self):
        return self._t

    def update(self, **kw):
        self.kw.update(kw)
        for k, v in kw.iteritems():
            setattr(self, k, v)

    def __repr__(self):
        return (self.fmtstr % self.__class__.__name__).format(**self.__dict__)

    def __str__(self):
        return self.__repr__() + ", with:\n%s" % pformat(self.__dict__)

    def __call__(self, t, *args, **kw):
        return self.interp(t, *args, **kw)

    def __getitem__(self, indx):
        if not isinstance(indx, tuple):
            indx = (indx,)
        return self.__class__(self._x[indx], self._t[indx[0]], *self.args, **self.kw)

    def __array__(self):
        return self.x

    def _wrap(self, x):
        return self.__class__(x, self._t, *self.args, **self.kw)

    def __array_wrap__(self, x):
        return self._result_class(x, self._t, *self.args, **self.kw)

    def _cmp_factory(op):
        def cmp(self, other, out=None):
            if isinstance(other, SignalBase):
                other = other.x
            return op(self.x, other, out)
        return cmp

    __lt__ = _cmp_factory(np.less)
    __le__ = _cmp_factory(np.less_equal)
    __eq__ = _cmp_factory(np.equal)
    __ne__ = _cmp_factory(np.not_equal)
    __ge__ = _cmp_factory(np.greater_equal)
    __gt__ = _cmp_factory(np.greater)

    def cat(self, other, axis=1):
        x, y = self._x, other._x
        x = x.reshape(x.shape + (1,)*(axis + 1 - x.ndim))
        y = y.reshape(y.shape + (1,)*(axis + 1 - y.ndim))
        return self._wrap(ma.concatenate((x, y), axis=axis))

    def copy(self):
        return self._wrap(self._x.copy())

    def astype(self, dtype):
        return self._wrap(self._x.astype(dtype))

    def masked(self, mask):
        return self._wrap(ma.masked_array(self._x, mask))

    def unmasked(self):
        if not isinstance(self._x, ma.masked_array):
            return self
        else:
            return self._wrap(self._x.data)

    def filled(self, fill=np.nan):
        if not isinstance(self._x, ma.masked_array):
            return self
        else:
            return self._wrap(self._x.filled(fill))


class Signal(SignalBase):
    fmtstr = "%s {name} with shape {shape}"

    def __init__(self, x, t=None, *args, **kw):
        SignalBase.__init__(self, x, t, *args, **kw) 

        self.number = kw.get('number', -1)
        self.name = kw.get('name', "")
        self.label = kw.get('label', self.name)
        self.type = kw.get('type', "")
        self.units = kw.get('units', "")
        self.tunits = kw.get('tunits', "s")

    def shift_t(self, dt):
        return self.__class__(self.x, self.t + dt, *self.args, **self.kw)

    def to_ms(self):
        assert self.tunits == 's', "tunits must be seconds"
        kw = self.kw.copy()
        kw['tunits'] = 'ms'
        return self.__class__(self.x, 1e3*self.t, *self.args, **kw)

    def _op_factory(op):
        def wapply(self, other):
            if isinstance(other, SignalBase):
                other = other.x
            return self.__array_wrap__(op(self.x, other))

        def iapply(self, other):
            if isinstance(other, SignalBase):
                other = other.x
            op(self.x, other, self.x)
            return self
        return wapply, iapply

    __add__, __iadd__ = _op_factory(np.add)
    __sub__, __isub__ = _op_factory(np.subtract)
    __mul__, __imul__ = _op_factory(np.multiply)
    __div__, __idiv__ = _op_factory(np.divide)
    __pow__, __ipow__ = _op_factory(np.power)
    
    def __neg__(self):
        return self.__array_wrap__(-self.x)

    def nonneg(self):
        return self.masked(self < 0)

    def nonpos(self):
        return self.masked(self > 0)

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
        x = ma.masked_invalid(self.x)
        return x.min(), x.max()

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

    def interp_PP(self, t, masked=False):
        t = np.atleast_1d(t)
        return self.__class__(self.PP(t, masked), t, **self.kw)

    def interp(self, t, masked=False):
        t = np.atleast_1d(t)
        i0 = np.searchsorted(self.t, t, 'right') - 1
        im, iM = 0, self.shape[0] - 2
        outl, outr = i0 < im, i0 > iM
        i0[outl], i0[outr] = im, iM
        i1 = i0 + 1
        S = self[np.concatenate((i0, i1)).reshape((2, -1))]
        x0, x1 = S.x
        t0, t1 = S.t
        x = x0 + (x1 - x0) * ((t - t0) / (t1 - t0)).reshape(self.bcast)
        if masked:
            x = ma.masked_array(x, outl | outr)
        return self._result_class(x, t, *self.args, **self.kw)

    @classmethod
    def constfit(cls, x, y, deg=0):
        return y.mean()

    @classmethod
    def linfit(cls, x, y, deg=1):
        xm = x.mean()
        ym = y.mean()
        dx = x - xm
        dy = y - ym
        k = np.dot(dx, dy) / np.dot(dx, dx)
        return np.array((k, ym - k*xm))

    def ppolyfit(self, i0, i1, deg):
        if deg == 0:
            polyfit = self.constfit
        elif deg == 1:
            polyfit = self.linfit
        else:
            polyfit = np.polyfit

        N = len(i1)
        c = np.empty((deg + 1, N))
        S = self.filled()
        x = S.x
        t = S.t
        for i in xrange(N):
            s = slice(i0[i], i1[i])
            ti = t[s]
            c[:,i] = polyfit(ti - ti[0], x[s], deg)
        return c

    def as_PP(self, PP):
        c = self.ppolyfit(PP.i0, PP.i1, PP.c.shape[0] - 1)
        return PP.__class__(c, PP.x, **PP.kw)

    def detrend(self):
        k, d = self.linfit(self.t, self.filled().x)
        return self.__array_wrap__(self.x - k*self.t - d)

    def smooth(self, w=100, mode='gaussian'):
        if mode == 'gaussian':
            x = smooth(self.x, window_len=2*w+1)
        elif mode == 'median':
            x = self.x.copy()
            mediansmooth(x, w)
        else:
            raise ValueError("Unknown smooth mode %s" % mode)
        return self.__array_wrap__(x)

    def despike(self, w=2):
        return self.smooth(w, mode='median')

    def _move_factory(name):
        movefun = movefuns[name]
        def move(self, w=100):
            x = self.x
            x = np.concatenate((2*x[0]-x[w:0:-1], x, 2*x[-1]-x[-2:-w-2:-1]), axis=0)
            return self.__array_wrap__(movefun(x, 2*w+1)[2*w:])
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
        def delta(x):
            return np.concatenate((x[1:2]-x[:1], x[2:]-x[:-2], x[-1:]-x[-2:-1]))

        dx_dt = delta(self.x) / delta(self.t).reshape(self.bcast)
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

    def crossings(self, lvl, threshold=0):
        x0, x1 = self.x[:-1], self.x[1:]
        ind = np.flatnonzero(self.cross(lvl, x0, x1) | self.cross(lvl, x1, x0))
        if threshold == 0:
            ind0, ind1 = ind, ind + 1
        else:
            ind0, ind1 = self.group(ind, threshold)
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
        ylab = self.type
        if len(self.units) > 0:
            ylab += " (%s)" % self.units
        return ylab

    def plot(self, ax=None, **kw):
        ax = get_axes(ax, xlab=self.xlab, ylab=self.ylab)
        kw.setdefault('label', self.label)
        ax.plot(self.t, self.x, **kw)
        return ax

    @memoized_property
    def fs(self):
        return (self.t.size - 1) / (self.t[-1] - self.t[0])

    @memoized_property
    def f_Nyq(self):
        return self.fs / 2

    def filter(self, fm=None, fM=None):
        filter = Filter(fm, fM, self.f_Nyq)
        return self.__array_wrap__(filter(self.x))

    def specgram(self, ax=None, NFFT=2048, step=512, detrend='linear', **kw):
        spec = Specgram(NFFT, step, detrend)(self)
        spec.y *= 1e-3
        spec.yunits = 'kHz'
        return spec.plot(ax, **kw)

    def plot_over_specgram(self, ax=None, yloc=None, *args, **kw):
        ax = self.specgram(ax, *args, **kw)
        ax2 = ax.get_2nd_axes(ylab=self.ylab, yloc=yloc)
        self.plot(ax2)
        ax.figure.tight_layout()
        return ax2

    def xcorr(self, other, w=None, step=None, **kw):
        if other is self:
            S = self
        else:
            S = self.cat(other)
        if w is None:
            w = self.shape[0]
        if step is None:
            step = w
        return XCorr(w=w, step=step, **kw)(S)

    def autocorr(self, **kw):
        return self.xcorr(self, **kw)

    def lag(self, other, sign=1, **kw):
        """
        Returns C, where C.t is the amount of time by which 'self' is lagging 
        behind 'other'
        """
        C = self.xcorr(other, **kw)[0]
        return C[np.argmax(sign*C.x)]


class AmpSignal(Signal):
    def __init__(self, x, t=None, amp=Amp(), **kw):
        Signal.__init__(self, x, t, amp, **kw)
        self._result_class = Signal

        self.amp = amp
        self.dtype = kw.get('dtype', None)

    @property
    def x(self):
        return self.amp(np.asarray(self._x, dtype=self.dtype))

    def apply(self):
        return self._result_class(self.x, self.t, **self.kw)

    def __getitem__(self, indx):
        if not isinstance(indx, tuple):
            indx = (indx,)
        return self.__class__(self._x[indx], self._t[indx[0]], self.amp[indx[0]], **self.kw)
    
    def _op_factory(op, fallback_op):
        def wapply(self, other):
            try:
                return self.__class__(self._x, self._t, op(self.amp, other), **self.kw)
            except TypeError:
                return fallback_op(self, other)

        def iapply(self, other):
            self.amp = op(self.amp, other)
            self.args = (self.amp,)
            return self
        return wapply, iapply

    __add__, __iadd__ = _op_factory(operator.add, Signal.__add__)
    __sub__, __isub__ = _op_factory(operator.sub, Signal.__sub__)
    __mul__, __imul__ = _op_factory(operator.mul, Signal.__mul__)
    __div__, __idiv__ = _op_factory(operator.div, Signal.__div__)


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



