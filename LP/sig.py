import numpy as np
import numpy.ma as ma
import os
import h5py
import copy
import operator

from pdb import set_trace

import scipy.optimize as opt

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sm_pyplot import tight_figure

usetex = False

math_sel = tight_figure.MathSelector(usetex=usetex)

tfigure = tight_figure.pickable_linked_lod_tight_figure
figure = tight_figure.pickable_linked_lod_figure
#tfigure = tight_figure.pickable_linked_tight_figure
#figure = tight_figure.pickable_linked_figure
plot = plt.plot
ion = plt.ion

def get_tfig(*args, **kw):
    kw.setdefault('figure', tfigure)
    return tight_figure.get_fig(*args, **kw)

def get_fig(*args, **kw):
    kw.setdefault('figure', figure)
    return tight_figure.get_fig(*args, **kw)

def get_axes(*args, **kw):
    kw.setdefault('figure', figure)
    return tight_figure.get_axes(*args, **kw)


from mdsclient import *
from mediansmooth import *
from cookb_signalsmooth import smooth

from collections import MutableMapping, Iterable


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


class DictView(MutableMapping):
    def __init__(self, source, valid_keys):
        self.source, self.valid_keys = source, valid_keys

    def __getitem__(self, key):
        if key in self.valid_keys:
            return self.source[key]
        else:
            raise KeyError(key)

    def __len__(self):
        return len(self.valid_keys)

    def __iter__(self):
        for key in self.valid_keys:
            yield key

    def __setitem__(self, key, value):
        if key in self.valid_keys:
            self.source[key] = value
        else:
            raise KeyError(key)

    def __delitem__(self, key):
        self.valid_keys.remove(key)


class PiecewisePolynomial:
    def __init__(self, c, x, **kw):
        self.c, self.x, self.kw = c, x, kw

        if isinstance(self.c, (tuple, list)):
            self.c = np.concatenate([c[None] for c in self.c], axis=0)

        self.fill = kw.get('fill', None)
        self.i0 = kw.get('i0', np.arange(x.size))
        self.i1 = kw.get('i1', None)

        self.N = self.i0.size
        self.shape = self.c.shape[2:]

    def __getitem__(self, index):
        if not isinstance(index, tuple): index = (index,)
        index = (slice(None), slice(None)) + index
        return PiecewisePolynomial(self.c[index], self.x, **self.kw)

    def __call__(self, X, side='right'):
        xi = self.x[self.i0]
        ind = np.searchsorted(xi, X, side) - 1
        outl, outr = ind < 0, ind > self.N-2
        ind[outl], ind[outr] = 0, self.N-2
        
        dX = X - xi[ind]

        #Y = reduce(lambda Y, c: Y*dX + c[ind], self.c)
        
        c = self.c
        Y = c[0,ind].copy()
        for a in c[1:]:
            Y = Y*dX + a[ind]

        if self.fill is not None:
            Y[outl | outr] = self.fill
        return Y

    @memoized_property
    def T(self):
        return PiecewisePolynomial(self.c.swapaxes(2,3), self.x, **self.kw)

    @memoized_property
    def savefields(self):
        return DictView(self.__dict__, ('c', 'x', 'i0'))

    @memoized_property
    def pascal(self):
        n = self.c.shape[0]
        p = np.zeros((n,n), np.int)
        p[:,0] = 1
        for i in xrange(1, n):
            for j in xrange(1, i+1):
                p[i,j] = p[i-1,j-1] + p[i-1,j]
        return p

    def X(self, x):
        n = self.c.shape[0]
        
        def diag_ind(k):
            if k >= 0:
                i = k
                f = n*(n-k)
            else:
                i = (-k) * n
                f = n*n
            return np.arange(i, f, n+1)

        X = np.zeros((n,n))
        y = 1.
        for k in xrange(n):
            X.flat[diag_ind(-k)] = y
            y *= x
        return X

    def calc_c(self, dX, c, out):
        out[::-1] = np.dot(c[::-1], self.pascal*self.X(dX))

    def add_nodes_old(self, x):
        M = self.x.size
        X = np.r_[self.x, x]
        self.x, perm = np.unique(X, return_index=True)
        ind = perm[perm >= M]
        
        n, m = self.c.shape[0], X.size - M
        self.c = np.c_[self.c, np.zeros((n,m+1))]

        ind2 = np.searchsorted(X[:M], X[ind]) - 1

        for i, j in zip(ind, ind2):
            self.calc_c(X[i] - X[j], self.c[:,j], self.c[:,i])
        
        self.c = self.c[:,perm][:,:-1]

        self.N = self.x.size
        self.i0 = np.arange(self.N)

    def add_nodes(self, x):
        M = self.x.size
        X = np.r_[self.x, x]
        self.x, perm = np.unique(X, return_index=True)
        
        ind = np.flatnonzero(perm >= M)
        
        n, m = self.c.shape[0], X.size - M
        nans = np.zeros((n,m+1))
        nans.fill(np.nan)
        self.c = np.c_[self.c, nans][:,perm][:,:-1]

        for i in ind[(0 < ind) & (ind < self.c.shape[1])]:
            dX = self.x[i] - self.x[i-1]
            self.calc_c(dX, self.c[:,i-1], self.c[:,i])

        self.N = self.x.size
        self.i0 = np.arange(self.N)

    def _mask(self, w):
        ind0, ind1 = np.searchsorted(self.i0, w)
        return np.concatenate(map(np.arange, ind0, ind1))

    @staticmethod
    def cat(a, axis=0):
        return a[0].__array_wrap__(ma.concatenate(a, axis))

    def eval(self, w=None):
        try:
            i0 = self.i0[self._mask(w)]
        except:
            i0 = self.i0   

        il, ir = i0[:-1], i0[1:]
        yl, yr = self(self.x[il], 'right'), self(self.x[ir], 'left')

        shape = (il.size + ir.size,)
        i = self.cat((il[:,None], ir[:,None]), 1).reshape(shape)
        y = self.cat((yl[:,None], yr[:,None]), 1).reshape(shape + self.shape)
        return i, y

    def plot(self, ax=None, x=None, w=None):
        if x is None:
            x = self.x

        i, y = self.eval(w=w)
        x = x[i]
                
        ax = get_axes(ax)
        ax.plot(x, y)
        return ax


class IOH5:
    def __init__(self, h5name="test.h5"):
        self.h5name = h5name

    def save(self, d, compression="gzip"):
        f = h5py.File(self.h5name, "w")
        for key, val in d.iteritems():
            f.create_dataset(key, data=val, compression=compression)
        f.close()

    def load(self, d=None):
        f = h5py.File(self.h5name, "r")
        if d is None: 
            d = dict()
        for k in f:
            d[k] = f[k][:]
        f.close()
        return d


class IO:
    def __init__(self):
        pass

    def get_size(self, node):
        pass

    def get_node(self, node):
        pass

    def put_node(self, node, val):
        pass
    
    def load(self, nodes, more_nodes=()):
        if isinstance(nodes, str):
            nodes = (nodes,)
        self.nodes = nodes
        M = len(nodes) + len(more_nodes)
        N = self.get_size(nodes[0])

        dtype = [np.float32]*M
        x = np.empty(N, zip(nodes + more_nodes, dtype))

        for node in nodes:
            x[node][:] = self.get_node(node)

        return x

    def save(self, x, nodes=None):
        nodes = nodes or x.dtype.names
        for node in nodes:
            self.put_node(node, x[node])


class IOFile(IO):
    def __init__(self, shn=0, suffix="", subdir=""):
        self.shn, self.suffix, self.subdir = shn, suffix, subdir

        self.basepath = os.environ['DATAPATH']
        self.fullpath = os.path.join(self.basepath, self.subdir)
        self.fname = str(self.shn) + self.suffix + '.h5'
        
        self.h5name = os.path.join(self.fullpath, self.fname)

        IO.__init__(self)

    def get_size(self, node):
        return self._f[node].len()

    def get_node(self, node):
        return self._f[node].value

    def put_node(self, node, val):
        self._f.create_dataset(node, data=val, compression="gzip")

    def load(self, *args):
        self._f = h5py.File(self.h5name,"r")
        x = IO.load(self, *args)
        self._f.close()
        return x

    def save(self, *args):
        self._f = h5py.File(self.h5name,"w")
        IO.save(self, *args)
        self._f.close()


class TdiError(Exception):
    pass

class IOMds(IO):
    def __init__(self, shn=0, sock=None):
        self.shn = shn
        self.mdsserver = "localhost"
        self.mdsport = "8000"
        self.mdstree = None
        self.mdsfmt = "%s"
        self.datadeco = "data(%s)"
        self.timedeco = "dim_of(%s)"
        self.sizedeco = "size(%s)"

    @memoized_property
    def sock(self):
        s = mdsconnect(self.mdsserver + ':' + str(self.mdsport))
        if self.mdstree is not None:
            mdsopen(s, self.mdstree, self.shn)
        return s

    def _mdsstr(self, node):
        mdsfmt = self.mdsfmt
        if node == 't':
            node = self.nodes[self.nodes.index('t')-1]
            if node == 't':
                raise TdiError("Need other node name to obtain 't'")
            mdsfmt = self.timedeco % mdsfmt
        else:
            mdsfmt = self.datadeco % mdsfmt

        return mdsfmt % node

    def get_size(self, node):
        return self.mdsvalue(self.sizedeco % self._mdsstr(node))

    def get_node(self, node):
        return self.mdsvalue(self._mdsstr(node))

    def save(self, x):
        raise NotImplementedError("Saving to MDS not implemented")

    def mdsvalue(self, *args):
        ret = mdsvalue(self.sock, *args)
        if isinstance(ret, str) and ret.startswith("Tdi"):
            raise TdiError(ret)
        return ret


class Amp:
    def __init__(self, fact=1., offs=0., fixpoints=None):
        if fixpoints is not None:
            (x0, y0), (x1, y1) = fixpoints
            self.fact = (y1-y0)/(x1-x0)
            self.offs = y0 - self.fact*x0
        else:
            self.fact, self.offs = fact, offs

    def __call__(self, x):
        return x*self.fact + self.offs   # try x.__mult__ if x is object

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
        return "Amp instance (fact=%.3f, offs=%.3f)" % (self.fact, self.offs)


class Signal:
    def __init__(self, x, t, **kw):
        self.x, self.t, self.kw = x, t, kw
        
        self.number = kw.get('number', -1)
        self.name = kw.get('name', "")
        self.type = kw.get('type', None)
        self.units = kw.get('units', "")
        self.tunits = kw.get('tunits', "s")

    def __call__(self, t):
        return self.PP(t)

    def __getitem__(self, index):
        return self.__class__(self.x[index], self.t[index], **self.kw)

    def _op_factory(op, calcfun):
        def calcx(self, other):
            if isinstance(other, Signal):
                if other.units != self.units:
                    raise Exception("Unit mismatch")
                other = other.x
            
            return op(self.x, other)

        def calci(self, other):
            x = calcx(self, other)
            return self

        def calc(self, other):
            x = calcx(self, other)
            return self.__class__(x, self.t, **self.kw)
        
        return locals()[calcfun]

    __add__  = _op_factory(operator.add , 'calc')
    __sub__  = _op_factory(operator.sub , 'calc')
    __iadd__ = _op_factory(operator.iadd, 'calci')
    __isub__ = _op_factory(operator.isub, 'calci')
    __lt__   = _op_factory(operator.lt  , 'calcx')
    __le__   = _op_factory(operator.le  , 'calcx')
    __eq__   = _op_factory(operator.eq  , 'calcx')
    __ne__   = _op_factory(operator.ne  , 'calcx')
    __ge__   = _op_factory(operator.ge  , 'calcx')
    __gt__   = _op_factory(operator.gt  , 'calcx')

    def __mul__(self, other):
        return self.__class__(other*self.x, self.t, **self.kw)

    def __imul__(self, other):
        try:
            other.apply(self.x)
        except AttributeError:
            self.x[:] *= other
        return self

    @property
    def size(self):
        return self.x.size

    def copy(self):
        s = copy.copy(self)
        s.x = s.x.copy()
        return s

    def trim(self, s):
        self.x, self.t = self.x[s], self.t[s]
        return self

    def range(self):
        return np.nanmin(self.x), np.nanmax(self.x)

    def plot_range(self):
        x = np.sort(self.x)
        i = np.round(np.array([0.001, 0.999])*(x.size-1)).astype('i')
        return tuple(x[i])

    def norm_to_region(self, cnd):
        self.x[:] -= self.x[cnd].mean()
        return self

    @memoized_property
    def PP(self):
        x, t = self.x.astype('d'), self.t.astype('d')
        x0, x1 = x[:-1], x[1:]
        t0, t1 = t[:-1], t[1:]
        dx_dt = (x1 - x0)/(t1 - t0)
        return PiecewisePolynomial((dx_dt, x0), t)

    def smooth(self, w=100):
        self.x[:] = smooth(self.x, window_len=2*w+1)
        return self

    def mediansmooth(self, w=100):
        mediansmooth(self.x, w)
        return self

    def despike(self, w=2):
        return self.mediansmooth(w)

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

    def plot(self, ax=None):
        ax = get_axes(ax)
        ax.plot(self.t, self.x.T)
        return ax


class PeriodPhaseFinder:
    def __init__(self, x):
        self.x = x

    def find_f(self):
        nextpow2 = lambda x: 2**np.ceil(np.log2(x)).astype('i')
        Nfft = nextpow2(self.x.size)

        X = np.fft.fft(self.x-self.x.mean(), Nfft)
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

    def plot_plunges(self, ax=None):
        ax = get_axes(ax)
        Signal.plot(self, ax)
        i0, iM, i1 = self.t_ind
        ax.plot(self.t[i0], self.x[i0], 'r*')
        ax.plot(self.t[i1], self.x[i1], 'g*')
        ax.plot(self.t[iM], self.x[iM], 'm*')
        return ax


class VoltageSignal(Signal):
    def __init__(self, x, t, **kw):
        kw.setdefault('type', 'Voltage')
        kw.setdefault('units', 'V')

        self.min_ptp = kw.pop('min_ptp', 20)

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
        return self.x.ptp() > D and self.PPF.D > D

    def plot_sweeps(self, ax=None):
        ax = get_axes(ax)
        Signal.plot(self, ax)
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
        s.V = s.V[index]
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


class Digitizer:
    def __init__(self, shn=0, sock=None, name=""):
        self.shn, self.sock, self.name = shn, sock, name

        self.IO_mds = self.IO_file = None
        self.nodes = self.more_nodes = ()
        self.window = slice(None)
        self.amp = dict()

    @memoized_property
    def x(self):
        return self.load()

    def load_raw_mds(self):
        self.x = self.IO_mds.load(self.nodes, self.more_nodes)
        return self.x
        
    def load_raw_file(self):
        self.x = self.IO_file.load(self.nodes, self.more_nodes)
        return self.x

    def load_raw(self):
        try:
            self.load_raw_file()
        except IOError:
            self.load_raw_mds()
            self.save()
        return self.x
    
    def save(self):
        self.IO_file.save(self.x, self.nodes)

    def calib(self):
        for node in self.nodes:
            try:
                self.amp[node].apply(self.x[node])
            except KeyError:
                pass

    def _load_calib(self, loadfun):
        loadfun()
        self.calib()
        return self.x

    def load_mds(self):
        return self._load_calib(self.load_raw_mds)

    def load_file(self):
        return self._load_calib(self.load_raw_file)

    def load(self):
        return self._load_calib(self.load_raw)

    def calib_offset(self):
        self.load_raw()
        offs = [np.median(self.x[node]) for node in self.nodes]
        return offs

    def plot(self, fig=None):
        nodes = list(self.nodes)
        nodes.remove('t')
        n = len(nodes)

        fig = get_fig(fig, (n, 1), xlab='t (s)')
        
        t = self.x['t']
        for node, ax in zip(nodes, fig.axes):
            ax.plot(t, self.x[node])
            ax.set_ylabel(node)
        fig.canvas.draw()
        return fig



