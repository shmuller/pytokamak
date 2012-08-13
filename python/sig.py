import numpy as np
import os
import h5py
import copy

import scipy.optimize as opt

import matplotlib.pyplot as plt
import tight_figure
reload(tight_figure)

from mdsclient import *

figure = tight_figure.pickable_linked_figure
plot = plt.plot

from collections import Mapping

class DictView(Mapping):
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


class IO:
    def __init__(self):
        pass

    def get_size(self, node):
        pass

    def get_node(self, node):
        pass

    def put_node(self, node, val):
        pass
    
    def load(self, nodes):
        if isinstance(nodes, str):
            nodes = (nodes,)
        self.nodes = nodes
        M = len(nodes)
        N = self.get_size(nodes[0])

        dtype = [np.float32]*M
        dtype[nodes.index('t')] = np.float32
        x = np.empty(N, zip(nodes, dtype))

        for node in nodes:
            x[node][:] = self.get_node(node)

        return x

    def save(self, x):
        for node in x.dtype.names:
            self.put_node(node, x[node])


class IOFile(IO):
    def __init__(self, shn=0, suffix="", subdir=""):
        self.shn, self.suffix, self.subdir = shn, suffix, subdir

        pth = os.environ['DATAPATH']
        pth = os.path.join(pth, self.subdir)
        fname = str(self.shn) + self.suffix + '.h5'
        
        self.h5name = os.path.join(pth, fname)

        IO.__init__(self)

    def get_size(self, node):
        return self._f[node].len()

    def get_node(self, node):
        return self._f[node].value

    def put_node(self, node, val):
        self._f.create_dataset(node, data=val, compression="gzip")

    def load(self, nodes):
        self._f = h5py.File(self.h5name,"r")
    
        x = IO.load(self, nodes)

        self._f.close()
        return x

    def save(self, x):
        self._f = h5py.File(self.h5name,"w")

        IO.save(self, x)

        self._f.close()


class TdiError(Exception):
    pass

class IOMds(IO):
    def __init__(self, shn=0, sock=None):
        self.shn, self._sock = shn, sock
        self.mdsport, self.mdsfmt = "8000", ""
        self.datadeco = "data(%s)"
        self.timedeco = "dim_of(%s)"
        self.sizedeco = "size(%s)"

    def get_sock(self):
        if self._sock is None:
            self.set_sock()
        return self._sock

    def set_sock(self):
        self._sock = mdsconnect('localhost:' + str(self.mdsport))

    sock = property(get_sock, set_sock)

    def _mdsstr(self, node):
        mdsfmt = self.mdsfmt
        if node == 't':
            node = self.nodes[self.nodes.index('t')-1]
            if node == 't':
                raise TdiError("Need other node name to obtain 't'")
            mdsfmt = self.timedeco % mdsfmt
        else:
            mdsfmt = self.datadeco % mdsfmt

        return mdsfmt % (self.shn, node)

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
        return self.fact*x + self.offs

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


class Signal:
    def __init__(self, x, t, name="", type=None):
        self.x, self.t, self.name, self.type = x, t, name, type

    def __getitem__(self, index):
        return Signal(self.x[index], self.t[index], self.name, self.type)

    def __add__(self, other):
        return Signal(other+self.x, self.t, self.name, self.type)

    def __mul__(self, other):
        return Signal(other*self.x, self.t, self.name, self.type)

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

    def range(self):
        return self.x.min(), self.x.max()

    def plot_range(self):
        x = np.sort(self.x)
        i = np.round(np.array([0.001, 0.999])*(x.size-1)).astype('i')
        return tuple(x[i])

    def norm_to_region(self, cnd):
        self.x[:] -= self.x[cnd].mean()

    def deriv(self):
        delta = lambda x: np.r_[x[1]-x[0], x[2:]-x[:-2], x[-1]-x[-2]]
        dx_dt = delta(self.x)/delta(self.t)
        return dx_dt

    def crossings(self, lvl, threshold):
        def cross(lvl, x0, x1):
            return (x0 < lvl) & (lvl < x1)
        
        def group(ind):
            di = np.diff(np.r_[0,ind])
            j = np.flatnonzero(di > threshold*np.median(di))
            ind0, ind1 = ind[j], ind[np.r_[j[1:]-1, -1]] + 1
            return ind0, ind1

        x0, x1 = self.x[:-1], self.x[1:]
        cnd = cross(lvl, x0, x1) | cross(lvl, x1, x0)
        ind0, ind1 = group(np.flatnonzero(cnd))

        is_rising = self.x[ind0] < self.x[ind1]
        return ind0, ind1, is_rising

    def apply_fun(self, fun, ind0=0, ind1=None):
        slices = map(slice, ind0, ind1)
        res = [fun(self.x[s]) for s in slices]
        return np.array(res)
    
    def apply_argfun(self, argfun, ind0=0, ind1=None):
        return ind0 + self.apply_fun(argfun, ind0, ind1)
    
    def local_argmin(self, *args):
        return self.apply_argfun(np.argmin, *args)

    def local_argmax(self, *args):
        return self.apply_argfun(np.argmax, *args)

    def plot(self, newfig=True):
        if newfig: figure()
        plot(self.t, self.x.T)


class PeriodPhaseFinder:
    def __init__(self, x):
        self.x = x
        self.di = np.zeros(2)
        self.iE = None

    def find_f(self):
        nextpow2 = lambda x: 2**np.ceil(np.log2(x)).astype('i')
        Nfft = nextpow2(self.x.size)

        X = np.fft.fft(self.x-self.x.mean(), Nfft)
        iM = np.abs(X[1:Nfft/2+1]).argmax()+1
        f = np.float(iM)/Nfft
        return f

    def cumdist(self, p):
        d, i = p[0], p[1]
        M = self.x.size
        iE = np.round(np.arange(i,M-1,d/2)).astype('i')
        
        dx = np.diff(self.x[iE])
        D = -np.sqrt(dx.dot(dx)/dx.size)

        self.iE = iE
        return D

    def guess_d(self, p0):
        d0, i0 = p0[0], p0[1]
        N = 200
        dd = np.linspace(0.99*d0, 1.01*d0, N)
        DD = np.empty(N)
        for j in xrange(N):
            DD[j] = self.cumdist([dd[j],i0])
        d = dd[DD.argmin()]
        return d

    def chop_sweeps(self):
        f0 = self.find_f()
        d0 = 1./f0
        x0 = self.x[:np.round(d0)]
        i0 = min(x0.argmin(), x0.argmax())
        
        d = self.guess_d([d0,i0])

        p = opt.fmin(self.cumdist, np.array([d,i0]))

        D = self.cumdist([d0,i0])
        print D
        
        D = self.cumdist([d,i0])
        print D

        D = self.cumdist(p)
        print D

        return D

    def get_iE(self):
        return self.iE


class PositionSignal(Signal):
    def __init__(self, x, t, **kw):
        _kw = dict(
                name="",
                baseline_slice=slice(None, 1000), 
                lvl_fact=20,
                dist_threshold=1000)
        _kw.update(kw)
        self.__dict__.update(_kw)
        self.kw = DictView(self.__dict__, _kw.keys())
        self._t_ind = None

        Signal.__init__(self, x, t, self.name, 'Position')
        
    def __getitem__(self, index):
        return PositionSignal(self.x[index], self.t[index], **self.kw)

    def get_baseline(self):
        return self.x[self.baseline_slice]

    def get_crossings(self):
        x = self.get_baseline()
        xm, xs = x.mean(), x.std()
        lvl = xm + self.lvl_fact*xs
        return self.crossings(lvl, self.dist_threshold)

    def get_t_ind(self):
        if self._t_ind is None:
            self.set_t_ind()
        return self._t_ind

    def set_t_ind(self):
        ind0, ind1, is_rising = self.get_crossings()
        i0, i1 = ind0[is_rising], ind1[~is_rising]
        iM = self.local_argmax(i0, i1)
        self._t_ind = i0, iM, i1

    t_ind = property(get_t_ind, set_t_ind)

    def regions(self, fun=None):
        i0, iM, i1 = self.t_ind
        return map(fun, i0, i1)

    def get_slices(self):
        return self.regions(fun=slice)

    def get_mask(self):
        return np.concatenate(self.regions(fun=np.arange))

    def plot_plunges(self):
        self.plot(self)
        i0, iM, i1 = self.t_ind
        plot(self.t[i0], self.x[i0], 'r*')
        plot(self.t[i1], self.x[i1], 'g*')
        plot(self.t[iM], self.x[iM], 'm*')


class VoltageSignal(Signal):
    def __init__(self, x, t, name=""):
        Signal.__init__(self, x, t, name, 'Voltage')
        self.iE = None

    def __getitem__(self, index):
        return VoltageSignal(self.x[index], self.t[index], self.name)

    def chop_sweeps(self):
        PPF = PeriodPhaseFinder(self.x)
        PPF.chop_sweeps()
        self.iE = PPF.get_iE()

    def plot(self, newfig=True):
        Signal.plot(self, newfig)
        if self.iE is not None:
            plot(self.t[self.iE], self.x[self.iE], 'r+')


class CurrentSignal(Signal):
    def __init__(self, x, t, V=None, name="", C=None):
        Signal.__init__(self, x, t, name, 'Current')
        self.V, self.C = V, C

    def __getitem__(self, index):
        return CurrentSignal(self.x[index], self.t[index], self.name)

    def capa_pickup(self):
        cnd = self.t - self.t[0] < 0.05
        dV_dtc = self.V.deriv()[cnd]
        N = (dV_dtc*dV_dtc).sum()

        dI = self.copy()
        dI.norm_to_region(cnd)
        self.C = (dI.x[cnd]*dV_dtc).sum()/N

    def I_capa(self):
        if self.C is None:
            self.capa_pickup()
        return self.C*self.V.deriv()

    def I_corr(self):
        I_capa = self.I_capa()
        s = self.copy()
        s.x[:] -= I_capa
        return s


class Digitizer:
    def __init__(self, shn=0, sock=None, name=""):
        self.shn, self.sock, self.name = shn, sock, name

        self.IO_mds = self.IO_file = None
        self.nodes = ()
        self.window = slice(None)
        self.amp = None

    def load_raw_mds(self):
        self.x = self.IO_mds.load(self.nodes)
        return self.x
        
    def load_raw_file(self):
        self.x = self.IO_file.load(self.nodes)
        return self.x

    def load_raw(self):
        try:
            self.x = self.IO_file.load(self.nodes)
        except:
            self.x = self.IO_mds.load(self.nodes)
            self.save()
        return self.x
    
    def save(self):
        self.IO_file.save(self.x)

    def calib(self):
        for node in self.nodes:
            self.amp[node].apply(self.x[node])

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

    def plot(self):
        fig = figure()
        nodes = list(self.nodes)
        nodes.remove('t')
        n = len(nodes)
        t = self.x['t']
        for i, node in enumerate(nodes):
            ax = fig.add_subplot(n, 1, 1+i)
            ax.grid(True)
            ax.plot(t, self.x[node])



