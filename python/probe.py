import numpy as np
import os
import h5py
import copy

import matplotlib.pyplot as plt

import scipy.interpolate as interp
import scipy.optimize as opt
import scipy.odr as odr

from mdsclient import *

from ipdb import set_trace

import tight_figure
reload(tight_figure)

figure = tight_figure.pickable_linked_figure
plot = plt.plot
ion = plt.ion
draw = plt.draw
hold = plt.hold
xlim = plt.xlim
ylim = plt.ylim
xlabel = plt.xlabel
ylabel = plt.ylabel
grid = plt.grid

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


class MouseMotionObserverViewer:
    def __init__(self, observers, viewer, plotfun):
        self.observers, self.viewer = observers, viewer
        self.plotfun = plotfun

        self.observer_canvas = self.observers[0].figure.canvas
        self.viewer_canvas = self.viewer.figure.canvas
        
        self.cid = self.observer_canvas.mpl_connect('motion_notify_event', self.on_move)
        self.viewer_canvas.mpl_connect('resize_event', self.on_resize)
        self.viewer_canvas.mpl_connect('close_event', self.on_close)

        self.on_resize(None)

    def set_viewer_visible(self, value):
        for line in self.viewer.lines:
            line.set_visible(value)

    def on_resize(self, event):
        self.set_viewer_visible(False)
        self.viewer_canvas.draw()
        self.background = self.viewer_canvas.copy_from_bbox(self.viewer.bbox)
        self.set_viewer_visible(True)

    def on_move(self, event):
        if event.inaxes in self.observers:
            self.viewer_canvas.restore_region(self.background)
                
            self.plotfun(event)

            for line in self.viewer.lines:
                self.viewer.draw_artist(line)
            self.viewer_canvas.blit(self.viewer.bbox)

    def on_close(self, event):
        self.observer_canvas.mpl_disconnect(self.cid)


class Data(dict):
    def __init__(self, nodes, X):
        self.nodes, self.X = nodes, X
        dict.__init__(self, zip(nodes, X))

    def view(self, s):
        return dict(zip(self.nodes, self.X[:,s]))


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

        X = np.empty((M,N),'f')

        for i in xrange(M):
            X[i,:] = self.get_node(nodes[i])

        return Data(nodes, X)

    def save(self, x):
        for node, val in x.iteritems():
            self.put_node(node, val)


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
        self._f.create_dataset(node,data=val,compression="gzip")

    def load(self, nodes):
        self._f = h5py.File(self.h5name,"r")
    
        x = IO.load(self, nodes)

        self._f.close()
        return x

    def save(self, x):
        self._f = h5py.File(self.h5name,"w")

        IO.save(self, x)

        self._f.close()


class IOMds(IO):
    def __init__(self, shn=0, sock=None):
        self.shn, self.sock = shn, sock
        self.mdsport, self.mdsfmt = "8000", ""

    def _mdsstr(self, node):
        mdsfmt = self.mdsfmt
        if node == 't':
            node = self.nodes[self.nodes.index('t')-1]
            if node == 't':
                raise RuntimeError("Need other node name to obtain 't'")
            mdsfmt = 'dim_of(%s)' % mdsfmt

        return mdsfmt % (self.shn, node)

    def get_size(self, node):
        return mdsvalue(self.sock,'size(%s)' % self._mdsstr(node))

    def get_node(self, node):
        return mdsvalue(self.sock, self._mdsstr(node))

    def load(self, nodes):
        if self.sock is None:
            self.sock = mdsconnect('localhost:' + str(self.mdsport))
        x = IO.load(self, nodes)
        return x

    def save(self, x):
        raise NotImplementedError("Saving to MDS not implemented") 


class PiecewiseLinear:
    def __init__(self, y, x):
        self.y, self.x = y, x
        self.N = self.x.shape[0]
        self.shape = self.y.shape[2:]

    def __getitem__(self, index):
        if not isinstance(index, tuple): index = (index,)
        index = (slice(None), slice(None)) + index
        return PiecewiseLinear(self.y[index], self.x)

    def plot(self):
        shape = (self.N*2,)
        x = self.x.reshape(shape)
        y = self.y.reshape(shape + self.shape)
        plot(x, y)

    def plot2(self):
        shape = (self.N, 1)
        nanx = np.empty(shape) + np.nan
        nany = np.empty(shape + self.shape) + np.nan

        shape = (self.N*3,)
        x = np.concatenate((self.x, nanx), 1).reshape(shape)
        y = np.concatenate((self.y, nany), 1).reshape(shape + self.shape)
        plot(x, y)


class PiecewisePolynomial:
    def __init__(self, c, x, **kw):
        self.kw = {'fill': None, 'i0': np.arange(x.size), 'i1': None}
        self.kw.update(kw)
        for key, val in self.kw.iteritems():
            setattr(self, key, val)

        self.c, self.x = c, x
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

        Y = self.c[0,ind].copy()
        for a in self.c[1:]:
            Y = Y*dX + a[ind]

        if self.fill is not None:
            Y[outl | outr] = self.fill
        return Y

    @property
    def T(self):
        return PiecewisePolynomial(self.c.swapaxes(2,3), self.x, **self.kw)

    def plot(self, newfig=True, x=None):
        xi = self.x[self.i0]
        xl, xr = xi[:-1], xi[1:]
        yl, yr = self(xl, 'right'), self(xr, 'left')

        if x is not None:
            xi = x[self.i0]
            xl, xr = xi[:-1], xi[1:]

        shape = ((self.N-1)*2,)
        x = np.concatenate((xl[:,None], xr[:,None]), 1).reshape(shape)
        y = np.concatenate((yl[:,None], yr[:,None]), 1).reshape(shape + self.shape)
        
        if newfig: figure()
        return plot(x, y)


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


class Signal:
    def __init__(self, t, x, name="", type=None):
        self.t, self.x, self.name, self.type = t, x, name, type

    def __getitem__(self, index):
        return Signal(self.t[index], self.x[index], self.name, self.type)

    @property
    def size(self):
        return self.x.size

    def copy(self):
        s = copy.copy(self)
        s.x = s.x.copy()
        return s

    def trim(self, s):
        self.t, self.x = self.t[s], self.x[s]

    def range(self):
        return self.x.min(), self.x.max()

    def plot_range(self):
        x = np.sort(self.x)
        i = np.round(np.array([0.001, 0.999])*(x.size-1)).astype('i')
        return tuple(x[i])

    def norm_to_region(self, cnd):
        s = self.copy()
        s.x[:] -= self.x[cnd].mean()
        return s

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


class PositionSignal(Signal):
    def __init__(self, t, x, **kw):
        _kw = dict(
                name="",
                baseline_slice=slice(None, 1000), 
                lvl_fact=20,
                dist_threshold=1000)
        _kw.update(kw)
        self.__dict__.update(_kw)
        self.kw = DictView(self.__dict__, _kw.keys())
        self._t_ind = None

        Signal.__init__(self, t, x, self.name, 'Position')
        
    def __getitem__(self, index):
        return PositionSignal(self.t[index], self.x[index], **self.kw)

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
    def __init__(self, t, x, name=""):
        Signal.__init__(self, t, x, name, 'Voltage')
        self.iE = None

    def __getitem__(self, index):
        return VoltageSignal(self.t[index], self.x[index], self.name)

    def chop_sweeps(self):
        PPF = PeriodPhaseFinder(self.x)
        PPF.chop_sweeps()
        self.iE = PPF.get_iE()

    def plot(self, newfig=True):
        Signal.plot(self, newfig)
        if self.iE is not None:
            plot(self.t[self.iE], self.x[self.iE], 'r+')


class CurrentSignal(Signal):
    def __init__(self, t, x, V=None, name="", C=None):
        Signal.__init__(self, t, x, name, 'Current')
        self.V, self.C = V, C

    def __getitem__(self, index):
        return CurrentSignal(self.t[index], self.x[index], self.name)

    def capa_pickup(self):
        cnd = self.t - self.t[0] < 0.05
        dV_dtc = self.V.deriv()[cnd]
        N = (dV_dtc*dV_dtc).sum()

        dI = self.norm_to_region(cnd)
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


def wrap_fmin(fun, p0, x, y):
    def dy2(p):
        dy = fun(p, x) - y
        return dy.dot(dy)/dy.size

    return opt.fmin(dy2, p0, disp=False)

def wrap_odr(fun, p0, x, y):
    mod = odr.Model(fun)
    dat = odr.Data(x, y)
    o = odr.ODR(dat, mod, p0)
    out = o.run()
    return out.beta

engines = {'fmin': wrap_fmin, 'odr': wrap_odr}


class Fitter:
    def __init__(self, x, y, engine='fmin'):
        self.x, self.y = x, y

        self.OK = None
        self.X = self.Y = None
        self.P = self.P0 = None
        self.p = self.p0 = None

        self.set_engine(engine)

    # overload
    def set_OK(self):
        self.OK = True
        return self.OK

    def set_norm(self):
        self.X, self.Y = self.x, self.y
        return self.X, self.Y

    def set_unnorm(self):
        self.p, self.p0 = self.P, self.P0
        return self.p, self.p0

    def set_guess(self):
        self.P0 = 0.
        return self.P0

    def fitfun(self, P, X):
        pass

    # static
    def is_OK(self):
        if self.OK is None:
            self.set_OK()
        return self.OK

    def get_norm(self):
        if self.X is None:
            self.set_norm()
        return self.X, self.Y

    def get_unnorm(self):
        if self.p is None:
            self.set_unnorm()
        return self.p, self.p0

    def get_guess(self):
        if self.P0 is None:
            self.set_guess()
        return self.P0
    
    def set_engine(self, engine):
        self.engine = engines[engine]

    def fit(self, P0=None):
        if not self.is_OK():
            raise RuntimeError("Cannot fit data that failed is_OK() check")
        
        if P0 is None:
            P0 = self.get_guess()
        self.P = self.engine(self.fitfun, P0, self.X, self.Y)
        self.set_unnorm()

        return self.p

    def eval_guess_norm(self, X):
        return self.fitfun(self.P0, X)

    def eval_guess(self, x):
        return self.fitfun(self.p0, x)

    def eval_norm(self, X):
        return self.fitfun(self.P, X)

    def eval(self, x):
        return self.fitfun(self.p, x)

    def get_XY(self):
        if self.is_OK(): 
            if self.P is None:
                self.fit()
            return self.X, (self.Y, self.eval_norm(self.X))
        else:
            return self.X, (self.Y,)

    def get_xy(self):
        if self.is_OK(): 
            if self.p is None:
                self.fit()
            return self.x, (self.y, self.eval(self.x))
        else:
            return self.x, (self.y,)

    def plot(self, ax=None, fun='get_xy', lines=None):
        x, y = getattr(self, fun)()
        sty = (dict(linewidth=1.5),)
        
        if ax is None:
            ax = figure().gca()

        if lines is None: lines = []
        Nl, Ny = len(lines), len(y)

        for li, yi in zip(lines, y):
            li.set_data(x, yi)
            li.set_visible(True)

        for li in lines[Ny:]:
            li.set_visible(False)

        if Nl > 0:
            color = lines[0].get_color()
        else:
            color = ax._get_lines.color_cycle.next()
            sty = (dict(),) + sty

        for yi, styi in zip(y[Nl:], sty):
            li, = ax.plot(x, yi, color=color, **styi)
            lines.append(li)

        return lines


def LP_unnormalize(P, Vm, VM, Im, IM):
    dV, dI = VM-Vm, IM-Im
    a = 2./dI
    b = -(Im*a+1)

    p = np.empty_like(P)
    p[0] = (P[0]-b)/a
    p[1] = (P[1]+P[2]*np.log(1.-b/P[0]))*dV+Vm
    p[2] = P[2]*dV
    return p

def LP_normalize(p, Vm, VM, Im, IM):
    dV, dI = VM-Vm, IM-Im
    a = 2./dI
    b = -(Im*a+1)

    P = np.empty_like(p)
    P[0] = a*p[0]+b
    P[1] = (p[1]-p[2]*np.log(1.-b/P[0])-Vm)/dV
    P[2] = p[2]/dV
    return P


class FitterIV(Fitter):
    def __init__(self, V, I, mask=None, cut_at_min=True, **kw):
        self.mask_ind = np.arange(V.size)
        if mask is not None:
            self.mask_ind = self.mask_ind[mask]

        self.sort_ind = V[self.mask_ind].argsort()
        self.ind = self.mask_ind[self.sort_ind]

        self.V, self.I = V[self.ind], I[self.ind]

        self.im = self.I.argmin()
        self.Vm, self.VM = self.V[0], self.V[-1]
        self.Im, self.IM = self.I[self.im], np.median(self.I[:self.I.size/2])
        self.dV = self.VM - self.Vm
        self.dI = self.IM - self.Im

        self.cut_at_min = cut_at_min
        self.M = self.V.size

        Fitter.__init__(self, self.V, self.I, **kw)

    def get_ind(self):
        return self.ind[:self.M]

    def set_OK(self):
        def medianstd(x):
            xm = np.median(x)
            dx = x-xm
            xs = np.sqrt(dx.dot(dx)/dx.size)
            return xm, xs

        N = self.I.size
        Ilm, Ils = medianstd(self.I[:N/10])
        Irm, Irs = medianstd(self.I[N*9/10:])

        cnd1 = Ilm > Ils
        cnd2 = Irm < -2*Ils
        self.OK = cnd1 & cnd2
        return self.OK

    def set_norm(self):
        if self.cut_at_min:
            self.M = self.im+1
        self.X = (self.V[:self.M].astype('d') - self.Vm)/self.dV
        self.Y = (self.I[:self.M].astype('d') - self.Im)/self.dI*2 - 1
        return self.X, self.Y

    def set_unnorm(self):
        self.p0 = LP_unnormalize(self.P0, self.Vm, self.VM, self.Im, self.IM)
        self.p  = LP_unnormalize(self.P , self.Vm, self.VM, self.Im, self.IM)
        return self.p, self.p0

    def set_guess(self):
        def find_i0(x):
            return np.flatnonzero(np.diff(np.sign(x)))[-1]

        V, I = self.get_norm()
        I0 = 1.
        i0 = find_i0(I)
        Vf = (V[i0]+V[i0+1])/2
        Te = (V[self.im]-Vf)/np.log(1-I[self.im]/I0)

        self.P0 = np.array([I0, Vf, Te])
        return self.P0

    def fitfun(self, P, X):
        return P[0]*(1.-np.exp((X-P[1])/P[2]))

    def fit(self):
        Fitter.fit(self)
        save = self.X, self.Y
        Y0 = 0.
        while True:
            self.M *= 0.95
            self.X, self.Y = self.X[:self.M], self.Y[:self.M]
            P_old = self.P
            Fitter.fit(self, P0=self.P)
            
            if np.any(self.P > P_old) or (self.eval_norm(self.X[-1]) > Y0):
                self.P = P_old
                break

        self.X, self.Y = save
        self.set_unnorm()
        return self.p


class IVChar:
    def __init__(self, V, I, mask=None, **kw):
        self.V, self.I = V, I
        self.fitter_IV = FitterIV(self.V.x, self.I.x, **kw)

    def fit(self, out=None):
        if out is None:
            out = np.empty(3)
        try:
            out[:] = self.fitter_IV.fit()
        except RuntimeError:
            out[:] = np.nan
        return out
    
    def mask(self, out=None):
        if out is None:
            out = np.empty(self.I.size, bool)
        out.fill(False)
        if self.fitter_IV.p is not None:
            out[self.fitter_IV.get_ind()] = True
        return out

    def plot(self, ax=None, fun='get_xy', lines=None):
        return self.fitter_IV.plot(ax, fun, lines)


class IVGroup:
    def __init__(self, V, II, s=slice(None), **kw):
        self.IV_char = np.empty(len(II), object)
        for j, I in enumerate(II):
            self.IV_char[j] = IVChar(V[s], I[s], **kw)

    def __getitem__(self, index):
        return self.IV_char[index]

    def fit(self, out=None):
        if out is None:
            out = np.empty((self.IV_char.size, 3))
        for p, IV_char in zip(out, self.IV_char):
            IV_char.fit(p)
        return out

    def mask(self, out=None):
        if out is None:
            out = np.empty((self.IV_char.size, self.IV_char[0].I.size), bool)
        for I, IV_char in zip(out, self.IV_char):
            IV_char.mask(I)
        return out

    def plot(self, ax=None, fun='get_xy'):
        if ax is None:
            ax = figure().gca

        cache = getattr(ax, 'lines_cache', None) 

        if cache is None:
            ax.lines_cache = [x.plot(ax, fun) for x in self.IV_char]
        else:
            ax.lines_cache = [x.plot(ax, fun, l) for x, l in zip(self.IV_char, cache)]


class IVSeries:
    def __init__(self, V, II, iE, **kw):       
        self.V, self.II, self.iE = V, II, iE
        N = iE.shape[0]
        self.siz = (N, len(II))
        self.ti = V.t[iE]

        self.V_range = V.plot_range()
        self.I_range = self._plot_range(II)

        self.IV_group = np.empty(N, object)
        for j in xrange(N):
            s = self._slice(j)
            self.IV_group[j] = IVGroup(V, II, s, **kw)

    def __getitem__(self, index):
        return self.IV_group[index]

    def _plot_range(self, II):
        I_range = np.array([I.plot_range() for I in II])
        return I_range[:,0].min(), I_range[:,1].max()

    def _slice(self, j):
        return slice(self.iE[j,0], self.iE[j,1]+1)

    def fit(self):
        out = np.empty(self.siz + (3,))
        for p, IV_group in zip(out, self.IV_group):
            IV_group.fit(p)

        i0 = np.r_[self.iE[:,0], self.iE[-1,1]]
        self.PP = PiecewisePolynomial(out[None], self.V.t, i0=i0)
        return self.PP

    def mask(self):
        out = np.zeros((len(self.II), self.II[0].size), bool)
        for j, IV_group in enumerate(self.IV_group):
            s = self._slice(j)
            IV_group.mask(out[:,s])
        return out

    def eval(self):
        V = self.V
        PP = self.PP(V.t)
        fitfun = self[0][0].fitter_IV.fitfun
        out = fitfun(PP.T, V.x)
        out[~self.mask()] = np.nan
        return out

    def plot(self):
        fig = figure()
        n = len(self.II)
        IIfit = self.eval()
        for i, I, Ifit in zip(xrange(n), self.II, IIfit):
            ax = fig.add_subplot(n, 1, 1+i)
            plot(I.t, I.x, I.t, Ifit)
            grid(True)
            ylabel(I.name)
        xlabel("t [s]")
        
    def animate(self, fun='get_xy'):
        ion()
        ax = figure().gca()
        if fun == 'get_xy':
            ax.set_xlim(self.V_range)
            ax.set_ylim(self.I_range)
        else:
            ax.set_xlim(( 0.0, 1.0))
            ax.set_ylim((-1.2, 1.2))

        for IV_group in self.IV_group:
            IV_group.plot(ax, fun)
            draw()



class FitterIV2(Fitter):
    def __init__(self, V, I, mask, p0, **kw):
        self.V, self.I = V[mask], I[mask]
        self.Vm, self.VM = self.V.min(), self.V.max()
        self.Im, self.IM = self.I.min(), self.I.max()
        self.dV = self.VM - self.Vm
        self.dI = self.IM - self.Im

    def set_norm(self):
        self.X = (self.V.astype('d') - self.Vm)/self.dV
        self.Y = (self.I.astype('d') - self.Im)/self.dI*2 - 1
        return self.X, self.Y

    def set_unnorm(self):
        self.p = LP_unnormalize(self.P , self.Vm, self.VM, self.Im, self.IM)
        return self.p


class IVSeries2:
    def __init__(self, V, II, PP, mask, iE=None, **kw):
        self.V, self.II, self.PP, self.mask, self.iE = V, II, PP, mask, iE

    def plot(self, i, j):
        s = slice(self.iE[i], self.iE[j]+1)
        V = self.V[s]
        I = self.II[0][s]
        pp = self.PP[s, 0]
        m = self.mask[0, s]

        t = V.t[m]
        V = V.x[m]
        I = I.x[m]
        pm = pp[m]

        print pp.shape

        def fitfun(p, V): 
            return p[0]*(1.-np.exp((V-p[1])/p[2]))

        a = (t-t[0])/(t[-1]-t[0])

        def fitfun2(p, V):
            Is = p[0] + a*(p[3]-p[0])
            Vf = p[1] + a*(p[4]-p[1])
            Te = p[2] + a*(p[5]-p[2])
            return Is*(1.-np.exp((V-Vf)/Te))

        Ifit = fitfun(pm.T, V)

        """
        p0 = pm.mean(0)
        p = wrap_fmin(fitfun, p0, V, I)
        Ifit2 = fitfun(p, V)
        """

        p0 = np.r_[pm[0], pm[-1]]
        p = wrap_fmin(fitfun2, p0, V, I)
        Ifit2 = fitfun2(p, V)

        figure()
        plot(t, I, t, Ifit, t, Ifit2)
        return p


class Probe:
    def __init__(self, shn=0, sock=None):
        self.shn, self.sock = shn, sock

        self.IO_mds = self.IO_file = None
        self.nodes = ()
        self.PP = self.IV_series = self.S = None

        self.xlab = "t [s]"
        self.ylab = ("Isat [A]", "Vf [V]", "Te [eV]")

    def __getitem__(self, index):
        if not isinstance(index, tuple):
            return self.S[index]
        else:
            return tuple([self.S[i] for i in index])

    def mapsig(self):
        pass

    def calib(self):
        pass
        
    def load_mds(self):
        self.x = self.IO_mds.load(self.nodes)
        self.mapsig()
        
    def load_file(self):
        self.x = self.IO_file.load(self.nodes)
        self.mapsig()

    def save(self):
        self.IO_file.save(self.x)

    def load(self, plunge=None, trim=True, calib=True, corr_capa=False):
        try:
            self.x = self.IO_file.load(self.nodes)
        except:
            self.x = self.IO_mds.load(self.nodes)
            self.save()

        self.mapsig()
        if trim: 
            self.trim(plunge)
        if calib:
            self.calib()
        if corr_capa:
            self.corr_capa()

    def get_type(self, type):
        istype = lambda x: x.type == type
        return filter(istype, self.S.itervalues())

    def plot_raw(self, **kw):
        if self.S is None:
            self.load(**kw)
        fig = figure()
        types = ['Current', 'Voltage', 'Position']
        
        for i, typ in enumerate(types):
            ax = fig.add_subplot(3, 1, 1+i)
            grid(True)
            ylabel(typ)
            for S in self.get_type(typ):
                S.plot(newfig=False)
        xlabel(self.xlab)

    def trim(self, plunge=None):
        S = self.get_type('Position')
        i0, iM, i1 = S[0].t_ind

        if plunge is None:
            s = np.concatenate(map(np.arange, i0, i1))
        else:
            s = slice(i0[plunge], i1[plunge])

        for S in self.S.itervalues():
            S.trim(s)

    def corr_capa(self):
        for I in self.get_type('Current'):
            I.x[:] -= I.I_capa()

    def calc_IV_series(self, n=1, **kw):
        V = self.S['V']
        II = self.get_type('Current')
        if V.iE is None:
            V.chop_sweeps()
        iE = np.c_[V.iE[:-n], V.iE[n:]]
        return IVSeries(V, II, iE, **kw)

    def analyze(self, **kw):
        self.load(**kw)
        self.IV_series = self.calc_IV_series(engine='fmin')
        self.PP = self.IV_series.fit()

        Isat = self.PP.c[0,:,:,0].T
        dens = np.e*np.sqrt(Isat[0]*Isat[1])
        Mach = 0.5*np.log(Isat[0]/Isat[1])
        
        c = self.PP.c.copy()
        c[0,:,0,0] = dens
        c[0,:,1,0] = Mach
        self.PP_Mach = PiecewisePolynomial(c, self.PP.x, **self.PP.kw)
        
    def analyze2(self, n=2, **kw):
        if self.IV_series is None:
            self.analyze(**kw)

        V = self.S['V']
        II = self.get_type('Current')
        self.IV_series2 = self.calc_IV_series(n=2, engine='fmin')
        self.PP2 = self.IV_series2.fit()

        """
        V = self.S['V']
        II = self.get_type('Current')
        PP = self.PP(V.t)
        mask = self.IV_series.mask()
        
        iE = self.IV_series.iE
        self.IV_series2 = IVSeries2(V, II, PP, mask, iE)
        """

    def plot(self, x=None, PP='PP'):
        if self.PP is None:
            self.analyze()

        if PP == 'PP':
            ylab = self.ylab
        else:
            ylab = ("n [au], M",) + self.ylab[1:]
        PP = getattr(self, PP)

        if x is None:
            xlab = self.xlab
        elif x == 'R':
            x = self['R'].x
            xlab = "R [mm]"

        fig = figure(figsize=(10,10))
        for i, pp in enumerate(PP.T):
            fig.add_subplot(3, 1, 1+i)
            
            pp.plot(False, x=x)
            grid(True)
            ylabel(ylab[i])
        xlabel(xlab)


        fig2 = figure(figsize=(6,5))
        ax = fig2.add_subplot(1, 1, 1)
        ti = self.IV_series.ti
        V_range, I_range = self.IV_series.V_range, self.IV_series.I_range
        corners = (V_range[0], I_range[0]), (V_range[1], I_range[1])
        #ax.dataLim.update_from_data_xy(corners, ignore=True)
        ax.set_xlim(V_range)
        ax.set_ylim(I_range)
        ax.grid(True)

        def plotfun(event):
            t_event = event.xdata
            c = (ti[:,0] <= t_event) & (t_event < ti[:,1])
            IV_group = self.IV_series.IV_group[c][0]

            IV_group.plot(ax, 'get_xy')

        self.MMOV = MouseMotionObserverViewer(fig.axes, ax, plotfun)

