import numpy as np
import os
import h5py
import copy

from matplotlib.pyplot import figure, plot, draw, hold, xlim, ylim
from tight_figure import tight_figure as figure

import scipy.optimize as opt
import scipy.odr as odr

from mdsclient import *

from ipdb import set_trace

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
    def __init__(self, shn=0, subdir=""):
        self.shn, self.subdir = shn, subdir

        pth = os.environ['DATAPATH']
        pth = os.path.join(pth, self.subdir)
        fname = str(self.shn) + '.h5'
        
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
                raise ValueError("Need other node name to obtain 't'")
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


class PeriodPhaseFinder():
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
        iE = np.round(np.arange(i,M,d/2)).astype('i')
        
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
        return self.__init__(self.t[index], self.x[index], self.name, self.type)

    def len(self):
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

    def plot(self, newfig=True):
        if newfig: figure()
        plot(self.t, self.x.T)


class PositionSignal(Signal):
    def __init__(self, t, x, name=""):
        Signal.__init__(self, t, x, name, 'Position')

    def __getitem__(self, index):
        return PositionSignal(self.t[index], self.x[index], self.name)

    def get_t_ind(self):
        i0 = np.argmax(self.x)
        i1 = np.where(self.x[i0:] < self.x[0])[0][0]
        i1 += i0
        return i0, i1


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
    def __init__(self, t, x, name="", V=None, C=None):
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

    return opt.fmin(dy2, p0)

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

    def fit(self):
        if not self.is_OK():
            raise RuntimeError("Cannot fit data that failed is_OK() check")
            
        self.P0 = self.get_guess()
        self.P = self.engine(self.fitfun, self.P0, self.X, self.Y)
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

    def plot(self, newfig=True, fun='get_xy', lines=None):
        x, y = getattr(self, fun)()
        
        if lines is None: lines = []
        Nl, Ny = len(lines), len(y)
        if newfig and Nl == 0: figure()

        for li, yi in zip(lines, y):
            li.set_data(x, yi)
            li.set_visible(True)

        for li in lines[Ny:]:
            li.set_visible(False)

        for yi in y[Nl:]:
            li, = plot(x, yi, linewidth=1.5)
            lines.append(li)

        return lines


class FitterIV(Fitter):
    def __init__(self, V, I, cut_at_min=True, **kwargs):
        ind = V.argsort()
        self.V, self.I = V[ind], I[ind]

        self.im = self.I.argmin()
        self.Vm, self.VM = self.V[0], self.V[-1]
        self.Im, self.IM = self.I[self.im], np.median(self.I[:I.size/2])
        self.dV = self.VM - self.Vm
        self.dI = self.IM - self.Im

        self.cut_at_min = cut_at_min

        Fitter.__init__(self, self.V, self.I, **kwargs)

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
            M = self.im+1
        else:
            M = self.V.size
        self.X = (self.V[:M].astype('d') - self.Vm)/self.dV
        self.Y = (self.I[:M].astype('d') - self.Im)/self.dI*2 - 1
        return self.X, self.Y

    def set_unnorm(self):
        def unnormalize(P, Vm, VM, Im, IM):
            dV, dI = VM-Vm, IM-Im
            a = 2./dI
            b = -(Im*a+1)

            p = np.empty_like(P)
            p[0] = (P[0]-b)/a
            p[1] = (P[1]+P[2]*np.log(1.-b/P[0]))*dV+Vm
            p[2] = P[2]*dV
            return p

        self.p0 = unnormalize(self.P0, self.Vm, self.VM, self.Im, self.IM)
        self.p  = unnormalize(self.P , self.Vm, self.VM, self.Im, self.IM)
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


class IVChar:
    def __init__(self, V, I, **kwargs):
        self.V, self.I = V, I
        self.fitter_IV = FitterIV(self.V.x, self.I.x, **kwargs)

    def fit(self, out=None):
        if out is None:
            out = np.empty(3)
        try:
            out[:] = self.fitter_IV.fit()
        except RuntimeError:
            out[:] = np.nan
        return out
        
    def plot(self, newfig=True, fun='get_xy', lines=None):
        return self.fitter_IV.plot(newfig, fun, lines)


class IVGroup:
    def __init__(self, V, II, s=slice(None), **kwargs):
        self.IV_char = np.empty(len(II), object)
        for j, I in enumerate(II):
            self.IV_char[j] = IVChar(V[s], I[s], **kwargs)

    def __getitem__(self, index):
        return self.IV_char[index]

    def fit(self, out=None):
        if out is None:
            out = np.empty((self.IV_char.size, 3))
        for p, IV_char in zip(out, self.IV_char):
            IV_char.fit(p)
        return out

    def plot(self, newfig=True, fun='get_xy', lines=None):
        if lines is None:
            if newfig: figure()
            return [x.plot(False, fun) for x in self.IV_char]
        else:
            return [x.plot(False, fun, l) for x, l in zip(self.IV_char, lines)]


class IVSeries:
    def __init__(self, V, II, **kwargs):
        self.V_range = V.plot_range()
        self.I_range = self._plot_range(II)

        if V.iE is None:
            V.chop_sweeps()
        N = V.iE.size-1
        self.siz = (N, len(II))

        self.IV_group = np.empty(N, object)
        for j in xrange(N):
            s = slice(V.iE[j], V.iE[j+1]+1)
            self.IV_group[j] = IVGroup(V, II, s, **kwargs)

    def __getitem__(self, index):
        return self.IV_group[index]

    def _plot_range(self, II):
        I_range = np.array([I.plot_range() for I in II])
        return I_range[:,0].min(), I_range[:,1].max()

    def fit(self):
        out = np.empty(self.siz + (3,))
        for p, IV_group in zip(out, self.IV_group):
            IV_group.fit(p)
        return out

    def plot(self, fun='get_xy'):
        figure()
        if fun == 'get_xy':
            xlim(self.V_range)
            ylim(self.I_range)
        else:
            xlim(( 0.0, 1.0))
            ylim((-1.2, 1.2))

        lines = None
        for IV_group in self.IV_group:
            lines = IV_group.plot(False, fun, lines)
            draw()


class Probe:
    def __init__(self, shn=0, sock=None):
        self.shn, self.sock = shn, sock

        self.IO_mds = self.IO_file = None
        self.nodes = ()

    def mapsig(self):
        pass
        
    def load_mds(self):
        self.x = self.IO_mds.load(self.nodes)
        self.mapsig()
        
    def load_file(self):
        self.x = self.IO_file.load(self.nodes)
        self.mapsig()

    def save(self):
        self.IO_file.save(self.x)

    def load(self):
        try:
            self.x = self.IO_file.load(self.nodes)
        except:
            self.x = self.IO_mds.load(self.nodes)
            self.save()

        self.mapsig()

    def get_type(self, type):
        istype = lambda x: x.type == type
        return filter(istype, self.S.itervalues())

    def plot(self):
        figure()
        hold(True)
        for S in self.S.itervalues():
            S.plot(newfig=False)

    def trim(self):
        S = self.get_type('Position')
        i0, i1 = S[0].get_t_ind()

        s = slice(i1)
        for S in self.S.itervalues():
            S.trim(s)

    def IV_series(self, **kwargs):
        V = self.S['V']
        II = self.get_type('Current')
        return IVSeries(V, II, **kwargs)


