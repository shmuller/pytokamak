import numpy as np
import os
import h5py
import copy

from matplotlib.pyplot import figure, plot, draw, hold, xlim, ylim
from tight_figure import tight_figure as figure

import scipy.optimize as opt

from mdsclient import *

from pdb import set_trace

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


class IVChar:
    def __init__(self, V, I):
        ind = V.x.argsort()
        self.V, self.I = V[ind], I[ind]
        self.Vm, self.VM = self.V.range()
        self.Im, self.IM = self.I.x.min(), np.median(self.I.x[:ind.size/2])
        self.dV = self.VM - self.Vm
        self.dI = self.IM - self.Im

    def get_IV(self):
        return self.V.x, self.I.x

    def get_IV_norm(self):
        V = (self.V.x - self.Vm)/self.dV
        I = (self.I.x - self.Im)/self.dI*2 - 1
        return V, I

    def plot(self, newfig=True, fun='get_IV'):
        get_IV = getattr(self, fun)
        if newfig: figure()
        line, = plot(*get_IV())
        return line

    def update(self, line, fun='get_IV'):
        get_IV = getattr(self, fun)
        line.set_data(*get_IV())


class IVGroup:
    def __init__(self, V, II, s=slice(None)):
        self.IV_char = np.empty(len(II), object)
        for j, I in enumerate(II):
            self.IV_char[j] = IVChar(V[s], I[s])

    def __getitem__(self, index):
        return self.IV_char[index]

    def plot(self, newfig=True, fun='get_IV'):
        if newfig: figure()
        return [IV_char.plot(False, fun) for IV_char in self.IV_char]

    def update(self, lines, fun='get_IV'):
        for IV_char, line in zip(self.IV_char, lines):
            IV_char.update(line, fun)


class IVSeries:
    def __init__(self, V, II):
        self.V_range = V.plot_range()
        self.I_range = self._plot_range(II)

        if V.iE is None:
            V.chop_sweeps()
        N = V.iE.size-1
        self.IV_group = np.empty(N, object)
        for j in xrange(N):
            s = slice(V.iE[j], V.iE[j+1]+1)
            self.IV_group[j] = IVGroup(V, II, s)

    def __getitem__(self, index):
        return self.IV_group[index]

    def _plot_range(self, II):
        I_range = np.array([I.plot_range() for I in II])
        return I_range[:,0].min(), I_range[:,1].max()

    def plot(self, fun='get_IV'):
        lines = self.IV_group[0].plot(True, fun)
        if fun == 'get_IV':
            xlim(self.V_range)
            ylim(self.I_range)
        elif fun == 'get_IV_norm':
            xlim(( 0.0, 1.0))
            ylim((-1.2, 1.2))
        draw()
        for IV_group in self.IV_group[1:]:
            IV_group.update(lines, fun)
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

    def IV_series(self):
        V = self.S['V']
        II = self.get_type('Current')
        return IVSeries(V, II)


