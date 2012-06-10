import numpy as np
import os
import h5py
import copy

from matplotlib.pyplot import figure, plot, show, hold
from tight_figure import tight_figure as figure

import scipy.optimize as opt
from mdsclient import *

from pdb import set_trace

class Data(dict):
    def __init__(self, nodes, X):
        self.nodes, self.X = nodes, X
        p = [(nodes[i], X[i]) for i in xrange(len(nodes))]
        dict.__init__(self, p)

    def view(self, s):
        p = [(self.nodes[i], self.X[i,s]) for i in xrange(len(self.nodes))]
        return dict(p)


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
        
        N = iE.size/2
        iE = iE[:2*N].reshape(N,2).T        

        dx = self.x[iE[0]] - self.x[iE[1]]
        D = -np.sqrt(dx.dot(dx)/N)

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

    def copy(self):
        s = copy.copy(self)
        s.x = s.x.copy()
        return s

    def trim(self, s):
        self.t, self.x = self.t[s], self.x[s]

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

    def get_t_ind(self):
        i0 = np.argmax(self.x)
        i1 = np.where(self.x[i0:] < self.x[0])[0][0]
        i1 += i0
        return i0, i1


class VoltageSignal(Signal):
    def __init__(self, t, x, name=""):
        Signal.__init__(self, t, x, name, 'Voltage')
        self.iE = None

    def chop_sweeps(self):
        PPF = PeriodPhaseFinder(self.x)
        PPF.chop_sweeps()
        self.iE = PPF.get_iE()

    def plot(self):
        Signal.plot(self)
        if self.iE is not None:
            plot(self.t[self.iE], self.x[self.iE], 'r+')


class CurrentSignal(Signal):
    def __init__(self, t, x, name="", V=None, C=None):
        Signal.__init__(self, t, x, name, 'Current')
        self.V, self.C = V, C

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

    def plot(self):
        figure()
        hold(True)
        for S in self.S.itervalues():
            S.plot(newfig=False)

    def trim(self):
        for S in self.S.itervalues():
            if isinstance(S, PositionSignal):
                i0, i1 = S.get_t_ind()
                break

        s = slice(i1)
        for S in self.S.itervalues():
            S.trim(s)



