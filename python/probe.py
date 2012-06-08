import numpy as np
import os
import h5py
from mdsclient import *

from matplotlib.pyplot import figure, plot, show
from tight_figure import tight_figure as figure

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

    def cut_tail(self):
        s = np.where(self.R < 0.99*self.R[0]);



    def plot(self):
        X = np.c_[self.I,self.R,self.V/100.]
        figure()
        plot(self.t,X)



