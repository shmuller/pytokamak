import numpy as np
import os
import h5py
from mdsclient import *

from matplotlib.pyplot import figure, plot, show
from tight_figure import tight_figure as figure

class IO:
    def __init__(self):
        pass

    def get_size(self, node):
        pass

    def get_node(self, node):
        pass

    def put_node(self, node, val):
        pass

    def set_view(self, s):
        nodes = self.x.keys()
        Xs = self.X[:,s]
        for i in xrange(len(nodes)):
            self.x[nodes[i]] = Xs[i]

    def load(self, nodes):
        self.nodes = nodes
        M = len(nodes)
        N = self.get_size(nodes[0])

        self.X = np.empty((M,N),'f')
        self.x = {}

        for i in xrange(M):
            node = nodes[i]
            self.X[i,:] = self.get_node(node)
            self.x[node] = self.X[i]

    def save(self):
        for node, val in self.x.iteritems():
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
    
        IO.load(self, nodes)

        self._f.close()
        return self.x

    def save(self):
        self._f = h5py.File(self.h5name,"w")

        IO.save(self)

        self._f.close()


class IOMds(IO):
    def __init__(self, shn=0, sock=None):
        self.shn, self.sock = shn, sock
        self.mdsport, self.mdsfmt = "8000", ""

    def mdsstr(self, node):
        return self.mdsfmt % (self.shn, node)

    def get_size(self, node):
        return mdsvalue(self.sock,'size(%s)' % self.mdsstr(node))

    def get_time(self, node):
        return mdsvalue(self.sock,'dim_of(%s)' % self.mdsstr(node))

    def get_node(self, node):
        if node == 't':
            node = self.nodes[self.nodes.index('t')-1]
            return self.get_time(node)
        else:
            return mdsvalue(self.sock,self.mdsstr(node))

    def load(self, nodes):
        if self.sock is None:
            self.sock = mdsconnect('localhost:' + str(self.mdsport))

        IO.load(self, nodes)

        return self.x

    def save(self):
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



