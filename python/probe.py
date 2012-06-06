import numpy as np
import os
import h5py
from mdsclient import *

from matplotlib.pyplot import figure, plot, show
from tight_figure import tight_figure as figure


class IOMds:
    def __init__(self, shn=0, sock=None):
        self.shn, self.sock = shn, sock
        self.mdsport, self.mdsfmt = "8000", ""

    def mdsstr(self, node):
        return self.mdsfmt % (self.shn, node)

    def getsig(self, node):
        return mdsvalue(self.sock,self.mdsstr(node))

    def gettim(self, node):
        return mdsvalue(self.sock,'dim_of(%s)' % self.mdsstr(node))

    def load(self, nodes):
        if self.sock is None:
            self.sock = mdsconnect('localhost:' + str(self.mdsport))

        t = self.gettim(nodes[0])
        x = {}
        for node in nodes:
            x[node] = self.getsig(node)

        return t, x

    def save(self, t, x):
        pass


class IOFile:
    def __init__(self, subdir="", shn=0):
        self.subdir, self.shn = subdir, shn

        pth = os.environ['DATAPATH']
        pth = os.path.join(pth, self.subdir)
        fname = str(self.shn) + '.h5'
        
        self.h5name = os.path.join(pth, fname)

    def load(self, nodes):
        f = h5py.File(self.h5name,"r")
    
        t = f['t'].value
        x = {}
        for node in nodes:
            x[node] = f[node].value
        
        f.close()
        return t, x

    def save(self, t, x):
        f = h5py.File(self.h5name,"w")

        f.create_dataset('t',data=t,compression="gzip")

        for node, val in x.iteritems():
            f.create_dataset(node,data=val,compression="gzip")

        f.close()


class Probe:
    def __init__(self, shn=0, sock=None):
        self.shn, self.sock = shn, sock

        self.IO_mds = self.IO_file = None
        self.nodes = ()
        
    def mapsig(self):
        pass

    def load_mds(self):
        self.t, self.x = self.IO_mds.load(self.nodes)
        self.mapsig()
        
    def load_file(self):
        self.t, self.x = self.IO_file.load(self.nodes)
        self.mapsig()

    def save(self):
        self.IO_file.save(self.t,self.x)

    def load(self):
        try:
            self.t, self.x = self.IO_file.load(self.nodes)
        except:
            self.t, self.x = self.IO_mds.load(self.nodes)
            self.save()

        self.mapsig()        

    def plot(self):
        X = np.c_[self.I1,self.I2,self.R,self.Vb/100.]
        figure()
        plot(self.t,X)



