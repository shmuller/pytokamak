import numpy as np
import os
from warnings import warn
from pprint import pformat

from collections import OrderedDict, Mapping

from mdsclient import *

from pytokamak.utils.utils import memoized_property, GeneratorDict
from pytokamak.utils.sig import Signal, AmpSignal, get_fig


class IOH5:
    def __init__(self, h5name="test.h5"):
        self.h5name = h5name

    def save(self, d, compression="szip"):
        with H5.File(self.h5name, "w") as f:
            for key, val in d.iteritems():
                f.create_dataset(key, data=val, compression=compression)

    def load(self, d=None):
        with H5.File(self.h5name, "r") as f:
            if d is None: 
                d = dict()
            for k in f:
                d[k] = f[k].value
        return d


class IO:
    def __init__(self, shn, *args, **kw):
        self.shn = shn

    def get_size(self, node, **kw):
        raise NotImplementedError

    def get_node(self, node, **kw):
        raise NotImplementedError

    def put_node(self, node, val):
        raise NotImplementedError
    
    def load(self, nodes, **kw):
        if isinstance(nodes, str):
            nodes = (nodes,)
        return OrderedDict([(k, self.get_node(k, **kw)) for k in nodes])

    def save(self, x, nodes=None):
        if nodes is None:
            nodes = x.keys()
        for node in nodes:
            self.put_node(node, x[node])


class IOFileH5(IO):
    def __init__(self, shn, suffix="", subdir="", group=""):
        IO.__init__(self, shn)
        self.suffix, self.subdir, self.group = suffix, subdir, group
        try:
            self.basepath = os.environ['DATAPATH']
        except KeyError:
            self.basepath = ""
        self.fullpath = os.path.join(self.basepath, self.subdir)
        self.fname = str(shn) + self.suffix + '.h5'
        
        self.h5name = os.path.join(self.fullpath, self.fname)

    def ensure_open(mode):
        def ensure_open_mode(fun):
            def open(self, *args, **kw):
                try:
                    x = fun(self, *args, **kw)
                except:
                    with H5.File(self.h5name, mode) as self._f:
                        x = fun(self, *args, **kw)
                return x
            return open
        return ensure_open_mode

    @ensure_open("r")
    def get_size(self, node, **kw):
        return self._f[self.group + '/' + node].len()

    @ensure_open("r")
    def get_node(self, node, **kw):
        return self._f[self.group + '/' + node].value

    @ensure_open("a")
    def put_node(self, node, val):
        name = self.group + '/' + node
        if name in self._f:
            del self._f[name]
        self._f.create_dataset(name, data=val, compression="szip")

    @ensure_open("a")
    def del_node(self, node):
        name = self.group + '/' + node
        if name in self._f:
            del self._f[name]
    
    @ensure_open("a")
    def del_all(self):
        del self._f[self.group]

    @ensure_open("r")
    def load(self, *args, **kw):
        return IO.load(self, *args, **kw)

    @ensure_open("w")
    def save(self, *args):
        IO.save(self, *args)


try:
    import h5py as H5
    #import hdf5_cffi as H5
except ImportError:
    warn("Could not import h5py, disabling cache files...")
    IOFile = IO
else:
    IOFile = IOFileH5


class IOMds(IO):
    def __init__(self, shn):
        IO.__init__(self, shn)
        self.mdsserver = "localhost"
        self.mdsport = "8000"
        self.mdstree = None
        self.mdsplaceholder = "$"
        self.mdsfmt = "%s"
        self.datadeco = "data(%s)"
        self.timedeco = "dim_of(%s)"
        self.sizedeco = "size(%s)"
        self.last_mdsbasestr = None

    @memoized_property
    def sock(self):
        s = mdsconnect(self.mdsserver + ':' + str(self.mdsport))
        if self.mdstree is not None:
            mdsopen(s, self.mdstree, self.shn)
        return s

    def get_mdsbasestr(self, node):
        return self.mdsfmt % node

    def _mdsslicestr(self, s):
        return '[%s:%s:%s]' % (str(s.start or '*'), str(s.stop or '*'), str(s.step or '*'))

    def get_mdsstr(self, node, s=None):
        if node == 't':
            mdsbasestr = self.last_mdsbasestr
            if mdsbasestr is None:
                raise RuntimeError("Need to load another node before 't'")
            mdsstr = self.timedeco % mdsbasestr
        else:
            self.last_mdsbasestr = mdsbasestr = self.get_mdsbasestr(node)
            mdsstr = self.datadeco % mdsbasestr
        if s is not None:
            mdsstr += self._mdsslicestr(s)
        return mdsstr

    def get_size(self, node, t0=None, t1=None, s=None):
        return self.mdsvalue(self.sizedeco % self.get_mdsstr(node, s), t0, t1)

    def get_node(self, node, t0=None, t1=None, s=None):
        return self.mdsvalue(self.get_mdsstr(node, s), t0, t1)

    def save(self, x):
        raise NotImplementedError("Saving to MDS not implemented")

    def _fix_args(self, mdsfmt, orig_args):
        """
        Replace mdsvalue placeholder (typically '$') by '*' for each argument
        that is None. Arguments than exceed the number of placeholders are ignored.
        """
        parts, args = [], []
        for arg in orig_args:
            head, sep, mdsfmt = mdsfmt.partition(self.mdsplaceholder)
            parts.append(head)
            if len(sep) == 0:
                break
            if arg is None:
                sep = '*'
            else:
                args.append(arg)
            parts.append(sep)

        parts.append(mdsfmt)
        mdsfmt = ''.join(parts)
        return mdsfmt, args

    def mdsvalue(self, mdsfmt, *args):
        mdsfmt, args = self._fix_args(mdsfmt, args)
        print mdsfmt
        ret = mdsvalue(self.sock, mdsfmt, *args)
        if isinstance(ret, str) and (ret.startswith("Tdi") or ret.startswith('%')):
            raise TdiError(ret, mdsfmt, args)
        return ret


class Digitizer(IO, Mapping):
    def __init__(self, shn=0, name="", t0=None, t1=None, s=None, 
                 nodes=(), tnode='t', tunits='s', IO_mds=None, IO_file=None):
        IO.__init__(self, shn)
        
        self.name, self.t0, self.t1, self.s = name, t0, t1, s
        self.nodes = nodes
        self.tnode, self.tunits = tnode, tunits
        self.IO_mds, self.IO_file = IO_mds, IO_file

        self.amp = dict()
        self.x = GeneratorDict(generator=self.get_node)

    @property
    def loaded_nodes(self):
        return [k for k in self.nodes if k in self.x]

    @memoized_property
    def window(self):
        return None

    def __repr__(self):
        fmtstr = "%s {name} for shot {shn}"
        return (fmtstr % self.__class__.__name__).format(**self.__dict__)

    def __str__(self):
        return self.__repr__() + ", with:\n%s" % pformat(self.__dict__)

    def keys(self):
        return self.nodes

    def __len__(self):
        return len(self.keys())

    def __iter__(self):
        for k in self.keys():
            yield k

    def get_size(self, node, **kw):
        try:
            return self.IO_file.get_size(node, **kw)
        except (NotImplementedError, IOError, KeyError):
            kw.setdefault('s', self.s)
            return self.IO_mds.get_size(node, **kw)

    def get_node(self, node, **kw):
        try:
            x = self.IO_file.get_node(node, **kw)
        except (NotImplementedError, IOError, KeyError):
            kw.setdefault('t0', self.t0)
            kw.setdefault('t1', self.t1)
            kw.setdefault('s', self.s)
            x = self.IO_mds.get_node(node, **kw)
            self.put_node(node, x)
        return x

    @memoized_property
    def window(self):
        return None

    def make_sig(self, node, x, t):
        if self.window is not None:
            x = x[self.window]
            t = t[self.window]
        if node in self.amp:
            amp = self.amp[node]
        else:
            amp = None
        return AmpSignal(x, t, amp=amp, dtype=np.float64, name=node)

    def __getitem__(self, node):
        if node not in self.nodes:
            raise KeyError(node)
        return self.make_sig(node, self.x[node], self.x[self.tnode])
        
    def put_node(self, node, val):
        try:
            self.IO_file.put_node(node, val)
        except NotImplementedError:
            pass

    def load(self):
        # visit all signals to load them
        self.values()
        return self
 
    def plot(self, nodes=None, fig=None):
        if nodes is None:
            nodes = self.loaded_nodes

        fig = get_fig(fig, (len(nodes), 1), xlab=self[nodes[0]].xlab, ylab=nodes)

        for node, ax in zip(nodes, fig.axes):
            self[node].plot(ax)
        fig.canvas.draw()
        return fig        

