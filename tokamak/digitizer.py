import numpy as np
import os
from pprint import pformat

from collections import OrderedDict, Mapping

import h5py as H5
#import hdf5_cffi as H5

from mdsclient import *

from sm_pyplot.tight_figure import get_fig

from utils.sig import memoized_property, Signal


class IOH5:
    def __init__(self, h5name="test.h5"):
        self.h5name = h5name

    def save(self, d, compression="gzip"):
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
    def __init__(self):
        pass

    def get_size(self, node, **kw):
        pass

    def get_node(self, node, **kw):
        pass

    def put_node(self, node, val):
        pass
    
    def load(self, nodes, **kw):
        if isinstance(nodes, str):
            nodes = (nodes,)
        return OrderedDict([(k, self.get_node(k, **kw)) for k in nodes])

    def save(self, x, nodes=None):
        if nodes is None:
            nodes = x.keys()
        for node in nodes:
            self.put_node(node, x[node])


class IOFile(IO):
    def __init__(self, shn, suffix="", subdir="", group=""):
        self.shn, self.suffix, self.subdir, self.group = shn, suffix, subdir, group

        self.basepath = os.environ['DATAPATH']
        self.fullpath = os.path.join(self.basepath, self.subdir)
        self.fname = str(self.shn) + self.suffix + '.h5'
        
        self.h5name = os.path.join(self.fullpath, self.fname)
        IO.__init__(self)

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
        self._f.create_dataset(name, data=val, compression="gzip")

    @ensure_open("r")
    def load(self, *args, **kw):
        return IO.load(self, *args, **kw)

    @ensure_open("w")
    def save(self, *args):
        IO.save(self, *args)


class TdiError(Exception):
    def __init__(self, err, mdsfmt, args):
        self.err, self.mdsfmt, self.args = err, mdsfmt, args

    def __str__(self):
        return self.err + '\nmdsfmt:\n' + self.mdsfmt + '\nargs:\n' + pformat(self.args)


class IOMds(IO):
    def __init__(self, shn):
        self.shn = shn
        self.mdsserver = "localhost"
        self.mdsport = "8000"
        self.mdstree = None
        self.mdsplaceholder = "$"
        self.mdsfmt = "%s"
        self.datadeco = "data(%s)"
        self.timedeco = "dim_of(%s)"
        self.sizedeco = "size(%s)"
        self.last_node = None
        IO.__init__(self)

    @memoized_property
    def sock(self):
        s = mdsconnect(self.mdsserver + ':' + str(self.mdsport))
        if self.mdstree is not None:
            mdsopen(s, self.mdstree, self.shn)
        return s

    def _mdsslicestr(self, s):
        return '[%s:%s:%s]' % (str(s.start or '*'), str(s.stop or '*'), str(s.step or '*'))

    def _mdsstr(self, node, s=None):
        mdsfmt = self.mdsfmt
        if node == 't':
            node = self.last_node
            if node is None:
                raise RuntimeError("Need to load another node before 't'")
            mdsfmt = self.timedeco % mdsfmt
        else:
            self.last_node = node
            mdsfmt = self.datadeco % mdsfmt
        if s is not None:
            mdsfmt += self._mdsslicestr(s)
        return mdsfmt % node

    def get_size(self, node, t0=None, t1=None, s=None):
        return self.mdsvalue(self.sizedeco % self._mdsstr(node, s), t0, t1)

    def get_node(self, node, t0=None, t1=None, s=None):
        return self.mdsvalue(self._mdsstr(node, s), t0, t1)

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
    def __init__(self, shn=0, name="", s=None, nodes=(), tnode='t', tunits='s',
            IO_mds=None, IO_file=None):
        self.shn, self.name, self.s = shn, name, s
        self.nodes, self.tnode, self.tunits = nodes, tnode, tunits
        self.IO_mds, self.IO_file = IO_mds, IO_file

        self.window = None
        self.perm = dict()
        self.amp = dict()

        if tnode in nodes:
            self.all_nodes = nodes
        else:
            self.all_nodes = nodes[:1] + (tnode,) + nodes[1:]

        IO.__init__(self)

    def __repr__(self):
        fmtstr = "%s {name} for shot {shn}"
        return (fmtstr % self.__class__.__name__).format(**self.__dict__)

    def __str__(self):
        return self.__repr__() + ", with:\n%s" % pformat(self.__dict__)

    @memoized_property
    def x(self):
        return self.load()

    def assignal(self, *args, **kw):
        return Signal(*args, tunits=self.tunits, **kw)

    def __getitem__(self, indx):
        return self.assignal(self.x[indx], self.x[self.tnode], name=indx)

    def keys(self):
        return [k for k in self.x.keys() if k != self.tnode]

    def __len__(self):
        return len(self.keys())

    def __iter__(self):
        for k in self.keys():
            yield k

    def load_raw_mds(self, **kw):
        self.x = self.IO_mds.load(self.all_nodes, **kw)
        for k in self.x.keys():
            self.x[k] = self.x[k].astype(np.float32)
        return self.x
        
    def load_raw_file(self, **kw):
        self.x = self.IO_file.load(self.all_nodes, **kw)
        return self.x

    def get_size(self, node, **kw):
        try:
            return self.IO_file.get_size(node, **kw)
        except (IOError, KeyError):
            kw.setdefault('s', self.s)
            return self.IO_mds.get_size(node, **kw)

    def get_node(self, node, **kw):
        try:
            return self.IO_file.get_node(node, **kw)
        except (IOError, KeyError):
            kw.setdefault('s', self.s)
            val = self.IO_mds.get_node(node, **kw).astype(np.float32)
            self.put_node(node, val)
            return val

    def put_node(self, node, val):
        return self.IO_file.put_node(node, val)
        
    def load_raw(self, **kw):
        self.x = IO.load(self, self.all_nodes, **kw)
        return self.x
    
    def save(self):
        self.IO_file.save(self.x, self.all_nodes)

    def calib(self):
        for node in self.perm:
            self.x[node][:] = self.x[node].transpose(self.perm[node])
        for node in self.amp:
            self.amp[node].apply(self.x[node])
        if self.window is not None:
            for node in self.all_nodes:
                self.x[node] = self.x[node][self.window]

    def _load_calib_factory(name):
        def load_calib(self, **kw):
            getattr(self, name)(**kw)
            self.calib()
            return self.x
        return load_calib

    load_mds  = _load_calib_factory('load_raw_mds')
    load_file = _load_calib_factory('load_raw_file')
    load      = _load_calib_factory('load_raw')

    def calib_offset(self, **kw):
        self.load_raw(**kw)
        offs = [median(self.x[node]) for node in self.nodes]
        return offs

    def plot(self, nodes=None, fig=None):
        if nodes is None:
            nodes = self.keys()

        fig = get_fig(fig, (len(nodes), 1), xlab=self[nodes[0]].xlab, ylab=nodes)

        for node, ax in zip(nodes, fig.axes):
            self[node].plot(ax)
        fig.canvas.draw()
        return fig        

