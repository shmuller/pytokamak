import numpy as np
from time import time
from collections import MutableMapping, Iterable, OrderedDict

_tstart_stack = []

def tic():
    _tstart_stack.append(time())

def toc(fmt="Elapsed: %s s"):
    print fmt % (time() - _tstart_stack.pop())


class Timer(object):
    def __init__(self, name=None):
        if name:
            self.fmt = '[%s] Elapsed: %%s' % name
        else:
            self.fmt = 'Elapsed: %s'
        self.tstart = time()

    def toc(self):
        print self.fmt % (time() - self.tstart)

    def __enter__(self):
        pass

    def __exit__(self, type, value, traceback):
        self.toc()


def ensure_tuple(d, *keys):
    for k in keys:
        if not isinstance(d[k], tuple):
            d[k] = (d[k],)


class memoized_property(object):
    """A read-only @property that is only evaluated once."""
    def __init__(self, fget, doc=None):
        self.fget = fget
        self.__doc__ = doc or fget.__doc__
        self.__name__ = fget.__name__

    def __get__(self, obj, cls):
        if obj is None:
            return self
        obj.__dict__[self.__name__] = result = self.fget(obj)
        return result


class ArrayView(np.ndarray):
    def __new__(subtype, x, fields):
        dtype = {f: x.dtype.fields[f] for f in fields}
        return np.ndarray.__new__(subtype, x.shape, dtype,
                                  buffer=x, strides=x.strides)


class DictView(MutableMapping):
    def __init__(self, base, valid_keys):
        self.base, self.valid_keys = base, valid_keys

    def __getitem__(self, key):
        if key in self.valid_keys:
            return self.base[key]
        else:
            raise KeyError(key)

    def __len__(self):
        return len(self.valid_keys)

    def __iter__(self):
        for key in self.valid_keys:
            yield key

    def __setitem__(self, key, value):
        if key in self.valid_keys:
            self.base[key] = value
        else:
            raise KeyError(key)

    def __delitem__(self, key):
        self.valid_keys.remove(key)

    def __repr__(self):
        return self.__class__.__name__ + " with: " + self.valid_keys.__repr__()


class GeneratorDict(dict):
    def __init__(self, generator):
        self.generator = generator
        dict.__init__(self)

    def __missing__(self, indx):
        self[indx] = result = self.generator(indx)
        return result


def recursive_dictcopy(d):
    if isinstance(d, dict):
        d = d.copy()
        for k, v in d.iteritems():
            d[k] = recursive_dictcopy(v)
    return d


class rdict(dict):
    def copy(self):
        # Can't use dict.copy(), since it returns dict not rdict
        d = rdict(self)
        for k, v in d.iteritems():
            if isinstance(v, rdict):
                d[k] = v.copy()
        return d

    def mod(self, **kw):
        for k, v in kw.iteritems():
            keys = k.split('_')
            cur = self
            for key in keys[:-1]:
                cur = cur[key]
            cur[keys[-1]] = v
        return self

    def rep(self, **kw):
        return self.copy().mod(**kw)


class Container(Iterable):
    def __init__(self):
        self.x = OrderedDict()

    def __getitem__(self, indx):
        return self.x[indx]

    def __iter__(self):
        return self.x.itervalues()

    def __add__(self, other):
        s = self.__class__()
        s.x = self.x.copy()
        s.x.update(other.x)
        return s

    def __iadd__(self, other):
        self.x.update(other.x)
        return self

    def __repr__(self):
        return self.__class__.__name__ + " with: " + self.x.keys().__repr__()

    @staticmethod
    def _item(v, attr, cnd):
        x = np.array([getattr(v, attr)])
        if cnd(v):
            return x
        else:
            return np.empty((0,), x.dtype)

    def collect_as_list(self, attr, cnd=lambda v: True): 
        return np.concatenate([v.collect_as_list(attr, cnd) if isinstance(v, Container)
            else self._item(v, attr, cnd) for v in self.x.itervalues()])

    def collect_as_dict(self, attr):
        return {k: getattr(v, attr) for k, v in self.x.iteritems()}


class BoundingBox:
    def __init__(self, x0, x1=None):
        if x1 is None:
            self.x0, self.x1 = x0
        else:
            self.x0, self.x1 = x0, x1

    def copy(self):
        return self.__class__(self.x0.copy(), self.x1.copy())

    def isin(self, x):
        return np.all((self.x0 <= x) & (x <= self.x1))


def dict_pop(d, **kw):
    return {k: d.pop(k, kw[k]) for k in kw.keys()}

def dict_pop_if_present(d, keys):
    return {k: d.pop(k) for k in set(d.keys()) & set(keys)}


