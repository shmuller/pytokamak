import numpy as np

import textwrap

from collections import OrderedDict

from probe import PositionSignal, VoltageSignal, CurrentSignal

def ensure_tuple(d, *keys):
    for k in keys:
        if not isinstance(d[k], tuple):
            d[k] = (d[k],)
    
class Tip:
    def __init__(self, area, proj_area, number, pos, 
            V_keys=None, I_keys=None, name=None):
        self.area = area
        self.proj_area = proj_area
        self.number = number
        self.pos = pos
        self.V_keys = V_keys
        self.I_keys = I_keys
        self.name = name or 'tip%d' % self.number


class CylindricalTip(Tip):
    def __init__(self, r, z, *args, **kw):
        self.r, self.z = r, z
        area = 2.*np.pi*r*z
        proj_area = 2*r*z
        Tip.__init__(self, area, proj_area, *args, **kw)


class Head:
    def __init__(self, tips, R_keys=None):
        self.tips = tips
        self.R_keys = R_keys
        self.unique_sigs = None

    @staticmethod
    def _unique(x):
        s = set(x)
        s.discard(None)
        return list(s)

    @property
    def V_keys(self):
        return [tip.V_keys for tip in self.tips]

    @property
    def I_keys(self):
        return [tip.I_keys for tip in self.tips]

    @property
    def all_keys(self):
        return [self.R_keys] + self.V_keys + self.I_keys

    @property
    def unique_V_keys(self):
        return self._unique(self.V_keys)

    @property
    def unique_I_keys(self):
        return self._unique(self.I_keys)

    @property
    def unique_keys(self):
        return [self.R_keys] + self.unique_V_keys + self.unique_I_keys
        
    @property
    def unique_types(self):
        return [PositionSignal] \
             + [VoltageSignal]*len(self.unique_V_keys) \
             + [CurrentSignal]*len(self.unique_I_keys)

    @property
    def unique_keys_types(self):
        return zip(self.unique_keys, self.unique_types)

    def unique_signals(self, get_x, t):
        return {k: S(get_x(k), t) for k, S in self.unique_keys_types}

    """
    def signals(self, get_x, t):
        unique_sigs = self.unique_signals(get_x, t)

        S = OrderedDict(R=)

        for tip in self.tips:
            V = tip.V_keys or get_x(V)
            x 

            V = V or VoltageSignal(get_x(V), t)
            I = I or CurrentSignal(get_x(I), t, V=V)
                

            if keyI is None:
                S[V.name] = V
            else:
                I = self.unique_sigs[keyI]
                I.name = 'I%d' % i
                I.V = V
                S[I.name] = I
    """

    def get_tip_number_by_position(self, pos):
        for tip in self.tips:
            if tip.pos == pos:
                return tip.number


def recursive_dictcopy(d):
    if isinstance(d, dict):
        d = d.copy()
        for k, v in d.iteritems():
            d[k] = recursive_dictcopy(v)
    return d


class Shot:
    def __init__(self, comment="", **kw):
        self.comment = comment

        self.descr = textwrap.dedent(kw.pop('descr', ''))
        self.stars = kw.pop('stars', '*')

        self.attrs = ('expt', 'shn', 'dig', 'head', 'amp_default', 
                'lines', 'times', 'posit')
        
        kw.setdefault('times', ())
        kw.setdefault('posit', ())

        for attr in self.attrs:
            setattr(self, attr, kw.pop(attr, None))

        ensure_tuple(self.__dict__, 'times', 'posit')

        self.amp_default = self.amp_default.copy()
        for k in set(self.amp_default.keys()) & set(kw.keys()):
            self.amp_default[k] = kw.pop(k)

        self.lines = recursive_dictcopy(self.lines)
        for k, v in kw.iteritems():
            line, what, key = k.split('_')
            self.lines[line][what][key] = v
                
    def copy(self, comment="", **kw):
        for attr in self.attrs:
            kw.setdefault(attr, getattr(self, attr))

        return self.__class__(comment, **kw)

    def get(self, line, what, key):
        if key is None:
            return None
        elif what == 'amp':
            return self.lines[line]['amp'].get(key, self.amp_default[key])
        else:
            return self.lines[line][what][key]

    def __repr__(self):
        return "%d %-5s: %s" % (self.shn, self.stars, self.comment)

    def print_descr(self):
        print self.descr

    def mapsig(self, x, line):
        def get_x(key):
            return x[self.get(line, 'mapping', key)].astype('d')

        self.unique_sigs = {k: get_x(k) for k in self.head.unique_keys}

        t = x['t']
        R = self.unique_sigs[self.head.R_keys]
        S = OrderedDict(R=PositionSignal(R, t, name='R'))

        for tip in self.head.tips:
            i, keyV, keyI = tip.number, tip.V_keys, tip.I_keys
            if keyV is None:
                V = None
            else:
                V = VoltageSignal(self.unique_sigs[keyV], t, number=i, name='V%d' % i)
            if keyI is None:
                S[V.name] = V
            else:
                I = CurrentSignal(self.unique_sigs[keyI], t, V=V, number=i, name='I%d' % i) 
                S[I.name] = I

        return S

    def mapsig2(self, x, line):
        def get_x(key):
            return x[self.get(line, 'mapping', key)].astype('d')

        t = x['t']
        self.unique_sigs = self.head.unique_signals(get_x, t)

        S = OrderedDict(R=self.unique_sigs[self.head.R_keys])

        for tip in self.head.tips:
            i, keyV, keyI = tip.number, tip.V_keys, tip.I_keys
            if keyV is None:
                V = None
            else:
                V = self.unique_sigs[keyV]
                V.name = 'V%d' % i
            if keyI is None:
                S[V.name] = V
            else:
                I = self.unique_sigs[keyI]
                I.name = 'I%d' % i
                I.V = V
                S[I.name] = I

        return S

    def calib(self, S, line):
        for k in self.unique_sigs.iterkeys():
            amp = self.get(line, 'amp', k)
            amp.apply(self.unique_sigs[k])


class Container:
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


class ShotContainer(Container):
    @property
    def shots(self):
        return self.collect_as_list('shn')

    @property
    def times(self):
        return self.collect_as_dict('times')

    @property
    def posit(self):
        return self.collect_as_dict('posit')

    @property
    def descr(self):
        return self.collect_as_dict('descr')

    @property
    def stars(self):
        return self.collect_as_dict('stars')

    def shots_with_min_stars(self, min_stars='***'):
        cnd = lambda v: v.stars >= min_stars
        return self.collect_as_list('shn', cnd)


class Experiment(ShotContainer):
    def __init__(self, date=None, campaign=None, ShotClass=Shot):
        self.date, self.campaign, self.ShotClass = date, campaign, ShotClass
        ShotContainer.__init__(self)

    def add(self, shn, *args, **kw):
        self.x[shn] = self.ShotClass(*args, expt=self, shn=shn, **kw)
        
    def rep(self, shn, shn0, comment="", **kw):
        self.x[shn] = self.x[shn0].copy(comment, shn=shn, **kw)

    def __repr__(self):
        s = ["  " + str(v) + "\n" for v in self.x.itervalues()]
        return self.date + ":\n" + np.array(s).tostring()


class Campaign(ShotContainer):
    def __init__(self, ExperimentClass=Experiment):
        self.ExperimentClass = ExperimentClass
        ShotContainer.__init__(self)
    
    def add_experiment(self, *args, **kw):
        E = self.ExperimentClass(*args, campaign=self, **kw)
        self.x[kw['date']] = E
        return E

    def find_shot(self, shn):
        for E in self.x.itervalues():
            try: 
                return E[shn]
            except KeyError:
                pass
    
    def __repr__(self):
        s = [str(v) + "\n" for v in self.x.itervalues()]
        return np.array(s).tostring()



