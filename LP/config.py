import numpy as np

import textwrap

from collections import OrderedDict

from sig import ensure_tuple, Container, PositionSignal, VoltageSignal, CurrentSignal

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
        
    def unique_signals(self, get_x):
        return {k: get_x(k) for k in self.unique_keys}

    def mapsig(self, get_x, t):
        self.unique_sigs = self.unique_signals(get_x)

        R = self.unique_sigs[self.R_keys]
        self.S = OrderedDict(R=PositionSignal(R, t, name='R'))

        nans = np.zeros_like(t)
        nans.fill(np.nan)

        def get_sig(key):
            try:
                return self.unique_sigs[key]
            except KeyError:
                return nans

        for tip in self.tips:
            i = tip.number
            x_V = get_sig(tip.V_keys)
            x_I = get_sig(tip.I_keys)

            V = VoltageSignal(x_V, t, number=i, name='V%d' % i)
            I = CurrentSignal(x_I, t, V=V, number=i, name='I%d' % i)
            self.S[tip.name] = I
        
        return self.S

    def calib(self, get_amp):
        for k in self.unique_sigs.iterkeys():
            amp = get_amp(k)
            amp.apply(self.unique_sigs[k])

    def norm_to_region(self, s):
        for tip in self.tips:
            S = self.S[tip.name]
            S.norm_to_region(s)

            # I_keys is None: floating potential
            if tip.I_keys is None:
                S.V.norm_to_region(s)

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

        return self.head.mapsig(get_x, x['t'].astype('d'))

    def calib(self, line):
        def get_amp(key):
            return self.get(line, 'amp', key)

        self.head.calib(get_amp)
        

class ShotContainer(Container):
    @property
    def shots(self):
        return self.collect_as_list('shn')

    def _collect_as_dict_factory(attr):
        @property
        def wrapper(self):
            return self.collect_as_dict(attr)
        return wrapper

    times = _collect_as_dict_factory('times')
    posit = _collect_as_dict_factory('posit')
    descr = _collect_as_dict_factory('descr')
    stars = _collect_as_dict_factory('stars')

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


class ShotNotFoundError(Exception):
    pass

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

        raise ShotNotFoundError("No config information for shot %d" % shn)
    
    def __repr__(self):
        s = [str(v) + "\n" for v in self.x.itervalues()]
        return np.array(s).tostring()



