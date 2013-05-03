import numpy as np

import textwrap

from pprint import pformat

from sig import ensure_tuple, Container

from matplotlib.path import Path
from matplotlib.patches import PathPatch

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

    def __repr__(self):
        fmtstr = "%s {number}, '{name}', {pos}, connected to '{V_keys}' and '{I_keys}'"
        return (fmtstr % self.__class__.__name__).format(**self.__dict__)

    def __str__(self):
        return self.__repr__() + ", with:\n%s" % pformat(self.__dict__)


class CylindricalTip(Tip):
    def __init__(self, r, z, *args, **kw):
        self.r, self.z = r, z
        area = 2.*np.pi*r*z
        proj_area = 2*r*z
        Tip.__init__(self, area, proj_area, *args, **kw)


class Head:
    def __init__(self, tips, R_keys=None, d=None):
        self.tips, self.R_keys, self.d = tips, R_keys, d

        zp, zm = d/2., -d/2.
        self.xy = np.array([(0., zp), (0., zm), (2.5, zm), (2.5, zp)])
   
    def get_tip_number_by_position(self, pos):
        for tip in self.tips:
            if tip.pos == pos:
                return tip.number

    def __repr__(self):
        return "%s with %d tips" % (self.__class__.__name__, len(self.tips))

    def __str__(self):
        return self.__repr__() + ":\n%s" % pformat(self.tips)

    def plot(self, ax, R, z, **kw):
        kw.setdefault('facecolor', 'k')
        kw.setdefault('edgecolor', 'none')

        xy = self.xy.copy()
        xy[:2,0] += R
        xy[:,1] += z
        pp = PathPatch(Path(xy), **kw)
        ax.add_patch(pp)
        return ax


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

    def __str__(self):
        return self.__repr__() + ":\n%s" % self.descr


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



