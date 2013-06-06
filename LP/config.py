import numpy as np
import numpy.ma as ma

import textwrap

from pprint import pformat

from sig import ensure_tuple, recursive_dictcopy, Container

from matplotlib.path import Path
from matplotlib.patches import PathPatch

class Tip:
    def __init__(self, area, proj_area, number, pos, label=None):
        self.area = area
        self.proj_area = proj_area
        self.number = number
        self.pos = pos
        self.name = 'tip%d' % number
        self.label = label or self.name

    def __repr__(self):
        fmtstr = "%s {number}, '{name}', {pos}"
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
    def __init__(self, tips, d=None):
        self.tips, self.d = tips, d

        zp, zm = d/2., -d/2.
        self.xy = np.array([(0., zp), (0., zm), (2.5, zm), (2.5, zp)])
   
    def get_tip_by_name(self, name):
        for tip in self.tips:
            if tip.name == name:
                return tip

    def get_tip_number_by_position(self, pos):
        for tip in self.tips:
            if tip.pos == pos:
                return tip.number

    def __repr__(self):
        return "%s with %d tips" % (self.__class__.__name__, len(self.tips))

    def __str__(self):
        return self.__repr__() + ":\n%s" % pformat(self.tips)

    def as_path_patch(self, R, z, **kw):
        kw.setdefault('facecolor', 'k')
        kw.setdefault('edgecolor', 'none')

        xy = self.xy.copy()
        if R is ma.masked:
            R = np.nan

        xy[:2,0] += R
        xy[:,1] += z
        return PathPatch(Path(xy), **kw)

    def plot(self, ax, R, z, **kw):
        ax.add_patch(self.as_path_patch(R, z, **kw))
        return ax


class Shot:
    def __init__(self, comment="", **kw):
        self.comment = comment

        self.descr = textwrap.dedent(kw.pop('descr', ''))
        self.stars = kw.pop('stars', '*')

        self.attrs = ('expt', 'shn', 'dig', 'head', 'amp_default', 
                'lines', 'times', 'posit', 'tipmap')
        
        kw.setdefault('times', ())
        kw.setdefault('posit', ())

        for attr in self.attrs:
            setattr(self, attr, kw.pop(attr, None))

        ensure_tuple(self.__dict__, 'times', 'posit')

        self.tipmap = self.tipmap.copy()

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
        s = ["  " + repr(v) + "\n" for v in self.x.itervalues()]
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
        s = [repr(v) + "\n" for v in self.x.itervalues()]
        return np.array(s).tostring()



