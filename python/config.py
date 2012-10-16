import numpy as np

from collections import OrderedDict

from probe import PositionSignal, VoltageSignal, CurrentSignal

class Tip:
    def __init__(self, area, proj_area, number, pos, V_keys=None, I_keys=None):
        self.area = area
        self.proj_area = proj_area
        self.number = number
        self.pos = pos
        self.keys = dict(V=V_keys, I=I_keys)


class CylindricalTip(Tip):
    def __init__(self, r, z, *args, **kw):
        self.r, self.z = r, z
        area = 2.*np.pi*r*z
        proj_area = 2*r*z
        Tip.__init__(self, area, proj_area, *args, **kw)


class Head:
    def __init__(self, tips, R_keys=None):
        self.tips = tips
        self.keys = dict(R=R_keys)


def recursive_dictcopy(d):
    if isinstance(d, dict):
        d = d.copy()
        for k, v in d.iteritems():
            d[k] = recursive_dictcopy(v)
    return d


class Shot:
    def __init__(self, comment="", **kw):
        self.comment = comment

        self.descr = kw.pop('descr', '')
        self.stars = kw.pop('stars', '')

        self.attrs = ('expt', 'shn', 'dig', 'head', 'amp_default', 'lines', 'times')
        for attr in self.attrs:
            setattr(self, attr, kw.pop(attr, None))

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

    def all_keys(self):
        keys = [self.head.keys['R']]
        for tip in self.head.tips:
            keys.extend([tip.keys['V'], tip.keys['I']])
        return keys

    def unique_keys(self):
        unique_keys = np.unique(self.all_keys())
        if unique_keys[0] is None:
            unique_keys = unique_keys[1:]
        return list(unique_keys)

    def __repr__(self):
        return "%d: %s" % (self.shn, self.comment)

    def mapsig(self, x, line):
        self.unique_sigs = {k: x[self.get(line, 'mapping', k)].astype('d') 
                for k in self.unique_keys()}

        t = x['t']
        R = self.unique_sigs[self.head.keys['R']]
        S = dict(R=PositionSignal(R, t, name='R'))

        for tip in self.head.tips:
            i, keyV, keyI = tip.number, tip.keys['V'], tip.keys['I']
            if keyV is None:
                V = None
            else:
                V = VoltageSignal(self.unique_sigs[keyV], t, name='V%d' % i)
            if keyI is None:
                S['V%d' % i] = V
            else:
                S['I%d' % i] = CurrentSignal(self.unique_sigs[keyI], t, V=V, name='I%d' % i)

        return S

    def calib(self, S, line):
        for k in self.unique_sigs.iterkeys():
            amp = self.get(line, 'amp', k)
            amp.apply(self.unique_sigs[k])


class ShotContainer:
    def __init__(self):
        self.x = OrderedDict()

    def __getitem__(self, indx):
        return self.x[indx]

    @property
    def times(self):
        return {k: v.times for k, v in self.x.iteritems()}

    @property
    def descr(self):
        return {k: v.descr for k, v in self.x.iteritems()}

    @property
    def stars(self):
        return {k: v.stars for k, v in self.x.iteritems()}


class Experiment(ShotContainer):
    def __init__(self, date=None, campaign=None, ShotClass=Shot):
        self.date, self.campaign, self.ShotClass = date, campaign, ShotClass
        ShotContainer.__init__(self)
    
    def add(self, shn, *args, **kw):
        self.x[shn] = self.ShotClass(*args, expt=self, shn=shn, **kw)
        
    def rep(self, shn, shn0, comment="", **kw):
        self.x[shn] = self.x[shn0].copy(comment, shn=shn, **kw)

    @property
    def shots(self):
        return np.array(self.x.keys())

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

    @property
    def shots(self):
        return np.concatenate([E.shots for E in self.x.itervalues()])

    def find_shot(self, shn):
        for E in self.x.itervalues():
            try: 
                return E[shn]
            except KeyError:
                pass
    
    def __repr__(self):
        s = [str(v) + "\n" for v in self.x.itervalues()]
        return np.array(s).tostring()



