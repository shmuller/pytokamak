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


class Map:
    def __init__(self, R, VI=None, V=None, I=None):
        self.R = R
        if VI is None:
            self.VI = np.array(zip(V, I), dtype=object)
        else:
            self.VI = np.asarray(VI, dtype=object)
        self.V, self.I = self.VI[:,0], self.VI[:,1]

    def __repr__(self):
        return repr(self.VI)

    def copy(self):
        return Map(self.R, self.VI.copy())

    def unique(self, attr):
        a = np.unique(getattr(self, attr))
        if a[0] is None:
            return a[1:]
        else:
            return a

    def replace(self, attr, key, value):
        a = getattr(self, attr)
        a[a == key] = value


class Shot:
    def __init__(self, comment="", expt=None, shn=None, dig=None, 
                 mapping=None, amp=None, alt_dig=None):

        if not isinstance(mapping, dict):
            mapping = {dig: mapping}
        if not isinstance(amp, dict):
            amp = {dig: amp}

        self.comment = comment
        self.expt = expt
        self.shn = shn
        self.dig = dig
        self.mapping = mapping
        self.amp = amp

        if alt_dig is not None:
            self.add_dig(**alt_dig)

    def add_dig(self, dig, mapping=None, amp=None):
        self.mapping[dig] = mapping
        self.amp[dig] = amp

    def copy(self, comment="", expt=None, shn=None):
        if expt is None: 
            expt = self.expt
        if shn is None:
            shn = self.shn
        return self.__class__(comment=comment, expt=expt, shn=shn, dig=self.dig, 
                              mapping=self.mapping, amp=self.amp)

    def __repr__(self):
        return "%d: %s" % (self.shn, self.comment)

    def mapsig(self, x, dig):
        mapping = self.mapping[dig]

        t = x['t']
        R = x[mapping.R].astype('d')
        S = dict(R=PositionSignal(R, t, name='R'))

        unique_sigs = {k: x[k].astype('d') for k in mapping.unique('VI')}

        for i, (mapV, mapI) in enumerate(mapping.VI, start=1):
            if mapV is None:
                V = None
            else:
                V = VoltageSignal(unique_sigs[mapV], t, name='V%d' % i)
            if mapI is None:
                S[i] = V
            else:
                S[i] = CurrentSignal(unique_sigs[mapI], t, V, name='I%d' % i)

        return S

    def calib(self, S, dig):
        amp = self.amp[dig]
        S['R'] *= amp.R
        for i, (ampV, ampI) in enumerate(amp.VI, start=1):
            if isinstance(S[i], CurrentSignal):
                S[i] *= ampI
                V = S[i].V
            else:
                V = S[i]
            if V is not None and ampV is not None:
                V *= ampV


class Shot2:
    def __init__(self, comment="", **kw):
        self.comment = comment

        self.attrs = ('expt', 'shn', 'dig', 'head', 'amp_default', 'lines')
        for attr in self.attrs:
            setattr(self, attr, kw.pop(attr, None))

        self.amp_default = self.amp_default.copy()
        for k in set(self.amp_default.keys()) & set(kw.keys()):
            self.amp_default[k] = kw.pop(k)
                
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
                S[i] = V
            else:
                S[i] = CurrentSignal(self.unique_sigs[keyI], t, V, name='I%d' % i)

        return S

    def calib(self, S, line):
        for k in self.unique_sigs.iterkeys():
            amp = self.get(line, 'amp', k)
            amp.apply(self.unique_sigs[k])


class Experiment:
    def __init__(self, date=None, campaign=None, ShotClass=Shot2):
        self.date, self.campaign, self.ShotClass = date, campaign, ShotClass
        self.x = OrderedDict()

    def __getitem__(self, indx):
        return self.x[indx]

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


class Campaign:
    def __init__(self, ExperimentClass=Experiment):
        self.ExperimentClass = ExperimentClass
        self.x = OrderedDict()
        
    def __getitem__(self, indx):
        return self.x[indx]

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



