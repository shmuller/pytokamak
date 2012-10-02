import numpy as np

from collections import OrderedDict

from probe import PositionSignal, VoltageSignal, CurrentSignal

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


class Experiment:
    def __init__(self, date=None, campaign=None, ShotClass=Shot):
        self.date, self.campaign, self.ShotClass = date, campaign, ShotClass
        self.x = OrderedDict()

    def __getitem__(self, indx):
        return self.x[indx]

    def add(self, shn, *args, **kw):
        self.x[shn] = self.ShotClass(*args, expt=self, shn=shn, **kw)
        
    def rep(self, shn, shn0, comment=""):
        self.x[shn] = self.x[shn0].copy(comment, shn=shn)

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



