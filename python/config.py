import numpy as np

from collections import OrderedDict

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


class Experiment:
    def __init__(self, date=None, campaign=None, ShotClass=None):
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



