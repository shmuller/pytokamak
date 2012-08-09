import numpy as np

from collections import OrderedDict

import probe


def_window = slice(2048, None)
def_mapping = dict(R='VOL3', V='VOL1', I1='CUR1', I2='CUR2', VF='VOL2')


Amp = probe.Amp
ampUnity = Amp(fact=1.)
ampInv   = Amp(fact=-1.)

fixpoints = (-1.8767, -106), (3.8011, 336)
ampR  = Amp(fixpoints=fixpoints)
ampV  = Amp(fact=100., offs=-183.76)
ampVF = Amp(fact=100.)

Preamp1 = {
    2: Amp(fact=1.032).inv(), 
    5: Amp(fact=2.580).inv()}

Preamp2 = {
    2: Amp(fact=1.936).inv(), 
    5: Amp(fact=4.840).inv()}

CurrentProbe1 = {
    10: Amp(fact=0.5*10/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10)}

CurrentProbe2 = {
    10: Amp(fact=0.5*10/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10)}

    
class Shot:
    def __init__(self, comment="", expt=None, shn=None,
            window=def_window, mapping=def_mapping, amp=None,
            ampR  = ampR, 
            ampV  = ampV, 
            ampI1 = CurrentProbe1[20],
            ampI2 = CurrentProbe2[20],
            ampVF = ampVF):

        if amp is None:
            amp = dict(R=ampR, V=ampV, I1=ampI1, I2=ampI2, VF=ampVF)

        self.comment = comment
        self.expt = expt
        self.shn = shn
        self.window = window
        self.mapping = mapping
        self.amp = amp

    def copy(self, comment="", expt=None, shn=None):
        if expt is None: 
            expt = self.expt
        if shn is None:
            shn = self.shn
        return Shot(comment=comment, expt=expt, shn=shn,
            window=self.window, mapping=self.mapping, amp=self.amp)

    def __repr__(self):
        return "%d: %s" % (self.shn, self.comment)


class Experiment:
    def __init__(self, date=None):
        self.date = date
        self.x = OrderedDict()

    def __getitem__(self, indx):
        return self.x[indx]

    def add(self, shn, *args, **kw):
        self.x[shn] = Shot(*args, expt=self, shn=shn, **kw)
        
    def rep(self, shn, shn0, comment=""):
        self.x[shn] = self.x[shn0].copy(comment, shn=shn)

    @property
    def shots(self):
        return np.array(self.x.keys())

    def __repr__(self):
        s = ["  " + str(v) + "\n" for v in self.x.itervalues()]
        return self.date + ":\n" + np.array(s).tostring()


class Campaign:
    def __init__(self):
        self.x = OrderedDict()
        
    def __getitem__(self, indx):
        return self.x[indx]

    def add_experiment(self, *args, **kw):
        E = Experiment(*args, **kw)
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


campaign = Campaign()

############################################
E = campaign.add_experiment(date="20120405")

E.add(27684, "", 
             ampI1 = ampInv*CurrentProbe1[10]*Preamp1[5],
             ampI2 = ampInv*CurrentProbe2[10]*Preamp2[5])

E.rep(27685, 27684)
E.rep(27686, 27684, "Both current probes on tip 1")

E.add(27687, "Changed direction of 2nd current probe",
             ampI1 = ampInv*CurrentProbe1[10]*Preamp1[5],
             ampI2 = CurrentProbe2[10]*Preamp2[5])

E.rep(27688, 27687)
E.rep(27689, 27687, "First shot of Leena's experiment")
E.rep(27690, 27687)

E.add(27691, "Current measurement from 10 mA/div to 20 mA/div",
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
             ampI2 = CurrentProbe2[20]*Preamp2[5])

E.rep(27692, 27691, "All the way through, but signals saturated")

E.add(27693, "Current probes 50 mA/div",
             ampI1 = ampInv*CurrentProbe1[50]*Preamp1[5],
             ampI2 = CurrentProbe2[50]*Preamp2[5])

E.rep(27694, 27693, "HWM resumed experiment")
E.rep(27695, 27693, "Calibration after this shot")


############################################
E = campaign.add_experiment(date="20120621")

# Henrik Mayer
E.add(28232, "304 mm, 1.5 s",
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
             ampI2 = CurrentProbe2[20]*Preamp2[5])

E.rep(28233, 28232, "440 mm, 5.0 s, disrupted before")
E.rep(28234, 28232, "204 mm, 4.8 s, disrupted before")


############################################
E = campaign.add_experiment(date="20120622")

# Rachael McDermott
E.add(28239, "304 mm, 3.8 s, 20 mA/div", 
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
             ampI2 = CurrentProbe2[20]*Preamp2[5])

E.add(28240, "440 mm, 3.8 s, 50 mA/div",
             ampI1 = ampInv*CurrentProbe1[50]*Preamp1[5],
             ampI2 = CurrentProbe2[50]*Preamp2[5])

E.rep(28241, 28239, "440 mm, 3.8 s, 20 mA/div")
E.rep(28242, 28239, "440 mm, 3.2 s, 20 mA/div")
E.rep(28243, 28239, "440 mm, 3.2 s, 20 mA/div -> saturated I2")

# Tim Happel
E.rep(28244, 28243, "204 mm, 1.5 s, 20 mA/div -> 440 mm acc.")

E.add(28245, "204 mm, 1.0 s, 20 mA/div, preamps 2x (was 5x)",
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
             ampI2 = CurrentProbe2[20]*Preamp2[2])

E.rep(28246, 28245, "254 mm, 0.85 s, 20 mA/div, preamps 2x -> L-H transition!")

# Henrik Mayer
E.rep(28250, 28245, "304 mm, 4.0 s, 20 mA/div, preamps 2x")
E.rep(28251, 28245, "304 mm, 1.8 s, 20 mA/div, preamps 2x, (bias voltage less positive)")
E.rep(28252, 28245, "354 mm, 1.75 s, 20 mA/div, preamps 2x -> 1.8 s acc.")
E.rep(28253, 28245, "304 mm, 1.75 s, 20 mA/div, preamps 2x (repeat 28250)")
E.rep(28254, 28245, "304 mm, 1.75 s, 20 mA/div, preamps 2x (repeat 28251)")


"""
if shn < 28426:
    ampI1 *= ampInv

if shn < 27687:
    ampI2 *= ampInv
"""

