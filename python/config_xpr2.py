import numpy as np

from config import Experiment, Campaign

from probe import Amp, PositionSignal, VoltageSignal, CurrentSignal


ampUnity = Amp(fact=1.)
ampInv   = Amp(fact=-1.)

fixpoints = (-1.8767, -106), (3.8011, 336)
ampR_LPS  = Amp(fixpoints=fixpoints)
ampV_LPS  = Amp(fact=100., offs=-183.76)

fixpoints = (3.6812, -72), (7.0382, 170)
ampR_XPR = Amp(fixpoints=fixpoints)
ampV_XPR = Amp(fact=100., offs=-69.2227)

ampVF = Amp(fact=100.)


Preamp1 = {
     2: Amp(fact=1.032).inv(), 
     5: Amp(fact=2.580).inv()}

Preamp2 = {
     2: Amp(fact=1.936).inv(), 
     5: Amp(fact=4.840).inv()}

CurrentProbe1 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10), 
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10)}

CurrentProbe2 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10),
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10)}

CurrentProbe3 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10),
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10)}


class Tip:
    def __init__(self, area, proj_area, pos, V_links=None, I_links=None):
        self.area = area
        self.proj_area = proj_area
        self.pos = pos
        self.links = dict(V=V_links, I=I_links)

    def get_V(self, line, field):
        return None


class CylindricalTip(Tip):
    def __init__(self, r, z, *args, **kw):
        self.r, self.z = r, z
        area = 2.*np.pi*r*z
        proj_area = 2*r*z
        Tip.__init__(self, area, proj_area, *args, **kw)


class Head:
    def __init__(self, tips, R_links=None):
        self.tips = tips
        self.links = dict(R=R_links)


class TipXPR(CylindricalTip):
    def __init__(self, *args, **kw):
        CylindricalTip.__init__(self, 0.0005, 0.003, *args, **kw)


tip1 = TipXPR(pos='lower left', V_links='ampV', I_links='ampI1')
tip2 = TipXPR(pos='lower right', V_links='ampV', I_links='ampI2')
tip3 = TipXPR(pos='upper', V_links='ampVF', I_links=None)

head = Head(tips=(tip1, tip2, tip3), R_links='ampR')


amp_default = dict(
            ampR  = None, 
            ampV  = None, 
            ampI1 = CurrentProbe1[20],
            ampI2 = CurrentProbe2[20],
            ampI3 = CurrentProbe3[20],
            ampVF = ampVF)

amp_XPR = dict(
            ampR = ampR_XPR, 
            ampV = ampV_XPR)

amp_LPS = dict(
            ampR = ampR_LPS, 
            ampV = ampV_LPS)

amp_mapping_XPR = dict(
            ampR  = 'S5', 
            ampV  = 'S1', 
            ampI1 = 'S4',
            ampI2 = 'S2',
            ampI3 = 'S6',
            ampVF = 'S6')

amp_mapping_LPS = dict(
            ampR  = 'VOL3', 
            ampV  = 'VOL1',
            ampI1 = 'CUR1',
            ampI2 = 'CUR2',
            ampI3 = 'VOL2',
            ampVF = 'VOL2')

lines_XPR = dict(amp=amp_XPR, mapping=amp_mapping_XPR)
lines_LPS = dict(amp=amp_LPS, mapping=amp_mapping_LPS)

lines = dict(XPR=lines_XPR, LPS=lines_LPS)

def_LPS = dict(dig='LPS', head=head, amp_default=amp_default, lines=dict(LPS=lines_LPS))

class Shot2:
    def __init__(self, comment="", expt=None, shn=None, dig=None,
            head=None, amp_default=None, lines=None, **kw):
        self.comment = comment
        self.expt = expt
        self.shn = shn
        self.dig = dig
        self.head = head

        self.amp_default = amp_default.copy()
        for k in self.amp_default.keys():
            try:
                self.amp_default[k] = kw.pop(k)
            except KeyError:
                pass

        self.lines = lines
    
    def copy(self, comment="", expt=None, shn=None):
        if expt is None: 
            expt = self.expt
        if shn is None:
            shn = self.shn
        return self.__class__(comment, expt, shn, self.dig, 
                              self.head, self.amp_default, self.lines)

    def get(self, line, what, key):
        if key is None:
            return None
        elif what == 'amp':
            return self.lines[line]['amp'].get(key, self.amp_default[key])
        else:
            return self.lines[line][what][key]

    def all_keys(self):
        keys = [self.head.links['R']]
        for tip in self.head.tips:
            keys.extend([tip.links['V'], tip.links['I']])
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
        R = self.unique_sigs[self.head.links['R']]
        S = dict(R=PositionSignal(R, t, name='R'))

        for i, tip in enumerate(self.head.tips, start=1):
            keyV, keyI = tip.links['V'], tip.links['I']
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


class ExperimentAUG(Experiment):
    def __init__(self, *args, **kw):
        Experiment.__init__(self, *args, ShotClass=Shot2, **kw)


class CampaignAUG(Campaign):
    def __init__(self, *args, **kw):
        Campaign.__init__(self, *args, ExperimentClass=ExperimentAUG, **kw)


campaign = CampaignAUG()


############################################
E = campaign.add_experiment(date="20120405")

E.add(27684, "",
             ampI1 = ampInv*CurrentProbe1[10]*Preamp1[5],
             ampI2 = ampInv*CurrentProbe2[10]*Preamp2[5], **def_LPS)

E.rep(27685, 27684)
E.rep(27686, 27684, "Both current probes on tip 1")

E.add(27687, "Changed direction of 2nd current probe",
             ampI1 = ampInv*CurrentProbe1[10]*Preamp1[5],
             ampI2 = CurrentProbe2[10]*Preamp2[5], **def_LPS)

E.rep(27688, 27687)
E.rep(27689, 27687, "First shot of Leena's experiment")
E.rep(27690, 27687)

E.add(27691, "Current measurement from 10 mA/div to 20 mA/div",
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
             ampI2 = CurrentProbe2[20]*Preamp2[5], **def_LPS)

E.rep(27692, 27691, "All the way through, but signals saturated")

E.add(27693, "Current probes 50 mA/div",
             ampI1 = ampInv*CurrentProbe1[50]*Preamp1[5],
             ampI2 = CurrentProbe2[50]*Preamp2[5], **def_LPS)

E.rep(27694, 27693, "HWM resumed experiment")
E.rep(27695, 27693, "Calibration after this shot")



