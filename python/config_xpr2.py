import numpy as np

from config import CylindricalTip, Head, Campaign

from probe import Amp


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


class TipXPR(CylindricalTip):
    def __init__(self, *args, **kw):
        CylindricalTip.__init__(self, 0.0005, 0.003, *args, **kw)


tip1 = TipXPR(pos='lower left', V_keys='ampV', I_keys='ampI1')
tip2 = TipXPR(pos='lower right', V_keys='ampV', I_keys='ampI2')
tip3 = TipXPR(pos='upper', V_keys='ampVF', I_keys=None)

head = Head(tips=(tip1, tip2, tip3), R_keys='ampR')


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

campaign = Campaign()


############################################
E = campaign.add_experiment(date="20120405")

E.add(27684, "",
             ampI1 = ampInv*CurrentProbe1[10]*Preamp1[5],
             ampI2 = ampInv*CurrentProbe2[10]*Preamp2[5], **def_LPS)

E.rep(27685, 27684)
E.rep(27686, 27684, "Both current probes on tip 1")

E.rep(27687, 27686, "Changed direction of 2nd current probe",
             ampI2 = CurrentProbe2[10]*Preamp2[5])

E.rep(27688, 27687)
E.rep(27689, 27687, "First shot of Leena's experiment")
E.rep(27690, 27687)

E.rep(27691, 27690, "Current measurement from 10 mA/div to 20 mA/div",
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
             ampI2 = CurrentProbe2[20]*Preamp2[5])

E.rep(27692, 27691, "All the way through, but signals saturated")

E.rep(27693, 27692, "Current probes 50 mA/div",
             ampI1 = ampInv*CurrentProbe1[50]*Preamp1[5],
             ampI2 = CurrentProbe2[50]*Preamp2[5])

E.rep(27694, 27693, "HWM resumed experiment")
E.rep(27695, 27693, "Calibration after this shot")



