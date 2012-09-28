import numpy as np

from config import Map, Shot, Experiment, Campaign

from config import CylindricalTip, Head, Shot2

import probe

DictView = probe.DictView

Amp = probe.Amp
ampUnity = Amp(fact=1.)
ampInv   = Amp(fact=-1.)

# LPS defaults
mapping_LPS = Map(R='VOL3', V=['VOL1', 'VOL1', 'VOL2'], I=['CUR1', 'CUR2', None])
mapping_LPS_allI = Map(R='VOL3', V=['VOL1', 'VOL1', 'VOL1'], I=['CUR1', 'CUR2', 'VOL2'])

fixpoints = (-1.8767, -106), (3.8011, 336)
ampR_LPS  = Amp(fixpoints=fixpoints)
ampV_LPS  = Amp(fact=100., offs=-183.76)

def_LPS = dict(dig='LPS', mapping=mapping_LPS, ampR=ampR_LPS, ampV=ampV_LPS)
allI_LPS = dict(dig='LPS', mapping=mapping_LPS_allI, ampR=ampR_LPS, ampV=ampV_LPS)

# XPR defaults
mapping_XPR = Map(R='S5', V=['S1', 'S1', 'S6'], I=['S4', 'S2', None])
mapping_XPR_allI = Map(R='S5', V=['S1', 'S1', 'S1'], I=['S4', 'S2', 'S6'])

fixpoints = (3.6812, -72), (7.0382, 170)
ampR_XPR = Amp(fixpoints=fixpoints)
ampV_XPR = Amp(fact=100., offs=-69.2227)

def_XPR = dict(dig='XPR', mapping=mapping_XPR, ampR=ampR_XPR, ampV=ampV_XPR)
allI_XPR = dict(dig='XPR', mapping=mapping_XPR_allI, ampR=ampR_XPR, ampV=ampV_XPR)

# LPS and XPR simultaneously
def_XPR_LPS = dict(alt_dig=def_LPS)
def_XPR_LPS.update(def_XPR)

allI_XPR_LPS = dict(alt_dig=allI_LPS)
allI_XPR_LPS.update(allI_XPR)


ampVF = Amp(fact=100.)

def_amp  = Map(R='ampR', V=['ampV', None, ampVF], I=['ampI1', 'ampI2', None])
allI_amp = Map(R='ampR', V=['ampV', None, None] , I=['ampI1', 'ampI2', 'ampI3'])

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


tip1 = TipXPR(pos='lower left', V_links='ampV', I_links='ampI1')
tip2 = TipXPR(pos='lower right', V_links='ampV', I_links='ampI2')
tip3 = TipXPR(pos='upper', V_links='ampVF', I_links=None)

head = Head(tips=(tip1, tip2, tip3), R_links='ampR')


amp_settings = dict(
            ampR  = None, 
            ampV  = None, 
            ampI1 = CurrentProbe1[20],
            ampI2 = CurrentProbe2[20],
            ampI3 = CurrentProbe3[20],
            ampVF = ampVF)

amp_mapping_XPR = dict(
            ampR  = 'S5', 
            ampV  = 'S1', 
            ampI1 = 'S4',
            ampI2 = 'S2',
            ampI3 = None,
            ampVF = 'S6')

amp_mapping_LPS = dict(
            ampR  = 'VOL3', 
            ampV  = 'VOL1',
            ampI1 = 'CUR1',
            ampI2 = 'CUR2',
            ampI3 =  None ,
            ampVF = 'VOL2')

amp_mappings = dict(
            XPR = amp_mapping_XPR,
            LPS = amp_mapping_LPS)


shot2 = Shot2(head, amp_settings, ampR=ampR_XPR, amp_mappings=amp_mappings)


class ShotAUG(Shot):
    def __init__(self, comment="", expt=None, shn=None, dig=None, 
            mapping=None, amp_template=def_amp, amp=None, alt_dig=None,
            ampR  = None, 
            ampV  = None, 
            ampI1 = CurrentProbe1[20],
            ampI2 = CurrentProbe2[20],
            ampI3 = CurrentProbe3[20],
            ampVF = ampVF):

        self.amp_template = amp_template

        if amp is None:
            amp = self.amp_template.copy()
            if amp.R == 'ampR':
                amp.R = ampR
            for key, val in (('ampV', ampV), ('ampVF', ampVF)):
                amp.replace('V', key, val)
            for key, val in (('ampI1', ampI1), ('ampI2', ampI2), ('ampI3', ampI3)):
                amp.replace('I', key, val)

        #if amp is None:
        #    amp = Map(R=ampR, V=[ampV, None, ampVF], I=[ampI1, ampI2, None])

        Shot.__init__(self, comment, expt, shn, dig, mapping, amp, alt_dig)

    def add_dig(self, dig=None, mapping=None, amp=None, ampR=None, ampV=None):
        if amp is None:
            amp = self.amp[self.dig].copy()
            amp.R = ampR
            amp.V[self.amp_template.V == 'ampV'] = ampV

            #a = self.amp[self.dig]
            #amp = Map(R=ampR, V=[ampV, None, a.V[2]], I=a.I)

        Shot.add_dig(self, dig, mapping, amp)

    def copy(self, comment="", expt=None, shn=None):
        if expt is None: 
            expt = self.expt
        if shn is None:
            shn = self.shn
        return self.__class__(comment=comment, expt=expt, shn=shn, dig=self.dig, 
                              mapping=self.mapping, amp_template=self.amp_template, 
                              amp=self.amp)


class ExperimentAUG(Experiment):
    def __init__(self, *args, **kw):
        Experiment.__init__(self, *args, ShotClass=ShotAUG, **kw)


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


############################################
E = campaign.add_experiment(date="20120621")

# Henrik Mayer
E.add(28232, "304 mm, 1.5 s",
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
             ampI2 = CurrentProbe2[20]*Preamp2[5], **def_LPS)

E.rep(28233, 28232, "440 mm, 5.0 s, disrupted before")
E.rep(28234, 28232, "204 mm, 4.8 s, disrupted before")


############################################
E = campaign.add_experiment(date="20120622")

# Rachael McDermott
E.add(28239, "304 mm, 3.8 s, 20 mA/div", 
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
             ampI2 = CurrentProbe2[20]*Preamp2[5], **def_LPS)

E.add(28240, "440 mm, 3.8 s, 50 mA/div",
             ampI1 = ampInv*CurrentProbe1[50]*Preamp1[5],
             ampI2 = CurrentProbe2[50]*Preamp2[5], **def_LPS)

E.rep(28241, 28239, "440 mm, 3.8 s, 20 mA/div")
E.rep(28242, 28239, "440 mm, 3.2 s, 20 mA/div")
E.rep(28243, 28239, "440 mm, 3.2 s, 20 mA/div -> saturated I2")

# Tim Happel
E.rep(28244, 28243, "204 mm, 1.5 s, 20 mA/div -> 440 mm acc.")

E.add(28245, "204 mm, 1.0 s, 20 mA/div, preamps 2x (was 5x)",
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
             ampI2 = CurrentProbe2[20]*Preamp2[2], **def_LPS)

E.rep(28246, 28245, "254 mm, 0.85 s, 20 mA/div, preamps 2x -> L-H transition!")

# Henrik Mayer
E.rep(28250, 28245, "304 mm, 4.0 s, 20 mA/div, preamps 2x")
E.rep(28251, 28245, "304 mm, 1.8 s, 20 mA/div, preamps 2x, (bias voltage less positive)")
E.rep(28252, 28245, "354 mm, 1.75 s, 20 mA/div, preamps 2x -> 1.8 s acc.")
E.rep(28253, 28245, "304 mm, 1.75 s, 20 mA/div, preamps 2x (repeat 28250)")
E.rep(28254, 28245, "304 mm, 1.75 s, 20 mA/div, preamps 2x (repeat 28251)")


############################################
E = campaign.add_experiment(date="20120712")

E.add(28379, "Fixed probe @2564.05 mm", 
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
             ampI2 = CurrentProbe2[20]*Preamp2[2], **def_XPR_LPS)

E.rep(28380, 28379, "Fixed probe @2569.05 mm -> no data")
E.rep(28381, 28379, "-200 V bias -> Kepco breaks in at 0.5 s")
E.rep(28382, 28379, "@2569.10 mm, sweep -> data all the way to 6 s")
E.rep(28383, 28379, "@2589.08 mm, DC offset with small sweep -> worked")


############################################
E = campaign.add_experiment(date="20120713")

E.add(28389, "Fixed probe @2569.05 mm, -180 V with sweep", 
             ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
             ampI2 = CurrentProbe2[20]*Preamp2[2], **def_XPR_LPS)

E.rep(28390, 28389, "-80 V / 150 V sweep at 100 Hz")
E.rep(28394, 28389)


E.add(28395, "Turn 2nd current probe", 
             ampI1 = CurrentProbe1[20]*Preamp1[2],
             ampI2 = CurrentProbe2[20]*Preamp2[2], **def_XPR_LPS)

E.rep(28403, 28395)
E.rep(28404, 28395)
E.rep(28405, 28395)
E.rep(28406, 28395)
E.rep(28407, 28395)


############################################
E = campaign.add_experiment(date="20120717")

E.add(28419, "Fcn gen. 20 Vpp, +8 VDC, 0.5 kHz (saturates in Isat regime)", 
             ampI1 = CurrentProbe1[5],
             ampI2 = CurrentProbe2[5], **def_XPR_LPS)

E.rep(28420, 28419, "Fcn gen. 12 Vpp, 4 VDC")
E.rep(28421, 28419)
E.rep(28422, 28419)
E.rep(28423, 28419)
E.rep(28424, 28419, "Asym. waveform, 7 Vpp, 1 VDC, 200 Hz")

E.add(28425, "6.9 Vpp, 1 VDC, 100 Hz, 1 mA/div", 
             ampI1 = CurrentProbe1[1],
             ampI2 = CurrentProbe2[1], **def_XPR_LPS)

E.rep(28426, 28425, "First data! Kepco breaks in in Isat")
E.rep(28427, 28425, "Change sweep pars -> Kepco still breaks")
E.rep(28428, 28425, "Back to other fcn gen. -> saturated at 1 mA/div")

E.add(28429, "Back to 5 mA/div -> no data", 
             ampI1 = CurrentProbe1[5],
             ampI2 = CurrentProbe2[5], **def_XPR_LPS)

E.add(28434, "20 mA/div, 0.1 kHz, all 3 tips on bias voltage", 
             amp_template = allI_amp,
             ampI1 = CurrentProbe1[20],
             ampI2 = CurrentProbe2[20], 
             ampI3 = CurrentProbe3[20], **allI_XPR_LPS)

E.rep(28435, 28434, "0.5 kHz, plunge at 1 s")
E.rep(28436, 28434)


############################################
E = campaign.add_experiment(date="20120719")

E.add(28442, "0.5 kHz, 3rd pin VF", 
             ampI1 = CurrentProbe1[20],
             ampI2 = CurrentProbe2[20], **def_XPR_LPS)

E.rep(28444, 28442, "Max plunge at 1.75 s and 3.95 s")
E.rep(28445, 28442, "100 V Kepco")
E.rep(28446, 28442, "200 V Kepco, 12 Vpp, 5 VDC, 0.5 kHz")

# Francois Ryter
E.rep(28448, 28442, "1 kHz, 16 Vpp (digital fcn gen.), VDC from Kepco")
E.rep(28449, 28442, "0.5 kHz, reduced VDC slightly")
E.rep(28450, 28442, "2nd plunge to 1.6 s")
E.rep(28451, 28442, "Max penetration -> shot didn't run")
E.rep(28452, 28442, "Max penetration -> arcs")


############################################
E = campaign.add_experiment(date="20120720")

E.add(28455, "Acquisition with turned-off Kepco", 
             ampI1 = CurrentProbe1[20],
             ampI2 = CurrentProbe2[20], **def_XPR_LPS)

E.rep(28466, 28455, "0.5 kHz, 16 Vpp, Kepco offset just avoids saturation")
E.rep(28467, 28455)
E.rep(28468, 28455)

E.rep(28469, 28455, "H: 14 Vpp, max offset on Kepco")
E.rep(28472, 28455, "He: 3 plunges, max penetration")
E.rep(28473, 28455, "He again: 1 plunges at 150 mm")


############################################
E = campaign.add_experiment(date="20120726")

E.add(28504, "Calibration, no signals attached", 
             dig = 'XPR', mapping = mapping_XPR,
             ampR = ampUnity, ampV = ampUnity, ampVF = ampUnity,
             alt_dig = dict(dig='LPS', mapping=mapping_LPS, ampR=ampUnity, ampV=ampUnity),
             ampI1 = CurrentProbe1[20],
             ampI2 = CurrentProbe2[20])

E.rep(28507, 28504, "Calibration, 10 Vpp into 50 Ohm (+/-0.1 A)")
E.rep(28508, 28504, "Signal also on bias voltage")




