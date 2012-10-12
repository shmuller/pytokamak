import numpy as np

from config import CylindricalTip, Head, Campaign

from probe import Amp


ampUnity = Amp(fact=1.)
ampInv   = Amp(fact=-1.)

fixpoints = (-1.8767, -0.106), (3.8011, 0.336)
ampR_LPS  = Amp(fixpoints=fixpoints)
ampV_LPS  = Amp(fact=100., offs=-183.76)

fixpoints = (3.6812, -0.072), (7.0382, 0.170)
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
    50: Amp(fact=0.5*50/10),
  5000: Amp(fact=0.5*5000/10)}

CurrentProbe2 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10),
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10),
  5000: Amp(fact=0.5*5000/10)}

CurrentProbe3 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10),
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10),
  5000: Amp(fact=0.5*5000/10)}


class TipXPR(CylindricalTip):
    def __init__(self, *args, **kw):
        CylindricalTip.__init__(self, 0.0005, 0.003, *args, **kw)


tip1 = TipXPR(number=1, pos='lower left', V_keys='ampV', I_keys='ampI1')
tip2 = TipXPR(number=2, pos='lower right', V_keys='ampV', I_keys='ampI2')
tip3 = TipXPR(number=3, pos='upper', V_keys='ampVF', I_keys=None)

head = Head(tips=(tip1, tip2, tip3), R_keys='ampR')

tip3I = TipXPR(number=3, pos='upper', V_keys='ampV', I_keys='ampI3')
headI = Head(tips=(tip1, tip2, tip3I), R_keys='ampR')

amp_default_unity = dict(
            ampR  = ampUnity, 
            ampV  = ampUnity, 
            ampI1 = ampUnity,
            ampI2 = ampUnity,
            ampI3 = ampUnity,
            ampVF = ampUnity)

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

mapping_XPR = dict(
            ampR  = 'S5', 
            ampV  = 'S1', 
            ampI1 = 'S4',
            ampI2 = 'S2',
            ampI3 = 'S6',
            ampVF = 'S6')

mapping_LPS = dict(
            ampR  = 'VOL3', 
            ampV  = 'VOL1',
            ampI1 = 'CUR1',
            ampI2 = 'CUR2',
            ampI3 = 'VOL2',
            ampVF = 'VOL2')

lines_XPR = dict(amp=amp_XPR, mapping=mapping_XPR)
lines_LPS = dict(amp=amp_LPS, mapping=mapping_LPS)

lines = dict(XPR=lines_XPR, LPS=lines_LPS)

def_LPS = dict(dig='LPS', head=head, amp_default=amp_default, lines=dict(LPS=lines_LPS))
def_XPR = dict(dig='XPR', head=head, amp_default=amp_default, lines=dict(XPR=lines_XPR))

def_XPR_LPS = dict(dig='XPR', head=head, amp_default=amp_default, lines=lines)

allI_XPR = dict(dig='XPR', head=headI, amp_default=amp_default, lines=dict(XPR=lines_XPR))

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
             ampI2 = CurrentProbe2[50]*Preamp2[5],
             stars = '*****')

E.rep(27694, 27693, "HWM resumed experiment", 
             stars='*****')
E.rep(27695, 27693, "Calibration after this shot", 
             stars='*****')


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

E.rep(28434, 28429, "20 mA/div, 0.1 kHz, all 3 tips on bias voltage", 
             head = headI,
             ampI1 = CurrentProbe1[20],
             ampI2 = CurrentProbe2[20], 
             ampI3 = CurrentProbe3[20])

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
             dig='XPR', head=head, amp_default=amp_default_unity, 
             lines=dict(XPR=dict(amp={}, mapping=mapping_XPR), 
                        LPS=dict(amp={}, mapping=mapping_LPS)),
             ampI1 = CurrentProbe1[20],
             ampI2 = CurrentProbe2[20])

E.rep(28507, 28504, "Calibration, 10 Vpp into 50 Ohm (+/-0.1 A)")
E.rep(28508, 28504, "Signal also on bias voltage")


############################################
E = campaign.add_experiment(date="20121011")

E.add(28633, "DAQ test", **def_XPR)

E.add(28634, "Sweep attached to 2x100 V Kepco pair, all tips on sweep",
             times = (950, 1850),
             ampI1 = CurrentProbe1[5000],
             ampI2 = CurrentProbe2[5000],
             ampI3 = CurrentProbe3[5000], descr = """
             Fuse blown on whole Kepco rack. No data""", **allI_XPR)

E.rep(28636, 28634, "Switch Kepco's off", descr = """
             No motion. No signals?""")

E.rep(28637, 28636, "Acquire trigger signals", descr = """
             No motion. No signals?""")

E.rep(28641, 28637, "TTL via LWL 1061 on channel 5",
             descr = "Nothing came through")

E.rep(28643, 28641, "Sine via fcn gen on channel 5",
             descr = "")

E.rep(28645, 28643, "Kepcos on separate trafo, sweep on",
             descr = "")

E.add(28646, "Change sensitity",
             times = (0.950, 1.850),
             ampI1 = CurrentProbe1[20],
             ampI2 = CurrentProbe2[20],
             ampI3 = CurrentProbe3[20], descr = """
             """, **allI_XPR)

E.rep(28647, 28646, "Repeat")

E.rep(28648, 28647, "Reset local timer for PPG TS06", 
        descr = """
        Worked, but tip3 apparently short circuits
        """)

E.rep(28649, 28648, "Go to three plunges behind the wall", 
        times = (1.0, 2.0, 3.0),
        descr = """Short circuit on tip3 on plunge 0, then
        UCSD Kepco trips, then current 3 follows voltage.
        """)

E.rep(28650, 28649, "Take tip 3 off bias voltage, 3rd plunge slower")

E.rep(28651, 28650, "Position signal on S8", 
        XPR_mapping_ampR='S8')



