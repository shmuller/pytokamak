######################
# 2012 SPRING CAMPAIGN
######################

from config_xpr_common import *
from config import Campaign

campaign = Campaign()

tip1 = TipXPR(number=1, pos='lower left', V_keys='ampV', I_keys='ampI1')
tip2 = TipXPR(number=2, pos='lower right', V_keys='ampV', I_keys='ampI2')
tip3 = TipXPR(number=3, pos='upper', V_keys='ampVF', I_keys=None)

head = Head(tips=(tip1, tip2, tip3), R_keys='ampR')

tip3I = TipXPR(number=3, pos='upper', V_keys='ampV', I_keys='ampI3')
headI = Head(tips=(tip1, tip2, tip3I), R_keys='ampR')


############################################
E = campaign.add_experiment(date="20120405")

E.add(27684, "Plunge to 10 cm", 
        times = 1.2, 
        posit = 0.10,
        head = head,
        ampI1 = ampInv*CurrentProbe1[10]*Preamp1[5],
        ampI2 = ampInv*CurrentProbe2[10]*Preamp2[5], 
        descr = "Nice L-mode data", 
        stars = '**', **def_LPS)

E.rep(27685, 27684, "Repeat",
        descr = "Very similar to last shot",
        stars = '**')

E.rep(27686, 27685, "Technical shot: Both current probes on tip 1",
        posit = 0.13,
        descr = "Plunge in L-mode to 13 cm. Current measurements agree well.",
        stars = '')

E.rep(27687, 27686, "Changed direction of 2nd current probe",
        posit = 0.15,
        ampI2 = CurrentProbe2[10]*Preamp2[5],
        descr = "Plunge in L-mode to 15 cm. Arcs on way out.",
        stars = '***')

E.rep(27688, 27687, "2.5 cm at 3.9 s",
        times = 3.9,
        posit = 0.025,
        descr = "Plunge only to 2.5 cm.",
        stars = '*')

E.rep(27689, 27688, "First shot of Leena's experiment",
        posit = 0.10,
        descr = "Good L-mode measurements, but not far in.",
        stars = '**')

E.rep(27690, 27689, "Go to 20 cm",
        posit = 0.20,
        descr = "Made it to private flux region. Slightly saturated signals.",
        stars = '****')

E.rep(27691, 27690, "Current measurement from 10 mA/div to 20 mA/div",
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
        ampI2 = CurrentProbe2[20]*Preamp2[5],
        descr = 'Almost identical to last shot. Signals not saturated',
        stars = '****')

E.rep(27692, 27691, "All the way through, but signals saturated",
        posit = 0.34,
        descr = "Really nice measurements, but signals saturated on HFS",
        stars = "****")

E.rep(27693, 27692, "Current probes 50 mA/div, go to 20 cm",
        posit = 0.20,
        ampI1 = ampInv*CurrentProbe1[50]*Preamp1[5],
        ampI2 = CurrentProbe2[50]*Preamp2[5],
        descr = "Nice L-mode data.",
        stars = '*****')

E.rep(27694, 27693, "HWM resumed experiment", 
        times = 1.3,
        posit = 0.34,
        descr = """\
            Very nice L-mode data.
            Looks like the plunge is going through the core""",
        stars = '*****')

E.rep(27695, 27694, "Calibration after this shot", 
        times = 3.2,
        descr = """\
            Nice L-mode again. This time through private flux region.""",
        stars = '*****')


############################################
E = campaign.add_experiment(date="20120621")

# Hendrik Meyer
E.add(28232, "304 mm, 1.5 s",
        times = 1.65,
        posit = 0.20,
        head = head,
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
        ampI2 = CurrentProbe2[20]*Preamp2[5], 
        descr = """\
            Plunge before L-H transition. Nice L-mode data.
            Plasma shifted by 1.5 cm between inward/outward stroke.""",
        stars = "****", **def_LPS)

E.rep(28233, 28232, "440 mm, 5.0 s, disrupted before",
        times = 5.2,
        posit = 0.34,
        descr = "Disrupted before plunge.",
        stars = '')

E.rep(28234, 28233, "204 mm, 4.8 s, disrupted before",
        times = 4.9,
        posit = 0.10,
        descr = "Disrupted before plunge.",
        stars = '')


############################################
E = campaign.add_experiment(date="20120622")

# Rachael McDermott
E.add(28239, "304 mm, 3.8 s, 20 mA/div, I2 also on I3",
        times = 3.95,
        posit = 0.20,
        head = head,
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
        ampI2 = CurrentProbe2[20]*Preamp2[5],
        descr = "Very cold L-mode, very turbulent.",
        stars = '***', **def_LPS)

E.rep(28240, 28239, "440 mm, 3.8 s, 50 mA/div",
        times = 4.0,
        posit = 0.34,
        ampI1 = ampInv*CurrentProbe1[50]*Preamp1[5],
        ampI2 = CurrentProbe2[50]*Preamp2[5],
        descr = """\
            Plasma touches wall before plunge.
            Very low n and Te.""",
        stars = '**')

E.rep(28241, 28240, "440 mm, 3.8 s, 20 mA/div",
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
        ampI2 = CurrentProbe2[20]*Preamp2[5],
        descr = "Almost identical to last shot.",
        stars = '**')

E.rep(28242, 28241, "440 mm, 3.2 s, 20 mA/div",
        times = 3.4,
        descr = """\
            Very nice cold L-mode. Density LFS > HFS.
            I2 goes up on way out""",
        stars = '****')

E.rep(28243, 28242, "440 mm, 3.2 s, 20 mA/div -> saturated I2",
        descr = """\
            Very similar to last shot. Saturated signals on I2
            only concern way out.""",
        stars = '****')

# Tim Happel
E.rep(28244, 28243, "204 mm, 1.5 s, 20 mA/div -> 440 mm acc.",
        times = 1.7,
        descr = """\
            Too much power, lots of intermittent arcing.""",
        stars = '*')

E.rep(28245, 28244, "204 mm, 1.0 s, 20 mA/div, preamps 2x (was 5x)",
        times = 1.1,
        posit = 0.10,
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
        ampI2 = CurrentProbe2[20]*Preamp2[2],
        descr = "ELM-ing H-mode. Nice flow patterns during ELMs.",
        stars = "****")

E.rep(28246, 28245, "254 mm, 0.85 s, 20 mA/div, preamps 2x -> L-H transition!",
        times = 0.98,
        posit = 0.15,
        descr = """\
            L-H transition at 0.9 s, right before plunge.
            ELM-free H-mode on way in, ELMing H-mode on way out.""",
        stars = '*****')

# Hendrik Meyer
E.rep(28250, 28246, "304 mm, 4.0 s, 20 mA/div, preamps 2x",
        times = 4.15,
        posit = 0.20,
        descr = """\
            Caught L-H transition at 4.15 s, but no reliable biasing anymore.
            Kepco saturates due to high n. ELMing starts immediately.""",
        stars = '***')

E.rep(28251, 28250, "304 mm, 1.8 s, 20 mA/div, preamps 2x, (bias voltage less positive)",
        times = 1.95,
        descr = """\
            L-H transition at 1.95 s, pretty much at the dwell point in 
            the private flux region. Very high flows on LFS divertor leg 
            (Kepco saturates again). Flows in private flux region reorganize 
            during L-H transition""",
        stars = '****')

E.rep(28252, 28251, "354 mm, 1.75 s, 20 mA/div, preamps 2x -> 1.8 s acc.",
        times = 1.96,
        posit = 0.25,
        descr = """\
            L-H transition at 1.97 s. X-point further inside.
            More arcing than on previous shot""",
        stars = '**')

E.rep(28253, 28252, "304 mm, 1.75 s, 20 mA/div, preamps 2x (repeat 28250)",
        times = 1.9,
        posit = 0.20,
        descr = """\
            L-H transition at 1.93. Saturation and arcing as in last shot.""",
        stars = '**')

E.rep(28254, 28253, "304 mm, 1.75 s, 20 mA/div, preamps 2x (repeat 28251)",
        descr = """\
            L-H transition at 1.935 s, when probe is in private flux region
            (on the way back). Flows in private flux region turn off immediately 
            at transition. H-mode starts with ELM.""",
        stars = '****')


######################
# NEW DAQ FROM HERE ON
######################

############################################
E = campaign.add_experiment(date="20120712")

E.add(28379, "Fixed probe @2564.05 mm", 
        head = head,
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
        ampI2 = CurrentProbe2[20]*Preamp2[2], 
        descr = "Probe very far out, almost no signals.",
        stars = '', **def_XPR_LPS)

E.rep(28380, 28379, "Fixed probe @2569.05 mm -> no data",
        descr = "No plasma.",
        stars = '')

E.rep(28381, 28380, "-200 V bias -> Kepco breaks in at 0.5 s",
        descr = "DC biasing doesn't work.",
        stars = '')

E.rep(28382, 28379, "@2569.10 mm, sweep -> data all the way to 6 s",
        descr = "OK data.",
        stars = '')

E.rep(28383, 28379, "@2589.08 mm, DC offset with small sweep -> worked",
        descr = "Timebase error.",
        stars = '')


############################################
E = campaign.add_experiment(date="20120713")

E.add(28389, "Fixed probe @2569.05 mm, -180 V with sweep", 
        head = head,
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
        ampI2 = CurrentProbe2[20]*Preamp2[2], 
        descr = "Data already on tape.",
        stars = '', **def_XPR_LPS)

E.rep(28390, 28389, "-80 V / 150 V sweep at 100 Hz",
        descr = "Timebase error.",
        stars = '')

E.rep(28394, 28390, "",
        descr = "Data already on tape.",
        stars = '')

E.rep(28395, 28394, "Turn 2nd current probe", 
        ampI1 = CurrentProbe1[20]*Preamp1[2],
        ampI2 = CurrentProbe2[20]*Preamp2[2],
        descr = """\
            Example for DC biasing not working.
            Position calibration seems OK.""",
        stars = '')

# Peter Lang, Francois Ryter
E.rep(28403, 28395, "5 plunges, DC biasing with small rect sweeps",
        times = (0.7, 1.7, 2.7, 3.7, 4.7),
        posit = (0.01, 0.01, 0.01, 0.01, 0.01),
        descr = "Effectively DC biasing. Probe at 1 cm.",
        stars = '*')

E.rep(28404, 28403, "repeat",
        descr = "Similar to last shot.",
        stars = '*')

E.rep(28405, 28404, "4 plunges with decreasing depth",
        times = (0.79, 1.77, 2.74, 4.72),
        posit = (0.10, 0.08, 0.06, 0.04),
        descr = "Second plunge caught some sort of transition.",
        stars = '**')

# Rachael McDermott
E.rep(28406, 28405, "4 plunges to 10 cm",
        times = (0.79, 2.77, 3.78, 4.80),
        posit = (0.10, 0.10, 0.10, 0.10),
        descr = "Low power low density shot, nice fluctuation data.",
        stars = '***')

E.rep(28407, 28406, "5 plunges to 10 cm",
        times = (0.79, 1.78, 2.78, 3.80, 4.79),
        posit = (0.10, 0.10, 0.10, 0.10, 0.10),
        descr = "Low power low density shot, nice fluctuation data.",
        stars = '****')


###############
# REVERSED IpBt
###############

############################################
E = campaign.add_experiment(date="20120717")

E.add(28419, "Fcn gen. 20 Vpp, +8 VDC, 0.5 kHz (let saturate), plunge at 2.7 s", 
        head = head,
        ampI1 = CurrentProbe1[5],
        ampI2 = CurrentProbe2[5], 
        descr = """\
            Died before plunge. Timebase error at 3 s.
            Fcn gen. saturation leads to overshoots and strong
            capacitive pickup.""",
        stars = '', **def_XPR_LPS)

E.rep(28420, 28419, "DAQ test, Fcn gen. 12 Vpp, 4 VDC", 
        descr = "Timebase error at 4 s. No plunges on record",
        stars = '')

E.rep(28421, 28420, "DAQ test",
        descr = "Timebase error at 2.5 s. No plunge on record",
        stars = '')

E.rep(28422, 28421, "DAQ test",
        descr = "Timebase error at 0.8 s. No plunge on record",
        stars = '')

E.rep(28423, 28422, "DAQ test: Plunge at 7 s",
        descr = "No timebase error.",
        stars = '')

E.rep(28424, 28423, "Asym. waveform, 7 Vpp, 1 VDC, 200 Hz, plunge at 7 s",
        descr = "No timebase error.",
        stars = '')

E.rep(28425, 28424, "6.9 Vpp, 1 VDC, 100 Hz, 1 mA/div, plunges at 1.8 and 4.0 s", 
        ampI1 = CurrentProbe1[1],
        ampI2 = CurrentProbe2[1],
        descr = "Plasma died before first plunge.",
        stars = '')

E.rep(28426, 28425, "First data! Kepco breaks in in Isat",
        times = (1.0, 4.0),
        posit = (0.05, 0.05),
        descr = "Not enough voltage swing for fit to work.",
        stars = '*')

E.rep(28427, 28426, "Change sweep pars -> Kepco still breaks",
        posit = (0.10, 0.05),
        descr = "Same as on last shot.",
        stars = '*')

E.rep(28428, 28427, "Back to other fcn gen (DC). -> saturated at 1 mA/div",
        times = 1.0,
        posit = 0.10,
        descr = "DC biasing worked, but saturation on I1.",
        stars = '*')

E.rep(28429, 28428, "Back to 5 mA/div -> no plunge", 
        ampI1 = CurrentProbe1[5],
        ampI2 = CurrentProbe2[5],
        stars = '')

E.rep(28434, 28429, "20 mA/div, 0.1 kHz, all 3 tips on bias voltage", 
        head = headI,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20], 
        ampI3 = CurrentProbe3[20], 
        descr = "Plasma died before first plunge",
        stars = '')

E.rep(28435, 28434, "0.5 kHz, plunge at 1 s",
        times = (1.08, 1.89, 4.09),
        posit = (0.17, 0.17, 0.17),
        descr = "Strong arc at first plunge. Plasma died before 2nd plunge",
        stars = '**')

E.rep(28436, 28435, "1 plunge at 4 s",
        times = 4.09,
        posit = 0.17,
        descr = "Plasma died before plunge",
        stars = '')


############################################
E = campaign.add_experiment(date="20120719")

# Leena's experiment in reversed IpBt
E.add(28442, "0.5 kHz, 3rd pin VF", 
        times = (1.86, 4.06),
        posit = (0.10, 0.10),
        head = head,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20], 
        descr = """\
            1st plunge nice L-mode, 2nd caught disruption on way out.
            Data on way in very similar.""",
        stars = '***', **def_XPR)

E.rep(28444, 28442, "Max plunge at 1.75 s and 3.95 s",
        times = (1.89, 4.09),
        posit = (0.17, 0.17),
        descr = """\
            OK L-mode data on both plunges, but Kepco breaks in
            a bit due to highish density.""",
        stars = '***')

E.rep(28445, 28444, "100 V Kepco",
        times = 1.89,
        posit = 0.17,
        descr = """\
            Nice data on way in, then I1 starts emitting and even saturates.
            100 V Kepco isn't enough. Vf on tip 3 agrees well with VF from 
            sweeps, then goes into emission, approaching Vp.""",
        stars = '****')

E.rep(28446, 28445, "200 V Kepco, 12 Vpp, 5 VDC, 0.5 kHz",
        descr = """\
            200 V Kepco worked much better due to low density, but not enough
            positive voltage for good fits everywhere. I1 goes up on second 
            plunge. Vf on tip3 has same features as on previous shot.""",
        stars = '*****')

# Francois Ryter
E.rep(28448, 28446, "1 kHz, 16 Vpp (digital fcn gen.), VDC from Kepco",
        times = (1.06, 1.86),
        posit = (0.10, 0.10),
        descr = "Unspectacular L-mode data on both plunges.",
        stars = '**')

E.rep(28449, 28448, "0.5 kHz, reduced VDC slightly",
        times = (1.07, 1.88),
        posit = (0.12, 0.12),
        descr = "Unspectacular L-mode data again",
        stars = '**')

E.rep(28450, 28449, "2nd plunge to 1.6 s",
        times = (1.07, 1.67),
        descr = "Also very similar to last two shots",
        stars = '**')

E.rep(28451, 28450, "Max penetration -> shot didn't run",
        times = (1.09, 1.69),
        posit = (0.17, 0.17),
        descr = "Plasma died before first plunge.",
        stars = '')

E.rep(28452, 28451, "Max penetration -> arcs",
        descr = """\
            Nice signals on way in, arcs on way out. Plunges close together
            triggered bug in VPE delay generator programming, causing a VPE
            after 2nd plunge. VF on tip 3 again goes into emission.""",
        stars = '***')


############################################
E = campaign.add_experiment(date="20120720")

E.add(28455, "Acquisition with turned-off Kepco", 
        head = head,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20], 
        descr = "Acquisition works.",
        stars = '', **def_XPR_LPS)

# Rachael McDermott
E.rep(28466, 28455, "0.5 kHz, 16 Vpp, Kepco offset just avoids saturation",
        times = (1.62, 2.62),
        posit = (0.12, 0.12),
        descr = "First in L-mode, second in dying plasma.",
        stars = '**')

# Thomas Puetterich
E.rep(28467, 28466, "Same settings",
        descr = """\
            L-H transitions at 1.565 and 2.670 s. 
            Way too much power, so signals are pretty bad.""",
        stars = '*')

E.rep(28468, 28467, "Plunges to 1.4 and 2.5 s, and only 5 cm",
        times = (1.50, 2.59),
        posit = (0.05, 0.05),
        descr = "Both plunges in high power L-mode before transition.",
        stars = '*')

##########
# HYDROGEN
##########

# Volker Rohde
E.rep(28469, 28468, "H: 14 Vpp, max offset on Kepco",
        descr = "Hydrogen: Both plunges OK and very similar.",
        stars = '**')

########
# HELIUM
########

# Marco Wischmeier
E.rep(28472, 28469, "He: 3 plunges, max penetration", 
        times = (1.09, 1.99, 2.89),
        posit = (0.17, 0.17, 0.17),
        descr = """\
            Good helium data on way in on first plunge, arcs on way out.
            Good flow measurements. Tip 3 goes into emission again.
            Good comparison between DAQs.""",
        stars = '****')

E.rep(28473, 28472, "He again: 1 plunges at 150 mm",
        times = 1.09,
        posit = 0.15,
        descr = """\
            OK data on way in. Not quite enough positive voltage for
            reliable fits, but data very similar to previous shot.""",
        stars = '***')


############################################
E = campaign.add_experiment(date="20120726")

E.add(28504, "Calibration, no signals attached", 
        dig='XPR', head=head, amp_default=amp_default_unity, 
        lines=dict(XPR=dict(amp={}, mapping=mapping_XPR), 
                   LPS=dict(amp={}, mapping=mapping_LPS)),
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        stars = '')

E.rep(28507, 28504, "Calibration, 10 Vpp into 50 Ohm (+/-0.1 A)",
        stars = '')

E.rep(28508, 28507, "Signal also on bias voltage",
        stars = '')


