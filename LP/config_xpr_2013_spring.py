######################
# 2013 SPRING CAMPAIGN
######################

from config_xpr_common import *
from config import Campaign

campaign = Campaign()

tip1 = TipXPR(number=1, pos='lower left', V_keys='ampV1', I_keys='ampI3')
tip2 = TipXPR(number=2, pos='lower right', V_keys='ampV1', I_keys='ampI1')
tip3 = TipXPR(number=3, pos='upper', V_keys='ampV2', I_keys='ampI2')

head = HeadXPR(tips=(tip1, tip2, tip3), R_keys='ampR')


tip1_20130130 = CylindricalTip(r=0.0005, z=0.00303,
        number=1, pos='lower left', V_keys='ampV1', I_keys='ampI3')
tip2_20130130 = CylindricalTip(r=0.0005, z=0.00304,
        number=2, pos='lower right', V_keys='ampV1', I_keys='ampI1')
tip3_20130130 = CylindricalTip(r=0.0005, z=0.00183,
        number=3, pos='upper', V_keys='ampV2', I_keys='ampI2')

head_20130130 = HeadXPR(tips=(tip1_20130130, tip2_20130130, tip3_20130130), R_keys='ampR')


tip1_20130306 = CylindricalTip(r=0.0005, z=0.0026,
        number=1, pos='lower left', V_keys='ampV1', I_keys='ampI3')
tip2_20130306 = CylindricalTip(r=0.0005, z=0.0026,
        number=2, pos='lower right', V_keys='ampV1', I_keys='ampI1')
tip3_20130306 = CylindricalTip(r=0.0005, z=0.002,
        number=3, pos='upper', V_keys='ampV2', I_keys='ampI2')

head_20130306 = HeadXPR(tips=(tip1_20130306, tip2_20130306, tip3_20130306), R_keys='ampR')


tip4_20130312 = CylindricalTip(r=0.0005, z=0.0052,
        number=4, pos='lower', V_keys='ampV1', I_keys='ampI4')

head_20130312 = HeadXPR(tips=head_20130306.tips + (tip4_20130312,), R_keys='ampR')


fact = 4 * 5.54630/27. / 2**16
offs = -47578968*fact - 0.105

amp_XPR_pos = dict(
            ampR  = Amp(fact=fact, offs=offs),
            ampV1 = Amp(fact=100., offs=-68.48),
            ampV2 = Amp(fact=100., offs=-68.13))

mapping_XPR_pos = dict(
            ampR  = 'Pos', 
            ampV1 = 'S1', 
            ampV2 = 'S3',
            ampI1 = 'S4',
            ampI2 = 'S2',
            ampI3 = 'S6',
            ampVF = 'S6',
            ampI4 = 'S7')

lines_XPR_pos = dict(amp=amp_XPR_pos, mapping=mapping_XPR_pos)

def_XPR_pos = dict(dig='XPR_pos', amp_default=amp_default, lines=dict(XPR=lines_XPR_pos))



############################################
E = campaign.add_experiment(date="20130124")

E.add(29289, "DAQ test, Mach DC, single 1 kHz, 13.5 Vpp, DAQ 5 s",
        head = head,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            Everything OK.""",
        stars = '', **def_XPR_pos)

E.rep(29290, 29289, "Std H-mode, 10 cm, 1 s, single 1 kHz sine, 13.5 Vpp",
        times = 1.0,
        posit = 0.1,
        descr = """\
            Went well.""",
        stars = '**')

# Andreas Burkhart
E.rep(29291, 29290, "single 1 kHz triangle, 12.5 Vpp, more negative",
        descr = """\
            Tip 1 arced. Sweeps on tip 3 look better. Density higher than
            on last shot.""",
        stars = '*')

# Gregor Birkenmeier
E.rep(29302, 29291, "600 kA, 2.5 T, 1.3 MW NBI, all the way through",
        times = 3.4,
        posit = 0.34,
        descr = """\
            Beam came on at 3.4 s, not 3.5 s. Mach tip 1 arcs almost 
            immediately, single tip shortly after.""",
        stars = '**')

E.rep(29303, 29302, "B = 1.8 T, 1.3 MW NBI, 3.3 s, all on sweeps, 13.5 Vpp, 1 kHz",
        times = 3.3,
        posit = 0.34,
        descr = """\
            Second 100 V power supply was off, so swept Mach probes are useless.
            However, good data from single tip all the way through, with L-H
            transition when the probe was at 22 cm. NO arcing. Probe did not
            appear to alter density rise, which appears to be stopped by large
            ELMs.""",
        stars = '***')

E.rep(29304, 29303, "B = 1.2 T, switch 2nd Kepco on, didn't run",
        times = 3.3,
        posit = 0.34,
        descr = "",
        stars = '')

E.rep(29305, 29304, "Repeat, didn't run",
        times = 3.3,
        posit = 0.34,
        descr = "",
        stars = '')

E.rep(29306, 29305, "1.4 T, 800 kW",
        times = 3.3,
        posit = 0.34,
        descr = "Good data. L-H transition at 3.432 s.",
        stars = '****')

# Stefan Muller
E.rep(29307, 29306, "Repeat with Mach at -200 V",
        times = 3.3,
        posit = 0.34,
        descr = "Great! L-H transition at 3.427 s",
        stars = '*****')

# Gregor Birkenmeier
E.rep(29308, 29307, "B = 3.2 T, 1.6 MW, all on sweeps, single more pos",
        times = 3.3,
        posit = 0.34,
        descr = "Overheat pretty early. Late L-H transition.",
        stars = '**')

E.rep(29309, 29308, "B = 2.5 T, 1.6 MW, all on sweeps, single more neg",
        times = 3.3,
        posit = 0.34,
        descr = "Better. Mach signals OK.",
        stars = '***')

E.rep(29310, 29309, "400 kA, B = 1.4 T, 0.8 MW, Mach -200 V, 50 ms later, sweep at 0.4 kHz",
        times = 3.35,
        posit = 0.34,
        descr = "Already in H-mode, but nice data all the way across.",
        stars = '****')

# Stefan Muller
E.rep(29311, 29310, "600 kA, B = 1.8 T, 1.0 MW, sweep at 0.7 kHz",
        times = 3.35,
        posit = 0.34,
        descr = """\
            Already in H-mode again, nice data almost until dwell. LFS density
            drop marks onset of coherent oscillations.""",
        stars = '*****')

E.rep(29312, 29311, "Repeat, plunge 50 ms earlier, sweep at 1 kHz",
        times = 3.30,
        posit = 0.34,
        descr = """\
            Arced again, L-mode data complimentory to H-mode data in previous shot.
            Single tip still got L-H transition on HFS! Strong drop in HFS density at
            transition""",
        stars = '*****')

# Gregor Birkenmeier
E.rep(29313, 29312, "2.1 MW at 3.0 s, so plunge at 2.7 s",
        times = 2.70,
        posit = 0.34,
        descr = "Didn't run",
        stars = '')

E.rep(29315, 29313, "Repeat 29313",
        times = 2.70,
        posit = 0.34,
        descr = """\
            Nice L-mode data up to dwell. 
            Then I2 and I3 creep up to saturation.""",
        stars = '****')



############################################
E = campaign.add_experiment(date="20130125")

E.add(29319, "Mach DC, single 1 kHz, 13.5 Vpp, DAQ 5 s",
        head = head,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            DAQ OK.""",
        stars = '', **def_XPR_pos)

# Daniel Carralero
E.rep(29320, 29319, "3 plunges, Mach at -200 V",
        times = (2.4, 3.1, 3.8),
        posit = (0.34, 0.34, 0.34),
        descr = """\
            Plasma was in H-mode and probe pulled it back to L-mode.
            Oscillations are present in all 3 H-mode phases.
            Unfortunately, probe arced before H-L transition.""",
        stars = '***')

E.rep(29321, 29320, "1 plunge, all on sweeps",
        times = 3.8,
        posit = 0.34,
        descr = """\
            H-L back-transition on LFS divertor leg, but with swept Mach probes.
            Arcs when entering the HFS.""",
        stars = '***')

E.rep(29322, 29321, "Less density",
        descr = """\
            Arcs similar to last shot. Still some H-mode like behavior on way in.
            Divertor currents and Tdiv modulated with core MHD.""",
        stars = '***')

E.rep(29323, 29322, "Less density now",
        descr = """\
            Very low density. No arcs, but not enough positive Voltage.
            No perturbation of the probe visible in DCN and Tdiv.""",
        stars = '****')

E.rep(29324, 29323, "n = 2.5e19, Mach at -200 V",
        descr = """\
            Tip 1 arcs at 28 cm. Nice fluctuation data before that.
            Profiles also OK up to the same positions.
            Extremely low fluctuation amplitude.
            Clear non-saturation of Isat on HFS.""",
        stars = '****')

E.rep(29325, 29324, "n = 6e19, All swept at 14 Vpp (more pos)",
        descr = """\
            Didn't run.""",
        stars = '')

E.rep(29326, 29325, "Repeat",
        descr = """\
            Nice high density data. Lower density on HFS!
            Too high Mach numbers again.
            I1 goes from 1.5 A on LFS to almost 0 on HFS, and back to 1.5 A on way out!
            Probably a very good demonstration that these low currents are real!""",
        stars = '*****')

# DAQ test
E.rep(29336, 29326, "DAQ test",
        descr = "",
        stars = '')

# Hans-Werner Mueller
E.rep(29337, 29336, "n = 1.5e19, 1 MA, 600 kW ECRH, same sweep settings, 1 plunge at 2.3 s",
        times = 2.3,
        posit = 0.34,
        descr = """\
            Very nice profiles. Ran a bit out of positive sweep voltage.
            Single tip creeps up again on way out.""",
        stars = '****')

E.rep(29338, 29337, "Mach at -200 V, single swept at 14.5 Vpp, 1 kHz",
        descr = """\
            Nice fluctuation data on Mach probe, no arcs.""",
        stars = '*****')

E.rep(29339, 29338, "Increase sweep amplitude to 15.5 Vpp, 1 kHz",
        descr = """\
            Good data again. 15.5 Vpp was too high, since Kepco starts
            oscillating. Similar to last shot, but Te significantly higher
            on LFS and density higher on HFS.""",
        stars = '*****')


############################################
E = campaign.add_experiment(date="20130131")

# New tips lengths after change on 20130130: 
#
# lower right: 3.04 mm
# lower left: 3.03 mm
# upper: 1.83 mm

E.add(29377, "DAQ test, all tips on sweeps at 14.5 Vpp",
        head = head_20130130,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            DAQ OK.""",
        stars = '', **def_XPR_pos)

# Peter Lang
E.rep(29378, 29377, "Test plunge at 0.9 s, to 12 cm, sweeps at 13.5 Vpp",
        times = 0.9,
        posit = 0.12,
        descr = """\
            All signals there. Arcs immediately after entering plasma, but
            normal measurements already before dwell.""",
        stars = '*')


############################################
E = campaign.add_experiment(date="20130205")

E.add(29400, "DAQ test, all tips on sweeps at 13.5 Vpp",
        head = head_20130130,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            DAQ OK.""",
        stars = '', **def_XPR_pos)

# CTS
E.rep(29401, 29400, "Cond. at 6.0 s, 10 cm, 50 mA/mV to avoid saturation",
        times = 6.0,
        posit = 0.1,
        ampI1 = CurrentProbe1[50],
        ampI2 = CurrentProbe2[50],
        ampI3 = CurrentProbe3[50], 
        descr = """\
            OK. Already in rampdown. No arcs.""",
        stars = '*')

E.rep(29402, 29401, "Go to 5.5 s, 25 cm",
        times = 5.5,
        posit = 0.25,
        descr = """Didn't run.""",
        stars = '')

E.rep(29403, 29402, "Repeat",
        descr = """\
            Better. Now some mini arcs.""",
        stars = '**')

E.rep(29404, 29403, "All the way thru",
        posit = 0.34,
        descr = """\
            Programming error led to probe not retracting, VPE at 6.02122 s.
            Data appears to be H-mode, high density, ELMs. Rough.""",
        stars = '***')

E.rep(29405, 29404, "25 cm",
        posit = 0.25,
        descr = """\
            No strong arcs. Still mini-arcs on tip 1 (still not clean).""",
        stars = '***')

E.rep(29406, 29405, "30 cm, 5.9 s, 20 mA/mV",
        times = 5.9,
        posit = 0.3,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            Good data, but plasma disrupts at dwell.""",
        stars = '***')

# G. Birkenmeier
E.rep(29407, 29406, "Repeat of 29315, 34 cm, 2.9 s, Mach -200 V",
        times = 2.9,
        posit = 0.34,
        descr = """\
            Unfortunately arc on Mach tip 1 already pretty early.
            Should be L-mode, but it's not very clear.""",
        stars = '***')

E.rep(29408, 29407, "Repeat with more power, 2.725 s, all on sweeps",
        times = 2.725,
        posit = 0.34,
        descr = """\
            Missed H-mode, very nice L-mode data on way in and out. No arcs.
            Mach numbers not higher than 1 on HFS. Profile match exactly
            on way in and out.""",
        stars = '*****')


############################################
E = campaign.add_experiment(date="20130208")

# Stefan Muller
E.add(29438, "DAQ test, single tip on sweeps at 13.5 Vpp, Mach at -200 V",
        head = head_20130130,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            """,
        stars = '', **def_XPR_pos)

E.rep(29439, 29438, "Ref. 29311, beams on at 1.1, 2.7, 4.3 s",
        times = (1.15, 2.75, 4.35),
        posit = (0.34, 0.34, 0.34),
        descr = """\
            Arcs before dwell, but OK data. Discharge had modes that severely
            impacted the H-mode.""",
        stars = '***')

E.rep(29440, 29439, "Repeat with single tip floating",
        times = (1.15, 2.75, 4.35),
        posit = (0.34, 0.34, 0.34),
        descr = """\
            Arc burns all the way through shot!""",
        stars = '*')

E.rep(29441, 29440, "Repeat with single tip at -150 V, 20 ms earlier,",
        times = (1.13, 2.73, 4.33),
        posit = (0.34, 0.34, 0.34),
        descr = """\
            Arc burns all the way through shot! Ok data on single tip.""",
        stars = '*')

# Standard H Mode
E.rep(29442, 29441, "Probe test, all on sweeps",
        times = 6.0,
        posit = 0.05,
        descr = """\
            Tip 1 resistive, arcs even before plasma.""",
        stars = '*')

# Tests
E.rep(29443, 29442, "No plunge, but standing at waiting position",
        times = (),
        posit = (),
        descr = """\
            Tip 1 resistive, arcs even before plasma.""",
        stars = '')


############################################
E = campaign.add_experiment(date="20130307")

# new probe head test
E.add(29676, "New probe head!",
        times = 6.5,
        posit = 0.05,
        head = head_20130306,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            All signals there. Lots of arcs due to dirt.""",
        stars = '', **def_XPR_pos)

E.rep(29681, 29676, "Same test as before.", 
        times = 6.8,
        descr = """\
            Plasma died before plunge.""",
        stars = '')

E.rep(29683, 29681, "Go to 1 s", 
        times = 1.0,
        descr = """\
            Good. No arcs.""",
        stars = '**')


############################################
E = campaign.add_experiment(date="20130308")

# arc protection tests
E.add(29690, "Single tip current signal passes thru arc protection box",
        head = head_20130306,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            Compared to 29681, current signal S2 jumped up by 4 mV 
            and is a bit noisier.""",
        stars = '', **def_XPR_pos)

E.rep(29691, 29690, "Single tip signal directly on S2",
        descr = "Offset remained, noise level back to previous levels.",
        stars = '')

E.rep(29692, 29691, "Single tip signal teed outside arc box, all tips on 4 A Kepco",
        descr = """\
            Voltage calibration correct. Lower noise level on S2 if teed before
            arc box.""",
        stars = '')

# Std H-mode
E.rep(29693, 29692, "Single tip back on 1 A Kepco",
        times = 1.0,
        posit = 0.15,
        descr = """\
            Good data. Arc box switched off prematurely, failed to switch back
            on several times, but eventually succeeded.""",
        stars = '**')

E.rep(29694, 29693, "Repeat, go all the way thru",
        times = 1.0,
        posit = 0.34,
        descr = """\
            Arc box switched off at different conditions again. Real arcs on I2.
            Otherwise good data.""",
        stars = '***')

# Martin Oberkofler
E.rep(29695, 29694, "Full plunge at 3.8 s, arc box on I3 triggered by I2",
        times = 3.8,
        posit = 0.34,
        descr = """\
            No arcs on Mach probe. Arc box switched off I3 for no apparent reason.""",
        stars = '***')

E.rep(29696, 29695, "No plunge. Arc box switches 9 V battery on S7",
        times = (),
        posit = (),
        descr = "Arc box did not switch.",
        stars = '')

E.rep(29697, 29696, "Plunge again, arc box switches 9 V battery on S8",
        times = 3.8,
        posit = 0.34,
        descr = "Nice data, no arcs. Arc box switched twice.",
        stars = '****')

E.rep(29698, 29697, "Repeat",
        times = 3.8,
        posit = 0.34,
        descr = "Arcs on single tip on way out.",
        stars = '****')


############################################
E = campaign.add_experiment(date="20130312")

# probe test
E.add(29699, "Arc box switches 4 A Kepco (Mach tips), triggered by I2 at 1.9 V",
        head = head_20130306,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = "Everything looks OK.",
        stars = '', **def_XPR_pos)

# Std H-mode
E.rep(29700, 29699, "15 cm at 1 s",
        times = 1.0,
        posit = 0.15,
        descr = "All signals good. No arcs. Arc box didn't switch.",
        stars = '**')

# Std H-mode
E.rep(29701, 29700, "Go to 25 cm",
        times = 1.0,
        posit = 0.25,
        descr = "Single tip arcs, I2 very high. Arc box didn't switch.",
        stars = '**')

# Std H-mode
E.rep(29702, 29701, "Go to 18 cm",
        times = 1.0,
        posit = 0.18,
        descr = """\
            Additionally 1 MW ECRH at 0.9 s. Arcs from I1, which is *not*
            monitored. Positive currents still present after leaving plasma,
            negatives missing. Curious.""",
        stars = '*')

# rt-ECCD FB for disruption avoidance
E.rep(29703, 29702, "Check if everything still works",
        times = 1.0,
        posit = 0.05,
        descr = "All signals still there.",
        stars = '**')

# Andrea Scarabosio
E.rep(29708, 29703, "Two full plunges, all on sweeps at 17.0 Vpp at 1 kHz",
        times = (2.0, 4.0),
        posit = (0.34, 0.34),
        descr = """\
            First plunge excellent, second in accidental gas puff. Arc protection
            switched. Currents flowing between Mach tips in disconnected mode!""",
        stars = '****')

E.rep(29709, 29708, "Repeat w/o density increase, 2nd plunge at 3.9 s (didn't run)",
        times = (2.0, 3.9),
        posit = (0.34, 0.34),
        descr = "Didn't run",
        stars = '')

E.rep(29710, 29709, "Repeat",
        times = (2.0, 3.9),
        posit = (0.34, 0.34),
        descr = """\
            Both plunges excellent.""",
        stars = '*****')

E.rep(29711, 29710, "More gas puff, single tip on VF, I3 on total current from 4 A Kepco",
        times = (2.0, 3.9),
        posit = (0.34, 0.34),
        descr = """\
            First plunge excellent. Plasma starving during 2nd plunge.""",
        stars = '*****')

E.rep(29713, 29711, "I3 on 50 mA/div, no plunge",
        times = (),
        posit = (),
        ampI3 = CurrentProbe3[50],
        descr = "",
        stars = '')

E.add(29715, "Even higher n, config as before, total Mach currents on I4",
        times = (2.0, 3.9),
        posit = (0.34, 0.34),
        head = head_20130312,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20],
        ampI4 = CurrentProbe4[50],
        descr = """\
            First plunge good on way in. Second plunge in starving plasma.""",
        stars = '****', **def_XPR_pos)

E.rep(29717, 29715, "Ramp n from medium to high",
        descr = """\
            Both plunges good on way in.""",
        stars = '****')

E.rep(29721, 29717, "Change sensitivity to avoid saturation",
        ampI1 = CurrentProbe1[50],
        ampI2 = CurrentProbe2[50],
        ampI3 = CurrentProbe3[50],
        ampI4 = CurrentProbe4[100],
        descr = """\
            Again both plunges good on way in.""",
        stars = '****')

# Patrick Simon
E.rep(29722, 29721, "Limiter shot for GAM studies, sweeps at 13.5 Vpp",
        times = 1.7,
        posit = 0.34,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20],
        ampI4 = CurrentProbe4[50],
        descr = "Almost nothing down there.",
        stars = '')


############################################
E = campaign.add_experiment(date="20130314")

# probe test
E.add(29727, "I1 on 1 A Kepco, I2 on 4 A Kepco (13.5 Vpp at 1 kHz), single tip VF",
        head = head_20130312,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        ampI4 = CurrentProbe4[50],
        descr = """\
            Voltage of 1 A Kepco on I1. 
            Arc box on I2 (clicking started at 3.08 Vpp).
            I4 measures from 4 A Kepco into arc box.""",
        stars = '', **def_XPR_pos)

# Std H-Mode
E.rep(29728, 29727, "Test plunge at 10 cm",
        times = 1.0,
        posit = 0.1,
        descr = "Didn't run.",
        stars = '')

E.rep(29729, 29728, "Repeat",
        times = 1.0,
        posit = 0.1,
        descr = "Wimpy plasma. Not enough el. current.",
        stars = '**')

# Martin Oberkofler
E.rep(29730, 29729, "No plunge",
        times = (),
        posit = (),
        descr = "",
        stars = '')

E.rep(29731, 29730, "Full stroke at 3.8 s",
        times = 3.8,
        posit = 0.34,
        descr = """\
            Very similar data to 29697. So it didn't matter whether the two tips
            are connected to the same or to different power supplies. Arc box on
            I2 triggered twice: I2 goes to 0 and I1 is unaffected, as hoped.""",
        stars = '*****')

# Stefan Muller
E.rep(29733, 29731, "Full stroke at 1.1 and 2.7 s for 1st and 2nd L-H transition",
        times = (1.1, 2.7),
        posit = (0.34, 0.34),
        descr = """\
            Arc box on I4 at 4.5 Vpp. Arc box didn't trigger.
            First plunge: L-H and H-L transition while probe is on HFS.
            Coherent oscillations appear only in H-mode.""",
        stars = '*****')

E.rep(29734, 29733, "Full stroke at 1.1 and 4.3 s for 1st and 3rd L-H transition",
        times = (1.1, 4.3),
        posit = (0.34, 0.34),
        descr = """\
            Arc box on I4 at 3.4 Vpp. Arc box didn't trigger.
            L-H transitions less clear than on previous shot.""",
        stars = '*****')

E.rep(29735, 29734, "Density 3e19, 3 strokes, 20 ms earlier",
        times = (1.12, 2.72, 4.32),
        posit = (0.34, 0.34, 0.34),
        descr = """\
            Arc box on I4 at 2.36 Vpp. Arc box triggered correctly!
            Clearer L-H transitions, but more problems with arcing, as expected.
            Arc box saved nice L-H transitions data on plunge 1.""",
        stars = '*****')

# Tests
E.rep(29745, 29735, "I4 measures current in ground on I2, no plunge",
        times = (),
        posit = (),
        descr = """\
            No currents on ground. (Found out later that the manipulator
            is not connected with the tokamak, so clearly no currents
            can flow in these grounds.)""",
        stars = '')


############################################
E = campaign.add_experiment(date="20130319")

# Remove ground connection and digital position signal from
# measurement rack, so that it is only grounded via the tips.
E.add(29755, "Mach -200 V, single tip swept at 1 kHz at 13.5 Vpp.",
        head = head_20130312,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        ampI4 = CurrentProbe4[50],
        descr = """\
            Arc box setting as on previous experiment on 2013/03/14.
            No surprises. Noise level probably even lower than usual.""",
        stars = '', **def_XPR_pos)

E.rep(29756, 29755, "Std H-mode with removed rack grounding, 15 cm at 1 s",
        times = 1.0,
        posit = 0.15,
        descr = """\
            Total measured current was 0!""",
        stars = '')

# Garrard Conway
E.rep(29780, 29756, "Grounds reconnected, full stroke at 3.8 s",
        times = 3.8,
        posit = 0.34,
        descr = """\
            No I-phase. Good data except for switched-off arcs on LFS divertor
            leg.""",
        stars = '****')

E.rep(29781, 29780, "Repeat 29779, plunge at 2.9 s",
        times = 3.9,
        posit = 0.34,
        descr = """\
            Probably I-phase. Looks like L-mode where coherent oscillations appear
            intermittently. No arcs on way in.""",
        stars = '*****')


############################################
E = campaign.add_experiment(date="20130322")

# Matthias Bernert
E.add(29811, "Mach -200 V, single tip swept at 1 kHz at 13.5 Vpp.",
        times = 0.5,
        posit = 0.05,
        head = head_20130312,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        ampI4 = CurrentProbe4[50],
        descr = "ELM data, but camera didn't see us.",
        stars = '**', **def_XPR_pos)

E.rep(29812, 29811, "Go in at 0.35 s, dwell for 200 ms",
        times = 0.35,
        posit = 0.05,
        descr = """\
            Large type I ELMs. Nice flow data! Low density, high Te""",
        stars = '****')

E.rep(29813, 29812, "Go to 7.5 cm, stay for 300 ms",
        times = 0.35,
        posit = 0.075,
        descr = """\
            Lots of coherent oscillations plus ELMs. Some arcs that were
            switched off nicely. Again low density, high Te.""",
        stars = '****')

E.rep(29814, 29813, "Go back to 5 cm, dwell for 200 ms",
        times = 0.35,
        posit = 0.05,
        descr = """\
            Good L-H transition data again. Fewer arcs.""",
        stars = '***')

E.rep(29815, 29814, "Repeat",
        times = 0.35,
        posit = 0.05,
        descr = """\
            Similar data again.""",
        stars = '***')

E.rep(29816, 29815, "Repeat",
        times = 0.35,
        posit = 0.05,
        descr = """\
            Similar again.""",
        stars = '***')

# Tim Happel
E.rep(29817, 29816, "All the way thru at 0.34 s",
        times = 3.4,
        posit = 0.34,
        descr = """\
            High density. LFS measurements good up to 15 cm, then arc box
            in constant action. On way out good data again.""",
        stars = '***')

E.rep(29818, 29817, "Go to 20 cm",
        times = 3.4,
        posit = 0.2,
        descr = """\
            L-mode on way in, apparently H-mode with coherent oscillations on
            way out. Arc box switches on LFS divertor leg. Nice L-H transition
            in PF region!""",
        stars = '*****')

E.rep(29819, 29818, "Plasma lower",
        times = 3.4,
        posit = 0.2,
        descr = """\
            Much worse than previous shot, since plasma was sitting more in
            divertor. Sum signal I4 does not go to 0 after shot, but
            I1 and I2 do!""",
        stars = '***')

# Test acquisition between 29819 and 29820: Sweeps everywhere.
# Test shot 03200: Everything seems OK

E.rep(29820, 29819, "USN, all the way thru",
        times = 3.4,
        posit = 0.34,
        descr = """\
            Good data!""",
        stars = '***')

E.rep(29821, 29820, "USN, all the way thru, twice",
        times = (2.4, 3.4),
        posit = (0.34, 0.34),
        descr = """\
            Good data, similar to last shot.""",
        stars = '***')



