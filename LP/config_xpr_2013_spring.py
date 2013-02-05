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
            ampVF = 'S6')

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
        descr = "Already in H-mode again, nice data almost until dwell.",
        stars = '****')

E.rep(29312, 29311, "Repeat, plunge 50 ms earlier, sweep at 1 kHz",
        times = 3.30,
        posit = 0.34,
        descr = """\
            Arced again, L-mode data complimentory to H-mode data in previous shot.""",
        stars = '***')

# Gregor Birkenmeier
E.rep(29313, 29312, "2.1 MW at 3.0 s, so plunge at 2.7 s",
        times = 2.70,
        posit = 0.34,
        descr = "Didn't run",
        stars = '')

E.rep(29315, 29313, "Repeat 29313",
        times = 2.70,
        posit = 0.34,
        descr = "Nice L-mode data all the way through.",
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
E.rep(29320, 29319, "3 plunges",
        times = (2.4, 3.1, 3.8),
        posit = (0.34, 0.34, 0.34),
        descr = """\
            Plasma was in H-mode and probe pulled it back
            to L-mode.""",
        stars = '**')

E.rep(29321, 29320, "1 plunge, all on sweeps",
        times = 3.8,
        posit = 0.34,
        descr = """\
            Arcs at PF region.""",
        stars = '**')

E.rep(29322, 29321, "Less density",
        descr = """\
            Arcs similar to last shot. Density wasn't really lower.""",
        stars = '**')

E.rep(29323, 29322, "Less density now",
        descr = """\
            Very low density. No arcs, but not enough positive Voltage.""",
        stars = '****')

E.rep(29324, 29323, "n = 2.5e19, Mach at -200 V",
        descr = """\
            Tip 1 arcs at 28 cm. Nice fluctuation data before that.
            Profiles also OK up to the same positions.""",
        stars = '***')

E.rep(29325, 29324, "n = 6e19, All swept at 14 Vpp (more pos)",
        descr = """\
            Didn'r run.""",
        stars = '')

E.rep(29326, 29325, "Repeat",
        descr = """\
            Nice high density data. Too high Mach numbers again.""",
        stars = '****')

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
            oscillating. Similar to last shot, but HFS density twice as much
            as on previous shot.""",
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
            Unfortunately arc on Mach tip 1 already pretty early.""",
        stars = '***')

E.rep(29408, 29407, "Repeat with more power, 2.725 s, all on sweeps",
        times = 2.725,
        posit = 0.34,
        descr = """\
            Missed H-mode, very nice L-mode data on way in and out.
            No arcs.""",
        stars = '****')


