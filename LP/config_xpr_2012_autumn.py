######################
# 2012 AUTUMN CAMPAIGN
######################

from config_xpr_common import *
from config import Campaign

campaign = Campaign()

# full range again, and cart is 5 mm further in than in spring
#fixpoints = (1.436, -0.101), (7.06, 0.341)

# the above is off by -2 cm, so take visible wall crossing from 28799 (full stroke)
fixpoints = (2.46, 0), (7.06, 0.34)
amp_XPR['ampR'] = Amp(fixpoints=fixpoints)

tip1 = TipXPR(number=1, pos='lower left', V_keys='ampV1', I_keys='ampI3')
tip2 = TipXPR(number=2, pos='lower right', V_keys='ampV1', I_keys='ampI1')
tip3 = TipXPR(number=3, pos='upper', V_keys='ampV1', I_keys='ampI2')

headI = HeadXPR(tips=(tip1, tip2, tip3), R_keys='ampR')

tip3sep = TipXPR(number=3, pos='upper', V_keys='ampV2', I_keys='ampI2')

headI_tip3sep = HeadXPR(tips=(tip1, tip2, tip3sep), R_keys='ampR')

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
E = campaign.add_experiment(date="20121011")

E.add(28633, "DAQ test",
        head = headI,
        stars = '', **def_XPR)

E.add(28634, "Sweep attached to 2x100 V Kepco pair, all tips on sweep, plunges behind wall",
        head = headI,
        times = (0.950, 1.850),
        ampI1 = CurrentProbe1[5000],
        ampI2 = CurrentProbe2[5000],
        ampI3 = CurrentProbe3[5000], 
        descr = "Fuse blown on whole Kepco rack. No data", 
        stars = '', **def_XPR)

E.rep(28636, 28634, "Switch Kepcos off", 
        descr = "No motion. No signals?", 
        stars = '')

E.rep(28637, 28636, "Acquire trigger signals", 
        descr = "No motion. No signals?",
        stars = '')

E.rep(28641, 28637, "TTL via LWL 1061 on channel 5",
        descr = "Nothing came through",
        stars = '')

E.rep(28643, 28641, "Sine via fcn gen on channel 5",
        descr = "",
        stars = '')

E.rep(28645, 28643, "Kepcos on separate trafo, sweep on",
        descr = "",
        stars = '')

E.rep(28646, 28645, "Change sensitity",
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = "", 
        stars = '')

E.rep(28647, 28646, "Repeat",
        stars = '')

E.rep(28648, 28647, "Reset local timer for PPG TS06", 
        descr = "Worked, but tip 3 apparently short circuits",
        stars = '')

E.rep(28649, 28648, "Go to three plunges behind the wall", 
        times = (1.0, 2.0, 3.0),
        descr = """\
            Short circuit on tip 3 on plunge 0, then
            UCSD Kepco trips, then current 3 follows voltage.""",
        stars = '')

E.rep(28650, 28649, "Take tip 3 off bias voltage, 3rd plunge slower",
        descr = "3rd plunge didn't come out, position signal remains noisy",
        stars = '')

E.rep(28651, 28650, "Position signal on S8", 
        XPR_mapping_ampR = 'S8',
        descr = "Position signal just as noisy",
        stars = '')


############################################
E = campaign.add_experiment(date="20121016")

E.add(28657, "No plunges",
        head = headI,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = "Noisy position signal", 
        stars = '', **def_XPR)

E.rep(28668, 28657, "No plunges",
        descr = "Noisy position signal, tip 3 ok",
        stars = '')

E.rep(28669, 28668, "First plasma plunge of the season", 
        times = 3.4, 
        descr = """\
            Too late, missed plasma. 
            Tip 3 current follows bias voltage with resistance of ~1 kOhm.""",
        stars = '')

E.rep(28670, 28669, "Plunge at beginning of second heating phase", 
        times = (2.35,), 
        descr = """\
            Probe in plasma for first time, many small arcs on way in.
            Caught L-H transition. Wiggle in VF, at transition?""",
        stars = '***')

# Matthias Willensdorfer
E.rep(28671, 28670, "3 plunges", 
        times = (1.3, 2.0, 3.15), 
        descr = """\
            OK, current goes up on 3rd plunge.
            Current on upper tip roughly equal to sum of Mach tips.""",
        stars = '***')

E.rep(28672, 28671, "repeat", 
        times = (1.3, 2.0, 3.15), 
        descr = "More arcs than on last shot",
        stars = '')

E.rep(28673, 28672, "Only one plunge", 
        times = (3.5,), 
        descr = "Went better",
        stars = '')

E.rep(28674, 28673, "Two plunges", 
        times = (1.9, 3.5), 
        descr = "2nd Kepco failed between shots",
        stars = '')

E.rep(28675, 28674, "Three plunges", 
        times = (0.9, 1.7, 3.1), 
        descr = "2nd Kepco failed again between shots",
        stars = '')


############################################
E = campaign.add_experiment(date="20121025")

E.add(28747, "After repair of short circuit, 3 plunges to 2 cm",
        times = (0.9, 1.7, 3.1),
        posit = (0.02, 0.02, 0.02),
        head = headI,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = "Cleaning arcs on 1st plunge, others OK", 
        stars = '', **def_XPR)

E.rep(28753, 28747, "No plunges, only test if Kepcos still work",
        stars = '')

# Tilman Lunt
E.rep(28754, 28753, "Gas puff imaging: Two plunges. Thermography on X-point",
        times = (1.6, 2.9),
        posit = (0.16, 0.16),
        descr = "Tips still dirty. OK data on both plunges.",
        stars = '**')

E.rep(28755, 28754, "All the way",
        times = (1.6, 2.7),
        posit = (0.34, 0.34),
        descr = "Shot didn't run",
        stars = '')

E.rep(28756, 28755, "Try again",
        descr = "Very nice L-mode data",
        stars = '****')


############################################
E = campaign.add_experiment(date="20121031")

# Leena
E.add(28794, "Ref. 27692, 1 plunge at 2.8 s, all the way through",
        times = 2.8,
        posit = 0.34,
        head = headI,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = "Nice data, but arcs and current limit", 
        stars = '****', **def_XPR)

E.rep(28795, 28794, "N2 FF",
        times = 4.0,
        descr = "Slightly better than last shot",
        stars = '****')

E.rep(28796, 28795, "More N2 FF (puff without limit, disruption)",
        descr = "Disruption before plunge",
        stars = '')

E.rep(28797, 28796, "Less N2 FF",
        descr = "Very low signal, but everything OK",
        stars = '*****')

E.rep(28798, 28797, "N2 FF 3.4e21/s",
        descr = "Same as last shot, except for density",
        stars = '*****')

E.rep(28799, 28798, "N2 FF 2e21/s (1 valve)",
        descr = """\
            Almost identical to 28795. 
            On way out, cable 3 ripped inside vacuum,
            leading to loss of lower-left tip on channel S6""",
        stars = '****')

############################################
E = campaign.add_experiment(date="20121106")

# Leena
E.add(28818, "1 plunge at 2.8 s, all the way through, NO flow measurement",
        times = 2.8,
        posit = 0.34,
        head = headI,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            No flow measurement due to ripped cable on 28799.
            Otherwise nice data, but arcs and current limit""",
        stars = '***', **def_XPR)

E.rep(28819, 28818, "Go only to X-point",
        descr = "OK data with upper tip up to X-point.",
        stars = '**')

E.rep(28820, 28819, "No plunge",
        descr = "",
        stars = '')


############################################
E = campaign.add_experiment(date="20121122")

# DC biasing testing - NO POSITION SIGNAL
E.add(28871, "DC biasing test: -200 V",
        times = 3.2,
        posit = 0.05,
        head = headI,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            Got only -100 V. DC voltage dropped as expected""",
        stars = '', **def_XPR)

E.rep(28872, 28871, "Only line 1 on bias voltage; line 1 not on V-divider",
        descr = "Line 1 signal departs from Isat after short time")

E.rep(28873, 28872, "Sweep, line 1 back on V-divider",
        descr = "Sweep signal visible")

E.rep(28874, 28873, "All tips on sweeps",
        descr = "Works.")

E.rep(28875, 28874, "Isolated tip at -200 V from 2nd Kepco",
        descr = "No plunge")

E.rep(28876, 28875, "Repeat, this time with plunge",
        descr = "Isat on single tip worked!")

E.rep(28877, 28876, "Now Mach tips on VDC, sweep on upper tip",
        descr = "Isat on single tip worked!")

E.rep(28878, 28877, "Kepco of single tip also on VDC - HST test, no plasma")

E.rep(28879, 28878, "5x9V batteries on single tip! - HST test, no plasma")

E.rep(28880, 28879, "5x9V batteries on single tip! Plasma", 
        descr = "Strong variations on battery voltage measurement!")

E.rep(28883, 28880, "5x9V batteries on single tip, no connection to probe", 
        descr = """\
            Voltage divided by 200, as expected, not 100. Some current measured
            on floating tip.""")

E.rep(28884, 28883, "Repeat, with single tip short circuited", 
        descr = "Same current as on last shot")

E.rep(28885, 28880, "Repeat of 28880, with 10 s acquisition", 
        descr = "Strong variations of V measurements again")

E.rep(28886, 28885, "DAQ out of rack, only V measurement via divider", 
        descr = """\
            No voltage, batteries were disconnected. Some of them destroyed.""")

E.rep(28887, 28886, "No voltage connected", 
        descr = "")

E.rep(28888, 28887, "-44 V measured on battery before shot", 
        descr = "")

E.rep(28889, 28888, "-200 V via Kepco", 
        descr = "")

E.rep(28890, 28889, "Current measurement on S2", 
        descr = "")

E.rep(28891, 28890, "-200 V via Kepco on S1, batteries on S6", 
        descr = "")


############################################
E = campaign.add_experiment(date="20121123")

# DC biasing testing - NO POSITION SIGNAL
E.add(28894, "Continue from yesterday: Only V signal again.",
        head = headI,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            Voltage very noisy but stable.""",
        stars = '', **def_XPR)

E.rep(28895, 28894, "Put current measurement back on", 
        descr = "Voltage looks the same as on last shot.")

E.rep(28896, 28895, "With batteries on S6, like 28891", 
        descr = "Voltage lower than on last shot.")

E.rep(28897, 28896, "Repeat", 
        descr = "")

E.rep(28898, 28897, "Ground from Kepco to digitizer", 
        descr = "")

E.rep(28899, 28898, "Connect Kepco's minus to ground (green), with current on S2", 
        descr = "")

E.rep(28900, 28899, "Digitizer now on rack transformer, plunge at 1 s",
        times = 1.,
        posit = 0.05,
        descr = "Not armed in time. Measured -275 V, low noise")

E.rep(28901, 28898, "Digitizer back on 2nd transformer, no plunge",
        descr = "Measured -200 V, a bit more noise.")

E.rep(28902, 28901, "Total rack on 2nd transformer.",
        descr = "Measured -200 V")

E.rep(28903, 28902, "Total rack on rack transformer.",
        descr = "Measured -200 V")

E.rep(28904, 28903, "Everything on 2nd transformer, Kepco -200 V w/ ground connection",
        times = 1.,
        posit = 0.05,
        descr = "Worked! Quite a long arc, but Voltage didn't collapse.")


############################################
E = campaign.add_experiment(date="20121127")

# Reversed IpBt
E.add(28911, "Single/double Kepco on single/Mach tip, no plasma",
        head = headI_tip3sep,
        ampV1 = ampVF,
        ampV2 = ampVF,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            All signals are there.""",
        stars = '', **def_XPR_pos)

E.rep(28912, 28911, "Voltage signals 1 on S1 and S3",
        descr = """\
            No Voltage measurement, no currents on Mach tips.
            Capacitive pickup on single tip""",
        stars = '')

E.add(28913, "All on swept double Kepco, V on S1 and S3, multimeter @0.01 Hz: -225 to +57 V",
        times = 1.,
        posit = 0.02,
        head = headI,
        ampV1 = ampVF,
        ampV2 = ampVF,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            All signals on digitizer jumped down by 0.7 V.
            Position signal no longer noisy.""",
        stars = '*', **def_XPR_pos)

E.rep(28914, 28913, "Plunge at 5 cm",
        posit = 0.05,
        descr = """\
            Good data. Zero levels on ADC as on last shot.""",
        stars = '**')

E.rep(28915, 28914, "Zero level test: Lemo2BNC adapter on S3 (ground diff input).",
        descr = """\
            Strong pickup on S3""",
        stars = '**')

E.rep(28916, 28915, "Lemo2BNC adapter on S8. Go to 10 cm",
        posit = 0.1,
        descr = """\
            Nice signals. Noise is now on S8, as expected.
            Nice demonstration that sum of Mach tips is equivalent to single tip.""",
        stars = '****')

E.rep(28917, 28916, "Isolated tip at -200 V (1A Kepco)",
        head = headI_tip3sep,
        descr = """\
            Good data, but 35 kHz pickup on single tip.
            Good electron branch on swept Mach tips.""",
        stars = '***')

E.rep(28918, 28917, "Both Kepco's at -200 V. Plasma died before plunge.",
        descr = """\
            Again 35 kHz on single tip.""",
        stars = '')

E.rep(28919, 28918, "All tips at -200 V (4A Kepco)",
        head = headI,
        descr = """\
            35 kHz went away. All signals DC and very good.""",
        stars = '****')

E.rep(28920, 28919, "Mach on DC (4A), single on AC (1A), plunge to 15 cm",
        head = headI_tip3sep,
        posit = 0.15,
        descr = """\
            Bad arc.""",
        stars = '*')

E.rep(28921, 28920, "Plunge to 10 cm. No TS06",
        posit = 0.10,
        stars = '')

E.rep(28922, 28921, "Repeat",
        descr = """\
            No arcs. Nice comparison between swept and DC biased tips.""",
        stars = '***')

E.rep(28923, 28922, "Repeat",
        descr = """\
            Similar to last shot.""",
        stars = '***')

E.rep(28924, 28923, "Repeat",
        descr = """\
            Again similar""",
        stars = '***')

E.rep(28925, 28924, "Add plunge at 8 s to 34 cm.",
        times = (1., 8.),
        posit = (0.1, 0.34),
        descr = """\
            Again similar. Plunge at 8 s reveals that PosH signal is 
            NOT necessary to resolve maximum position.""",
        stars = '***')

E.rep(28926, 28925, "One plunge to 10 cm again",
        times = 1.,
        posit = 0.1,
        descr = """\
            OK. Probably not enough el current.""",
        stars = '***')

E.rep(28927, 28926, "15 cm at 0.9 s, sweep ampl from 14.5 to 15.5 Vpp",
        times = 0.9,
        posit = 0.15,
        descr = """\
            OK.""",
        stars = '***')

# Rachael McDermott
E.rep(28928, 28927, "32 cm at 4.4 s, shift sweeps positively out of saturation",
        times = 4.4,
        posit = 0.32,
        descr = """\
            Very nice Mach fluctuation measurements all the way across!
            Too much positive Vbias. Saturation in current leads to asymmetry 
            in Isat phase at very low densities.""",
        stars = '***')

E.rep(28929, 28928, "Back to 14.5 Vpp, shifted out of saturation. Backwards plunge at 1 s",
        times = (1., 4.4),
        posit = (-0.1, 0.32),
        descr = """\
            Plasma didn't run until 4.4 s.""",
        stars = '')

E.rep(28930, 28929, "Go to 1.9 s",
        times = 1.9,
        posit = 0.32,
        descr = """\
            Nice fluctuation data. Better sweep settings.""",
        stars = '****')




############################################
E = campaign.add_experiment(date="20121130")

E.add(28960, "Go to 10 cm at 0.9 s. Mach on DC, single on AC",
        times = 0.9,
        posit = 0.1,
        head = headI_tip3sep,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            Worked. Digitizer offset now different again!!!""",
        stars = '***', **def_XPR_pos)

E.rep(28961, 28960, "Go to 12 cm",
        posit = 0.12,
        descr = """\
            Worked.""",
        stars = '***')

E.rep(28962, 28961, "Go to 15 cm",
        posit = 0.15,
        descr = """\
            Arc at maximum position.""",
        stars = '**')

# Rachael McDermott
E.rep(28963, 28962, "32 cm at 1.9 and 3.9 s",
        times = (1.9, 3.9),
        posit = (0.32, 0.32),
        descr = """\
            Arcs on both plunges.""",
        stars = '**')

E.rep(28964, 28963, "Sweeps also on 4A Kepco.",
        times = (1.9, 3.9),
        posit = (0.32, 0.32),
        descr = """\
            Only small arcs on both plunges.""",
        stars = '****')

E.rep(28965, 28964, "4A Kepco -200 V again, plunge at 2 s to 20 cm, expect H mode",
        times = 2.,
        posit = 0.20,
        descr = """\
            3 MW of ECRH power. Plunge already in H-mode. Lots of arcing.""",
        stars = '*')

E.rep(28966, 28965, "Move plunge to 1.85 s",
        times = 1.85,
        posit = 0.20,
        descr = """\
            Good data in L-mode. Arcs at L-H transition on way out.""",
        stars = '**')

E.rep(28967, 28966, "Repeat",
        descr = """\
            Similar to last shot.""",
        stars = '**')

E.rep(28968, 28967, "Go to 18 cm",
        posit = 0.18,
        descr = """\
            Mach tip hit by beam blip very early.
            Swept single tip seems OK for the whole plunge.""",
        stars = '**')

# Bt calibration
E.rep(28973, 28968, "Digitizer offset calibration",
        times = (),
        posit = (),
        descr = """\
            All Kepcos off,""",
        stars = '')

# Felix Reimold
E.rep(28974, 28968, "All tips on swept 4A Kepco, 15 cm at 0.6 s, bias 17 Vpp",
        times = 0.6,
        posit = 0.15,
        head = headI,
        descr = """\
            Arcs again. Too much positive bias.""",
        stars = '*')

E.rep(28975, 28974, "Bias 14.5 Vpp, at 0.5 s",
        times = 0.5,
        posit = 0.15,
        descr = """\
            Finally no arcs, but nothing interesting.""",
        stars = '**')



############################################
E = campaign.add_experiment(date="20121204")

# Rev IpBt
E.add(28983, "Bias Voltage off, offset calibration",
        head = headI,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = "",
        stars = '', **def_XPR_pos)

# Steffen Potzel
E.rep(28986, 28983, "all tips on 4A Kepco, sweep 14.5 Vpp, at 1.5 and 2.5 s",
        times = (1.0, 1.5),
        posit = (0.34, 0.34),
        descr = """\
            Accidentally programmed 1.0 and 1.5 s. Plunges too close together, VPE.
            Data OK on 1st plunge, but not enough electron current.""",
        stars = '***')

# Gregor Birkenmeier
E.rep(28987, 28986, "sweep 15.5 Vpp, at 1.0 and 3.0 s",
        times = (1.0, 2.0),
        posit = (0.2, 0.2),
        descr = """\
            2nd plunge was at 2.0 s. Good data on first plunge. Single tip
            emitting at end of first plunge, strong emission/arc on second.""",
        stars = '***')

E.rep(28988, 28987, "Single tip on 1A Kepco (save Mach signals)",
        times = (1.0, 3.0),
        posit = (0.2, 0.12),
        head = headI_tip3sep,
        descr = """\
            Caught L-H transition on 2nd plunge!""",
        stars = '****')

# Steffen Potzel
E.rep(28989, 28988, "All the way through at 1.0 and 2.5 s",
        times = (1.0, 2.5),
        posit = (0.34, 0.34),
        descr = """\
            Great data on both plunges!""",
        stars = '*****')


############################################
E = campaign.add_experiment(date="20121205")

# ICRF wall conditioning
E.add(28998, "Two plunges at 1.0 and 2.5 s to 20 cm",
        times = (1.0, 2.5),
        posit = (0.2, 0.2),
        head = headI_tip3sep,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = "Plunged, but no signals.",
        stars = '', **def_XPR_pos)

E.rep(28999, 28998, "Repeat",
        descr = "No data again.",
        stars = '')

E.rep(29002, 28999, "Repeat",
        descr = "Measured something, but bias voltage affected by RF",
        stars = '')

E.rep(29003, 29002, "Repeat",
        descr = "Bias Voltage now strongly affected by RF.",
        stars = '')

E.rep(29005, 29003, "No plunge, check that Kepcos are still working.",
        descr = "OK.",
        stars = '')





