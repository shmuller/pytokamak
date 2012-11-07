######################
# 2012 AUTUMN CAMPAIGN
######################

tip1 = TipXPR(number=1, pos='lower left', V_keys='ampV', I_keys='ampI3')
tip2 = TipXPR(number=2, pos='lower right', V_keys='ampV', I_keys='ampI1')
tip3 = TipXPR(number=3, pos='upper', V_keys='ampV', I_keys='ampI2')

headI = Head(tips=(tip1, tip2, tip3), R_keys='ampR')


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






