######################
# 2013 SPRING CAMPAIGN
######################

from config_xpr_common import *
from config import Campaign

campaign = Campaign()

tipmap = rdict(
        tip1 = rdict(V='ampV1', I='ampI3'),
        tip2 = rdict(V='ampV1', I='ampI1'),
        tip3 = rdict(V='ampV2', I='ampI2'),
        tip4 = rdict(V='ampV1', I='ampI4'))

tipmap_tip1_V3 = tipmap.rep(tip1_V='ampV3')

tip1_20130306 = TipXPR1(r=0.0005, z=0.0026)
tip2_20130306 = TipXPR2(r=0.0005, z=0.0026)
tip3_20130306 = TipXPR3(r=0.0005, z=0.0020)
tip4_20130306 = TipXPR4(r=0.0005, z=0.0052)

head_20130306 = HeadXPR(tips=(tip1_20130306, tip2_20130306, tip3_20130306))
head_20130306_4tips = HeadXPR(tips=head_20130306.tips + (tip4_20130306,))


fact = 4 * 5.54630/27. / 2**16
offs = -47578968*fact - 0.105

amp_XPR_pos = dict(
            ampR  = Amp(fact=fact, offs=offs),
            ampV1 = Amp(fact=100., offs=-68.48),
            ampV2 = Amp(fact=100., offs=-68.13),
            ampV3 = Amp(fact=100., offs=-68.13))

mapping_XPR_pos = dict(
            ampR  = 'Pos', 
            ampV1 = 'S1',
            ampI2 = 'S2',
            ampV2 = 'S3',
            ampI1 = 'S4',
            ampV3 = 'S5',  # hook up V3 where analog position signal was before
            ampI3 = 'S6',
            ampI4 = 'S7')

lines_XPR_pos = dict(amp=amp_XPR_pos, mapping=mapping_XPR_pos)

def_XPR_pos = dict(dig='XPR_pos', amp_default=amp_default, lines=dict(XPR=lines_XPR_pos))



############################################
E = campaign.add_experiment(date="20140227")

# Startup 2014
E.add(30254, "DAQ test: Sweeps at 14.5 Vpp on all tips",
        head = head_20130306_4tips,
        tipmap = tipmap,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        ampI4 = CurrentProbe4[50],
        descr = """\
            DAQ working.""",
        stars = '', **def_XPR_pos)

E.rep(30255, 30254, "2 plunges behind wall",
        times = (2.0, 4.0),
        posit = (-0.05, -0.05),
        descr = """\
            No plunges. Rack with trigger light emitter was switched off.""",
        stars = '')

E.rep(30256, 30255, "2 plunges behind wall, next try",
        times = (2.0, 4.0),
        posit = (-0.05, -0.05),
        descr = """\
            Worked""",
        stars = '')

E.rep(30257, 30256, "1 plunge at 6 s, 5 cm, 13.5 Vpp, no arc box",
        times = 6.0,
        posit = 0.05,
        descr = """\
            Not enough acquisition window.""",
        stars = '')

E.rep(30258, 30257, "Try again, with 8 s acquisition window",
        times = 6.5,
        posit = 0.05,
        descr = """\
            Plasma short.""",
        stars = '')

E.rep(30259, 30258, "No plunge",
        times = (),
        posit = (),
        descr = """\
            """,
        stars = '')

E.rep(30260, 30259, "1 plunge at 6.5 s",
        times = 6.5,
        posit = 0.07,
        descr = """\
            Everything seems OK.""",
        stars = '**')

# Hans-Werner Mueller
E.add(30265, "No plunge - hook up V3 on S5 (old analog position signal), 12.5 Vpp",
        head = head_20130306_4tips,
        tipmap = tipmap_tip1_V3,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        ampI4 = CurrentProbe4[50],
        descr = """\
            """,
        stars = '', **def_XPR_pos)

E.rep(30266, 30265, "1 plunge at 4.65 s, to 22 cm",
        times = 4.65,
        posit = 0.22,
        descr = """\
            Plasma disrupted immediately before probe plunge. No magnetics
            data acquired, so repeat.""",
        stars = '')

E.rep(30267, 30266, "Repeat, go to 15 cm",
        times = 4.65,
        posit = 0.15,
        descr = """\
            No plasma.""",
        stars = '')

E.rep(30268, 30267, "Repeat",
        times = 4.65,
        posit = 0.15,
        descr = """\
            Disrupted before plunge again.""",
        stars = '')

E.rep(30269, 30268, "Repeat",
        times = 4.65,
        posit = 0.15,
        descr = """\
            Not enough electron current.""",
        stars = '')


E = campaign.add_experiment(date="20140305")

# Hans-Werner Mueller
E.add(30276, "14.0 Vpp on all tips",
        times = 4.65,
        posit = 0.15,
        head = head_20130306_4tips,
        tipmap = tipmap_tip1_V3,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        ampI4 = CurrentProbe4[50],
        descr = """\
            Didn't make it to probe plunge.""",
        stars = '', **def_XPR_pos)

E.rep(30277, 30276, "Repeat",
        descr = """\
            Didn't make it to probe plunge again.""",
        stars = '')

E.rep(30278, 30277, "2 plunges, all the way through",
        times = (2.65, 4.75),
        posit = (0.34, 0.34),
        descr = """\
            Good data on both plunges. Small arcs switched off by arc box.
            Now too much electron current.""",
        stars = '***')

E.rep(30279, 30278, "13.5 Vpp, 2 kHz sine",
        times = (2.65, 4.75),
        posit = (0.34, 0.34),
        descr = """\
            Didn't run. 2 kHz seems too high for capacitive pickup.""",
        stars = '')

E.rep(30280, 30279, "13.5 Vpp, 1 kHz sine",
        descr = """\
            Didn't run.""",
        stars = '')

E.rep(30281, 30280, "Repeat, 13.5 Vpp, 1 kHz sine",
        descr = """\
            OK data.""",
        stars = '**')

E.rep(30282, 30281, "Momentum transport, 13.0 Vpp, 1 kHz sine",
        times = 5.25,
        posit = 0.2,
        descr = """\
            Strong H-mode. Lots of arcing, but all measurements seem OK.""",
        stars = '*')

E.rep(30283, 30282, "NTM shot, back to 13.5 Vpp, no plunge",
        descr = "No plunge",
        stars = '')

# Leena Aho-Mantila
E.rep(30284, 30283, "One plunge at 4.3 s",
        times = 4.3,
        posit = 0.34,
        descr = """\
            H-mode. Many arcs again, but fits OK.""",
        stars = '***')

E.rep(30285, 30284, "One plunge at 4.65 s",
        times = 4.65,
        posit = 0.34,
        descr = """\
            Worked nicely. Only small arcs.""",
        stars = '****')

E.rep(30286, 30285, "One plunge at 4.60 s",
        times = 4.60,
        posit = 0.34,
        descr = """\
            A bit more arcs than before but still good.""",
        stars = '***')

E.rep(30287, 30286, "Now with N2",
        times = 4.60,
        posit = 0.34,
        descr = """\
            Good again. Significant difference in Vf and Te between
            way in and out. Beam blip on way out after PF region.""",
        stars = '****')

E.rep(30288, 30287, "More N2",
        times = 4.60,
        posit = 0.34,
        descr = """\
            Great data. Now even more difference.""",
        stars = '****')

E.rep(30289, 30288, "Repeat, due to feedback oscillations",
        times = 4.60,
        posit = 0.34,
        descr = """\
            Density rise (I-phase?) prior to probe plunge. Good data.""",
        stars = '****')

E.rep(30290, 30289, "Repeat again",
        times = 4.60,
        posit = 0.34,
        descr = """\
            Disrupted at 3.5 s.""",
        stars = '')

# Hans-Werner Mueller
E.rep(30291, 30282, "Repeat 30282",
        times = 5.25,
        posit = 0.1,
        descr = """\
            OK, but arcs on each of 4 ELMs.""",
        stars = '**')


E = campaign.add_experiment(date="20140306")

# Hans-Werner Mueller
E.add(30303, "13.5 Vpp on all tips",
        times = 4.35,
        posit = 0.22,
        head = head_20130306_4tips,
        tipmap = tipmap_tip1_V3,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        ampI4 = CurrentProbe4[50],
        descr = """\
            Disrupted before plunge.""",
        stars = '', **def_XPR_pos)

E.rep(30304, 30303, "Repeat of 30268",
        times = 4.5,
        posit = 0.22,
        descr = """\
            Disrupted really close to plunge.""",
        stars = '')

# Stefan Muller
E.rep(30305, 30304, "Detachment fluctuations",
        times = (2.0, 3.25, 4.5),
        posit = (0.34, 0.34, 0.34),
        descr = """\
            Three great plunges, straight through X-point. Only minor arcing,
            switched off by arc box. Unfortunately, fluctuations naturally
            disappeared before 1st and 2nd plunge, and were also absent while
            the probe was on the HFS in the 3rd plunge. Bad luck!
            X-point config: EOC 1.24. Thomas Eich saved config as *_XPProbe.
            """,
        stars = '****')


