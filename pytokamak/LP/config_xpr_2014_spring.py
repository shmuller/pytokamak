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
            ampV2 = 'S3',
            ampI1 = 'S4',
            ampI2 = 'S2',
            ampI3 = 'S6',
            ampVF = 'S6',
            ampI4 = 'S7',
            ampV3 = 'S8')

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



