######################
# 2013 SPRING CAMPAIGN
######################

from config_xpr_common import *
from config import Campaign

campaign = Campaign()

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
E = campaign.add_experiment(date="20130124")

E.add(29289, "Probe test",
        head = headI_tip3sep,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        ampI3 = CurrentProbe3[20], 
        descr = """\
            """,
        stars = '', **def_XPR_pos)

E.rep(29290, 29289, "Probe test",
        times = (1.6, 3.9),
        posit = (0.17, 0.17),
        descr = """\
            """,
        stars = '***')

