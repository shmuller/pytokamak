####################
# DATA PRIOR TO 2012
####################

from config_xpr_common import *
from config import Campaign

campaign = Campaign()

class TipLPS(CylindricalTip):
    def __init__(self, *args, **kw):
        CylindricalTip.__init__(self, 0.00045, 0.002, *args, **kw)

tip1 = TipLPS(number=1, pos='lower left', V_keys='ampVF', I_keys='ampI2')
tip2 = TipLPS(number=2, pos='lower right', V_keys='ampVF', I_keys='ampI3')
tip3 = TipLPS(number=3, pos='upper', V_keys='ampVF', I_keys='ampI1')

head = HeadXPR(tips=(tip1, tip2, tip3), R_keys='ampR')

amp_LPS_old = dict(ampR = Amp(fact=0.004, offs=-2745*0.004))

mapping_LPS_old = dict(
            ampR  = 'XPOS', 
            ampV1 = 'VOL1',
            ampI1 = 'CUR1',
            ampI2 = 'CUR2',
            ampI3 = 'VOL3',
            ampVF = 'VOL2')

lines_LPS = dict(amp=amp_LPS_old, mapping=mapping_LPS_old)

def_LPS_old = dict(dig='LPS_old', amp_default=amp_default_unity, lines=dict(LPS=lines_LPS))


############################################
E = campaign.add_experiment(date="20050621")

E.add(20326, "Standard Ohmic, Maximos fitting demo",
        head = head,
        descr = """\
            Nice shot. Sinusoidal sweep. Isat slanting visible.""",
        stars = '****', **def_LPS_old)

E.rep(20335, 20326, "LSN with upshifted X-point", 
        times = 2.5,
        posit = 0.31,
        descr = "Locked mode. Some sort of transition visible", 
        stars = '**')

E.rep(20336, 20335, "Repeat", 
        descr = "Similar to last shot, maybe a bit better.", 
        stars = '***')

E.rep(20337, 20336, "Repeat", 
        descr = "Large oscillations, poor fits.", 
        stars = '**')

E.rep(20339, 20337, "Modify strike point position", 
        descr = "H-mode, arcs.", 
        stars = '**')

E.rep(20341, 20339, "Repeat with interrupted DC biasing", 
        descr = "H-mode with ELMs. Arcs during ELMs and on the way out.", 
        stars = '**')

E.rep(20342, 20341, "Repeat", 
        descr = "Again pretty bad.", 
        stars = '**')

E.rep(20343, 20342, "Repeat", 
        descr = "Again similar.", 
        stars = '**')


############################################
E = campaign.add_experiment(date="20060404")

E.add(21256, "Maximos high density data", 
        times = 3.5, 
        posit = 0.31,
        head = head,
        descr = "Shot didn't run.",
        stars = '', **def_LPS_old)

E.rep(21258, 21256, "Plunge all the way through",
        descr = "",
        stars = '')


############################################
E = campaign.add_experiment(date="20060407")

E.add(21288, "Standard Ohmic", 
        times = 2.1, 
        posit = 0.34,
        head = head,
        descr = "",
        stars = '*', **def_LPS_old)

E.rep(21303, 21288, "Andrea's shot", 
        times = 2.1, 
        posit = 0.34,
        descr = "Very nice data. Mach numbers exactly 1 on HFS", 
        stars = '*****')

E.rep(21305, 21303, "Plunge again",
        descr = "Exactly the same as last shot",
        stars = '*****')

E.rep(21306, 21305, "Higher density",
        descr = "Higher density, Mach ~ 1.5 on HFS",
        stars = '*****')


############################################
E = campaign.add_experiment(date="20060411")

E.add(21320, "", 
        times = 2.2,
        posit = 0.31,
        head = head,
        descr = """\
            Data quality OK. Mach signals get very high, while 
            single tip signal is very low.""",
        stars = '***', **def_LPS_old)

E.rep(21321, 21320, "Higher density", 
        descr = "",
        stars = '***')

E.rep(21322, 21321, "", 
        descr = "Some sort of transition on way out.",
        stars = '**')

E.rep(21325, 21322, "", 
        descr = "Signals saturated on way out.",
        stars = '**')

E.rep(21326, 21325, "", 
        descr = "Similar to last shot.",
        stars = '***')

E.rep(21327, 21326, "No plunge", 
        stars = '')





