####################
# DATA PRIOR TO 2012
####################

from config_xpr_common import *
from config import Campaign

campaign = Campaign()

class TipLPS(CylindricalTip):
    def __init__(self, *args, **kw):
        CylindricalTip.__init__(self, 0.00045, 0.002, *args, **kw)


tip1 = TipLPS(number=1, pos='lower left', V_keys='ampV1', I_keys='ampI2')
tip2 = TipLPS(number=2, pos='lower right', V_keys='ampV1', I_keys='ampI3')
tip3 = TipLPS(number=3, pos='upper', V_keys='ampV1', I_keys='ampI1')

head = HeadXPR(tips=(tip1, tip2, tip3))

tipmap = dict(
        tip1 = dict(V='ampV1', I='ampI2'),
        tip2 = dict(V='ampV1', I='ampI3'),
        tip3 = dict(V='ampV1', I='ampI1'))


tip1_V2 = TipLPS(number=1, pos='lower left', V_keys='ampV2', I_keys='ampI2')
tip2_V2 = TipLPS(number=2, pos='lower right', V_keys='ampV2', I_keys='ampI3')
tip3_V2 = TipLPS(number=3, pos='upper', V_keys='ampV2', I_keys='ampI1')

head_V2 = HeadXPR(tips=(tip1_V2, tip2_V2, tip3_V2))

tipmap_V2 = dict(
        tip1 = dict(V='ampV2', I='ampI2'),
        tip2 = dict(V='ampV2', I='ampI3'),
        tip3 = dict(V='ampV2', I='ampI1'))


amp_LPS_old = dict(ampR = Amp(fact=0.004, offs=-2745*0.004))

mapping_LPS_old = dict(
        ampR  = 'XPOS', 
        ampV1 = 'VOL1',
        ampI1 = 'CUR1',
        ampI2 = 'CUR2',
        ampI3 = 'VOL3',
        ampV2 = 'VOL2')

lines_LPS = dict(amp=amp_LPS_old, mapping=mapping_LPS_old)

def_LPS_old = dict(dig='LPS_old', amp_default=amp_default_unity, lines=dict(LPS=lines_LPS))


############################################
E = campaign.add_experiment(date="19990323")

E.add(11816, "",
        times = 1.2,
        posit = 0.13,
        head = head,
        tipmap = tipmap,
        descr = "XXX mapping incorrect!",
        stars = '**', **def_LPS_old)

E.rep(11817, 11816, "",
        times = 1.5,
        posit = 0.31,
        descr = """\ 
            (2005_Tsalas_JNM)""",
        stars = '**')


############################################
E = campaign.add_experiment(date="20040318")

E.add(18786, "No plunge",
        head = head_V2,
        tipmap = tipmap_V2,
        descr = "XXX mapping incorrect!",
        stars = '', **def_LPS_old)

E.rep(18787, 18786, "",
        times = 4.9,
        posit = 0.18,
        descr = """\
            """,
        stars = '**')

E.rep(18788, 18787, "",
        times = 3.4,
        posit = 0.18,
        descr = "Didn't run.",
        stars = '')

E.rep(18789, 18788, "",
        times = 3.4,
        posit = 0.18,
        descr = """\
            Probe influencing ELMs! (2005_Tsalas_JNM)""",
        stars = '*****')


############################################
E = campaign.add_experiment(date="20050330")

E.add(19951, "Standard Ohmic, sweeps",
        times = 3.9,
        posit = 0.31,
        head = head_V2,
        tipmap = tipmap_V2,
        descr = """\
            Nice swept data, some arcs. (2007_Tsalas_PPCF)""",
        stars = '****', **def_LPS_old)

E.rep(19952, 19951, "Higher density", 
        times = 2.6,
        posit = 0.31,
        descr = "Good data too. (2007_Tsalas_PPCF)", 
        stars = '****')


############################################
E = campaign.add_experiment(date="20050621")

E.add(20326, "Standard Ohmic, Maximos fitting demo",
        head = head_V2,
        tipmap = tipmap_V2,
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
        descr = "H-mode with type I ELMs, arcs. (2007_Tsalas_PPCF)", 
        stars = '****')

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
E = campaign.add_experiment(date="20050624")

E.add(20372, "Standard Ohmic",
        times = 3.4,
        posit = 0.31,
        head = head_V2,
        tipmap = tipmap_V2,
        descr = """\
            All tips DC biased. Nice flow data from Mach -2 to Mach 2.
            No arcs. (2007_Tsalas_PPCF)""",
        stars = '****', **def_LPS_old)

E.rep(20373, 20372, "Standard H-mode",
        times = 2.1,
        posit = 0.13,
        descr = """\
            Rectangular sweep.""",
        stars = '***')


############################################
E = campaign.add_experiment(date="20060316")

E.add(21194, "Standard Ohmic",
        times = 1.9,
        posit = 0.31,
        head = head_V2,
        tipmap = tipmap_V2,
        descr = """\
            All tips swept. Very nice flow data from Mach -1 to Mach 1.
            No arcs. (2007_Tsalas_PPCF)""",
        stars = '*****', **def_LPS_old)


############################################
E = campaign.add_experiment(date="20060404")

E.add(21258, "Maximos high density data",
        times = 3.4,
        posit = 0.31,
        head = head_V2,
        tipmap = tipmap_V2,
        descr = "Excellent profiles, low HFS density. (2007_Tsalas_PPCF)",
        stars = '*****', **def_LPS_old)


############################################
E = campaign.add_experiment(date="20060407")

E.add(21288, "Standard Ohmic, Mach tips DC biased, single tip floating", 
        times = 1.9, 
        posit = 0.31,
        head = head,
        descr = """\
            Excellent fluctuation data. Compare with 21194 for swept profiles.""",
        stars = '*****', **def_LPS_old)

E.rep(21303, 21288, "Andrea's shot", 
        times = 2.0, 
        posit = 0.31,
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
        head = head_V2,
        tipmap = tipmap_V2,
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


############################################
E = campaign.add_experiment(date="20060421")

E.add(21380, "Standard Ohmic",
        times = 3.6,
        posit = 0.31,
        head = head_V2,
        tipmap = tipmap_V2,
        LPS_mapping_ampI1='CUR2',
        LPS_mapping_ampI2='CUR1',
        descr = """\
            Tip 1 and tip 3 reversed.
            All tips swept. Very nice flow data from Mach -1 to Mach 1.
            No arcs. (2007_Tsalas_PPCF)""",
        stars = '*****', **def_LPS_old)



