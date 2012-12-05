####################
# DATA PRIOR TO 2012
####################

from config_xpr_common import *
from config import Campaign

campaign = Campaign()

tip1 = TipXPR(number=1, pos='lower left', V_keys='ampV1', I_keys='ampI2')
tip2 = TipXPR(number=2, pos='lower right', V_keys='ampV1', I_keys='ampI3')
tip3 = TipXPR(number=3, pos='upper', V_keys='ampVF', I_keys='ampI1')

headI = HeadXPR(tips=(tip1, tip2, tip3), R_keys='ampR')

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
E = campaign.add_experiment(date="20090101")

E.add(21288, "Standard Ohmic", 
        times = 2.1, 
        posit = 0.34,
        head = headI,
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

E.rep(21306, 21305, "Plunge again",
        descr = "Higher density, Mach ~ 1.5 on HFS",
        stars = '*****')



