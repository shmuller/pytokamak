####################
# DATA PRIOR TO 2012
####################

from config_xpr_common import *
from config import Campaign

campaign = Campaign()

fixpoints = (-1.8767, -0.106), (3.8011, 0.336)
amp_LPS['ampR'] = Amp(fixpoints=fixpoints)

tip1 = TipXPR(number=1, pos='lower left', V_keys='ampV', I_keys='ampI2')
tip2 = TipXPR(number=2, pos='lower right', V_keys='ampV', I_keys='ampI3')
tip3 = TipXPR(number=3, pos='upper', V_keys='ampVF', I_keys='ampI1')

headI = Head(tips=(tip1, tip2, tip3), R_keys='ampR')

mapping_LPS = dict(
            ampR  = 'VOL4', 
            ampV  = 'VOL1',
            ampI1 = 'CUR1',
            ampI2 = 'CUR2',
            ampI3 = 'VOL3',
            ampVF = 'VOL2')

lines_LPS = dict(amp=dict(), mapping=mapping_LPS)

def_LPS = dict(dig='LPS_old', amp_default=amp_default_unity, lines=dict(LPS=lines_LPS))


############################################
E = campaign.add_experiment(date="20090101")

E.add(21303, "Andrea's shot", 
        times = 1.2, 
        posit = 0.10,
        head = headI,
        descr = "Very nice data", 
        stars = '*****', **def_LPS)

E.rep(21304, 21303, "No plunge")

E.rep(21305, 21304, "Plunge again")

