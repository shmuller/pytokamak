import numpy as np

from config import CylindricalTip, Head, Campaign

from probe import Amp


ampUnity = Amp(fact=1.)
ampInv   = Amp(fact=-1.)

fixpoints = (-1.8767, -0.106), (3.8011, 0.336)
ampR_LPS  = Amp(fixpoints=fixpoints)
ampV_LPS  = Amp(fact=100., offs=-183.76)

fixpoints = (3.6812, -0.072), (7.0382, 0.170)
ampR_XPR = Amp(fixpoints=fixpoints)
ampV_XPR = Amp(fact=100., offs=-69.2227)

ampVF = Amp(fact=100.)


Preamp1 = {
     2: Amp(fact=1.032).inv(), 
     5: Amp(fact=2.580).inv()}

Preamp2 = {
     2: Amp(fact=1.936).inv(), 
     5: Amp(fact=4.840).inv()}

CurrentProbe1 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10), 
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10),
  5000: Amp(fact=0.5*5000/10)}

CurrentProbe2 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10),
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10),
  5000: Amp(fact=0.5*5000/10)}

CurrentProbe3 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10),
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10),
  5000: Amp(fact=0.5*5000/10)}


class TipXPR(CylindricalTip):
    def __init__(self, *args, **kw):
        CylindricalTip.__init__(self, 0.0005, 0.003, *args, **kw)


amp_default_unity = dict(
            ampR  = ampUnity, 
            ampV  = ampUnity, 
            ampI1 = ampUnity,
            ampI2 = ampUnity,
            ampI3 = ampUnity,
            ampVF = ampUnity)

amp_default = dict(
            ampR  = None, 
            ampV  = None, 
            ampI1 = CurrentProbe1[20],
            ampI2 = CurrentProbe2[20],
            ampI3 = CurrentProbe3[20],
            ampVF = ampVF)

amp_XPR = dict(
            ampR = ampR_XPR, 
            ampV = ampV_XPR)

amp_LPS = dict(
            ampR = ampR_LPS, 
            ampV = ampV_LPS)

mapping_XPR = dict(
            ampR  = 'S5', 
            ampV  = 'S1', 
            ampI1 = 'S4',
            ampI2 = 'S2',
            ampI3 = 'S6',
            ampVF = 'S6')

mapping_LPS = dict(
            ampR  = 'VOL3', 
            ampV  = 'VOL1',
            ampI1 = 'CUR1',
            ampI2 = 'CUR2',
            ampI3 = 'VOL2',
            ampVF = 'VOL2')

lines_XPR = dict(amp=amp_XPR, mapping=mapping_XPR)
lines_LPS = dict(amp=amp_LPS, mapping=mapping_LPS)

lines = dict(XPR=lines_XPR, LPS=lines_LPS)

def_LPS = dict(dig='LPS', amp_default=amp_default, lines=dict(LPS=lines_LPS))
def_XPR = dict(dig='XPR', amp_default=amp_default, lines=dict(XPR=lines_XPR))

def_XPR_LPS = dict(dig='XPR', amp_default=amp_default, lines=lines)

campaign = Campaign()


execfile("config_xpr_2012_spring.py")

execfile("config_xpr_2012_autumn.py")



