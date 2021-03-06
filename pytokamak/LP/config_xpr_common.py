import numpy as np

from config import CylindricalTip, Head

from pytokamak.utils.utils import rdict
from pytokamak.utils.sig import Amp

ampUnity = Amp(fact=1.)
ampInv   = Amp(fact=-1.)

ampVF = Amp(fact=100.)


Preamp1 = {
     2: Amp(fact=1.032).inv(), 
     5: Amp(fact=2.580).inv()}

Preamp2 = {
     2: Amp(fact=1.936).inv(), 
     5: Amp(fact=4.840).inv()}

'''
The current probes are set in mA/div. They source 10 mV/div into a 50 Ohm 
input impedance, up to a maximum of 1 V. Since the DAQ's impedance is high-Z, 
this is effectively 20 mV/div or a maximum of 2 V. The calibration factor is
thus:

fact (A/V) = settings (mA/div) / 20

At 20 mA/div, the calibration is 1-to-1.
'''
CurrentProbe1 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10), 
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10),
   100: Amp(fact=0.5*100/10),
  5000: Amp(fact=0.5*5000/10)}

CurrentProbe2 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10),
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10),
   100: Amp(fact=0.5*100/10),
  5000: Amp(fact=0.5*5000/10)}

CurrentProbe3 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10),
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10),
   100: Amp(fact=0.5*100/10),
  5000: Amp(fact=0.5*5000/10)}

CurrentProbe4 = {
     1: Amp(fact=0.5*1/10), # mA/mV = A/V (0.5 from missing 50 Ohm term.)
     5: Amp(fact=0.5*5/10),
    10: Amp(fact=0.5*10/10),
    20: Amp(fact=0.5*20/10),
    50: Amp(fact=0.5*50/10),
   100: Amp(fact=0.5*100/10),
  5000: Amp(fact=0.5*5000/10)}


class TipXPR(CylindricalTip):
    def __init__(self, r=0.0005, z=0.003, *args, **kw):
        CylindricalTip.__init__(self, r=r, z=z, *args, **kw)


class TipXPR1(TipXPR):
    def __init__(self, *args, **kw):
        TipXPR.__init__(self, number=1, pos='lower left', label='Mach tip 1', *args, **kw)

class TipXPR2(TipXPR):
    def __init__(self, *args, **kw):
        TipXPR.__init__(self, number=2, pos='lower right', label='Mach tip 2', *args, **kw)

class TipXPR3(TipXPR):
    def __init__(self, *args, **kw):
        TipXPR.__init__(self, number=3, pos='upper', label='Single tip', *args, **kw)

class TipXPR4(TipXPR):
    def __init__(self, *args, **kw):
        TipXPR.__init__(self, number=4, pos='lower', label='Mach tips sum (meas.)', *args, **kw)


class HeadXPR(Head):
    def __init__(self, tips, d=0.02):
        Head.__init__(self, tips, d)


amp_default_unity = dict(
            ampR  = ampUnity, 
            ampV1 = ampUnity, 
            ampV2 = ampUnity,
            ampV3 = ampUnity,
            ampI1 = ampUnity,
            ampI2 = ampUnity,
            ampI3 = ampUnity,
            ampI4 = ampUnity,
            ampVF = ampUnity)

amp_default = dict(
            ampR  = None, 
            ampV1 = None, 
            ampV2 = None,
            ampV3 = None,
            ampI1 = ampUnity,
            ampI2 = ampUnity,
            ampI3 = ampUnity,
            ampI4 = ampUnity,
            ampVF = ampVF)

amp_XPR = dict(
            ampR = None, 
            ampV1 = Amp(fact=100., offs=-69.2227))

amp_LPS = dict(
            ampR = None, 
            ampV1 = Amp(fact=100., offs=-183.76))

mapping_XPR = dict(
            ampR  = 'S5', 
            ampV1 = 'S1', 
            ampV2 = 'S3',
            ampI1 = 'S4',
            ampI2 = 'S2',
            ampI3 = 'S6',
            ampVF = 'S6')

mapping_LPS = dict(
            ampR  = 'VOL3', 
            ampV1 = 'VOL1',
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


