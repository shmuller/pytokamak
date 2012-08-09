import numpy as np

import scipy.interpolate as interp

if __name__ == "__main__":
    import matplotlib
    #matplotlib.use('TkAgg')
    matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as plt

import probe
#reload(probe)

import sig
#reload(sig)


IOMds = probe.IOMds
IOFile = probe.IOFile
Digitizer = probe.Digitizer
Amp = probe.Amp
Probe = probe.Probe

ampUnity = Amp(fact=1., offs=0.)
ampInv   = Amp(fact=-1., offs=0.)
amp12Bit = Amp(fact=10./4095, offs=-5.)
amp14Bit = Amp(fact=20./16383, offs=-10.)


class IOMdsAUG(IOMds):
    def __init__(self, *args, **kw):
        # augsignal(_shot, _diag, _signame, _experiment, _edition, 
        #   _t1, _t2, _oshot, _oedition, _qual)

        diag = kw.pop('diag', 'XPR')

        IOMds.__init__(self, *args, **kw)
        self.mdsport = "8001"
        self.mdsfmt = '_s = augsignal(%%d,"%s","%%s","AUGD",*,*,*,*,*,"raw")' % diag

        self.datadeco = '%s; word_unsigned(data(_s))'
        self.timedeco = '%s; dim_of(_s)'
        self.sizedeco = '%s; size(_s)'


class IOFileAUG(IOFile):
    def __init__(self, shn, diag='XPR'):
        IOFile.__init__(self, shn=shn, suffix="_"+diag, subdir="AUG")


class DigitizerXPR(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock)

        self.IO_mds = IOMdsAUG(shn, sock, diag='XPR')
        self.IO_file = IOFileAUG(shn, diag='XPR')
        self.nodes = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 't')

        self.amp = {node: amp14Bit for node in self.nodes}
        self.amp['t'] = ampUnity


class DigitizerLPS(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock)

        self.IO_mds = IOMdsAUG(shn, sock, diag='LPS')
        self.IO_file = IOFileAUG(shn, diag='LPS')
        self.nodes = ('CUR1', 'VOL1', 'CUR2', 'VOL2', 'VOL3', 'VOL4', 't')

        self.amp = {node: amp12Bit for node in self.nodes}
        self.amp['t'] = ampUnity

        for node in ('CUR1', 'VOL3'):
            self.amp[node] = ampInv*amp12Bit


class ProbeXPoint(Probe):
    def __init__(self, digitizer=None):
        Probe.__init__(self, digitizer)

        self.window = slice(None)
        self.b2V = self.mapping = self.amp = None

    def mapsig(self):
        x = self.x[self.window]

        t  = x['t']
        R  = x[self.mapping['R']].astype('d')
        V  = x[self.mapping['V']].astype('d')
        I1 = x[self.mapping['I1']].astype('d')
        I2 = x[self.mapping['I2']].astype('d')
        VF = x[self.mapping['VF']].astype('d')

        R = sig.PositionSignal(R, t, name='R')
        V = sig.VoltageSignal(V, t, name='V')

        self.S = {'R': R,
                  'V': V,
                  'I1': sig.CurrentSignal(I1, t, V, name='I1'),
                  'I2': sig.CurrentSignal(I2, t, V, name='I2'),
                  'VF': sig.VoltageSignal(VF, t, name='VF')}

    def calib(self):
        R, V, I1, I2, VF = self['R', 'V', 'I1', 'I2', 'VF']

        R  *= self.amp['R']
        V  *= self.amp['V']
        I1 *= self.amp['I1']
        I2 *= self.amp['I2']
        VF *= self.amp['VF']

        s = slice(-5000, None)
        I1.norm_to_region(s)
        I2.norm_to_region(s)
        VF.norm_to_region(s)

        self.S['It'] = sig.CurrentSignal(I1.x + I2.x, I1.t, V, name='It')

    def position_calib(self):
        R = self['R']
        sp = interp.UnivariateSpline(R.t, R.x, s=10)
        Rs = sp(R.t)
        return Rs

    def voltage_calib(self):
        self.load(trim=False, calib=False)
        return self['V'].x.mean()


class ProbeXPR(ProbeXPoint):
    def __init__(self, shn, sock=None):
        digitizer = DigitizerXPR(shn, sock)

        ProbeXPoint.__init__(self, digitizer)

        self.window = slice(None)
        self.mapping = dict(R='S5', V='S1', I1='S4', I2='S2', VF='S6')
        
        fixpoints = (3.6812, -72), (7.0382, 170)
        ampR = Amp(fixpoints=fixpoints)
        ampV = Amp(fact=100., offs=-69.2227)

        ampI1 = Amp(fact=0.5*20./10)  # mA/mV = A/V (0.5 from missing 50 Ohm term.)
        ampI2 = Amp(fact=0.5*20./10)
        ampVF = Amp(fact=100.)

        self.amp = dict(R=ampR, V=ampV, I1=ampI1, I2=ampI2, VF=ampVF)


class ProbeLPS(ProbeXPoint):
    def __init__(self, shn, sock=None):
        digitizer = DigitizerLPS(shn, sock)

        ProbeXPoint.__init__(self, digitizer)

        self.window = slice(2048, None)
        self.mapping = dict(R='VOL3', V='VOL1', I1='CUR1', I2='CUR2', VF='VOL2')

        fixpoints = (-1.8767, -106), (3.8011, 336)
        ampR = Amp(fixpoints=fixpoints)
        ampV = Amp(fact=100., offs=-183.76)

        ampI1 = Amp(fact=0.5*20./10)  # mA/mV = A/V (0.5 from missing 50 Ohm term.)
        ampI2 = Amp(fact=0.5*20./10)
        ampVF = Amp(fact=100.)

        amp1x5 = Amp(fact=2.58)
        amp2x5 = Amp(fact=4.84)

        ampI1 *= amp1x5.inv()
        ampI2 *= amp2x5.inv()

        if shn < 28426:
            ampI1 *= ampInv

        if shn < 27687:
            ampI2 *= ampInv

        self.amp = dict(R=ampR, V=ampV, I1=ampI1, I2=ampI2, VF=ampVF)


if __name__ == "__main__":
    shn = 28469

    XPR = ProbeXPR(shn=shn)
    XPR.analyze(plunge=0)

    XPR.plot()
    plt.show()



