import numpy as np

import scipy.interpolate as interp

if __name__ == "__main__":
    import matplotlib
    #matplotlib.use('TkAgg')
    matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as plt

import probe
#reload(probe)

import config
#reload(config)

IOMds = probe.IOMds
IOFile = probe.IOFile
Digitizer = probe.Digitizer
Amp = probe.Amp
Probe = probe.Probe

PositionSignal = probe.PositionSignal
VoltageSignal = probe.VoltageSignal
CurrentSignal = probe.CurrentSignal

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

        #f = 20./16283
        #offs = (8759, 8830, 8756, 8770, 8759, 8747, 8746, 8752, 0)
        #self.amp = {node: Amp(fact=f, offs=-f*o) for node, o in zip(self.nodes, offs)}
        self.amp = {node: amp14Bit.copy() for node in self.nodes}
        self.amp['t'] = ampUnity.copy()


class DigitizerLPS(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock)

        self.IO_mds = IOMdsAUG(shn, sock, diag='LPS')
        self.IO_file = IOFileAUG(shn, diag='LPS')
        self.nodes = ('CUR1', 'VOL1', 'CUR2', 'VOL2', 'VOL3', 'VOL4', 't')

        #f = 10./4095
        #offs = (0,0,0,0,0,0,0)
        #self.amp = {node: Amp(fact=f, offs=-f*o) for node, o in zip(self.nodes, offs)}
        self.amp = {node: amp12Bit.copy() for node in self.nodes}
        self.amp['t'] = ampUnity.copy()

        for node in ('CUR1', 'VOL3'):
            self.amp[node] *= ampInv


class ProbeXPoint(Probe):
    def __init__(self, digitizer=None):
        Probe.__init__(self, digitizer)
        self.config = None

    def mapsig(self):
        window, mapping = self.config.window, self.config.mapping

        x = self.x[window]

        t  = x['t']
        R  = x[mapping['R']].astype('d')
        V  = x[mapping['V']].astype('d')
        I1 = x[mapping['I1']].astype('d')
        I2 = x[mapping['I2']].astype('d')
        VF = x[mapping['VF']].astype('d')

        R = PositionSignal(R, t, name='R')
        V = VoltageSignal(V, t, name='V')

        self.S = {'R': R,
                  'V': V,
                  'I1': CurrentSignal(I1, t, V, name='I1'),
                  'I2': CurrentSignal(I2, t, V, name='I2'),
                  'VF': VoltageSignal(VF, t, name='VF')}

    def calib(self):
        amp = self.config.amp

        R, V, I1, I2, VF = self['R', 'V', 'I1', 'I2', 'VF']

        R  *= amp['R']
        V  *= amp['V']
        I1 *= amp['I1']
        I2 *= amp['I2']
        VF *= amp['VF']

        s = slice(-5000, None)
        I1.norm_to_region(s)
        I2.norm_to_region(s)
        VF.norm_to_region(s)

        self.S['It'] = CurrentSignal(I1.x + I2.x, I1.t, V, name='It')

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

        self.config = config.campaign.find_shot(shn)

        """
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
        """

if __name__ == "__main__":
    shn = 28469

    XPR = ProbeXPR(shn=shn)
    XPR.analyze(plunge=0)

    XPR.plot()
    plt.show()



