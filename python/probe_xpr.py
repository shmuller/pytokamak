import numpy as np

import scipy.interpolate as interp

if __name__ == "__main__":
    import matplotlib
    #matplotlib.use('TkAgg')
    matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as plt

import probe
#reload(probe)


class IOMdsAUG(probe.IOMds):
    def __init__(self, *args, **kw):
        # augsignal(_shot, _diag, _signame, _experiment, _edition, 
        #   _t1, _t2, _oshot, _oedition, _qual)

        diag = kw.pop('diag', 'XPR')

        probe.IOMds.__init__(self, *args, **kw)
        self.mdsport = "8001"
        self.mdsfmt = '_s = augsignal(%%d,"%s","%%s","AUGD",*,*,*,*,*,"raw")' % diag

        self.datadeco = '%s; word_unsigned(data(_s))'
        self.timedeco = '%s; dim_of(_s)'
        self.sizedeco = '%s; size(_s)'


class IOFileAUG(probe.IOFile):
    def __init__(self, shn, diag='XPR'):
        probe.IOFile.__init__(self, shn=shn, suffix="_"+diag, subdir="AUG")


class ProbeXPoint(probe.Probe):
    def __init__(self, shn, sock=None, diag='XPR', nodes=()):
        probe.Probe.__init__(self, shn, sock)

        self.IO_mds = IOMdsAUG(shn, sock, diag=diag)
        self.IO_file = IOFileAUG(shn, diag=diag)
        self.nodes = nodes

        self.window = slice(None)
        self.b2V = self.mapping = self.amp = None

    def mapsig(self):
        x = self.x
        s = self.window

        t  = x['t'][s]
        R  = self.b2V['R']*x[self.mapping['R']][s].astype('d')
        V  = self.b2V['V']*x[self.mapping['V']][s].astype('d')
        I1 = self.b2V['I1']*x[self.mapping['I1']][s].astype('d')
        I2 = self.b2V['I2']*x[self.mapping['I2']][s].astype('d')
        VF = self.b2V['VF']*x[self.mapping['VF']][s].astype('d')

        R = probe.PositionSignal(t, R, name='R')
        V = probe.VoltageSignal(t, V, name='V')

        self.S = {'R': R,
                  'V': V,
                  'I1': probe.CurrentSignal(t, I1, V, name='I1'),
                  'I2': probe.CurrentSignal(t, I2, V, name='I2'),
                  'VF': probe.VoltageSignal(t, VF, name='VF')}

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

        self.S['It'] = probe.CurrentSignal(I1.t, I1.x + I2.x, V, name='It')

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
        nodes = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 't')

        ProbeXPoint.__init__(self, shn, sock, diag='XPR', nodes=nodes)

        self.window = slice(None)
        self.mapping = dict(R='S5', V='S1', I1='S4', I2='S2', VF='S6')
        
        b2V = probe.Amp(fact=20./16383, offs=-10.)
        self.b2V = dict(R=b2V, V=b2V, I1=b2V, I2=b2V, VF=b2V)

        fixpoints = (3.6812, -72), (7.0382, 170)
        ampR = probe.Amp(fixpoints=fixpoints)
        ampV = probe.Amp(fact=100., offs=-69.2227)

        ampI1 = probe.Amp(fact=0.5*20./10)  # mA/mV = A/V (0.5 from missing 50 Ohm term.)
        ampI2 = probe.Amp(fact=0.5*20./10)
        ampVF = probe.Amp(fact=100.)

        self.amp = dict(R=ampR, V=ampV, I1=ampI1, I2=ampI2, VF=ampVF)


class ProbeLPS(ProbeXPoint):
    def __init__(self, shn, sock=None):
        nodes = ('CUR1', 'VOL1', 'CUR2', 'VOL2', 'VOL3', 'VOL4', 't')

        ProbeXPoint.__init__(self, shn, sock, diag='LPS', nodes=nodes)

        self.window = slice(2048, None)
        self.mapping = dict(R='VOL3', V='VOL1', I1='CUR1', I2='CUR2', VF='VOL2')

        inv = probe.Amp(fact=-1.)
        b2V = probe.Amp(fact=10./4095, offs=-5.)
        self.b2V = dict(R=inv*b2V, V=b2V, I1=inv*b2V, I2=b2V, VF=b2V)

        fixpoints = (-1.8767, -106), (3.8011, 336)
        ampR = probe.Amp(fixpoints=fixpoints)
        ampV = probe.Amp(fact=100., offs=-183.76)

        ampI1 = probe.Amp(fact=0.5*20./10)  # mA/mV = A/V (0.5 from missing 50 Ohm term.)
        ampI2 = probe.Amp(fact=0.5*20./10)
        ampVF = probe.Amp(fact=100.)

        amp1x5 = probe.Amp(fact=2.58)
        amp2x5 = probe.Amp(fact=4.84)

        ampI1 *= amp1x5.inv()
        ampI2 *= amp2x5.inv()

        if shn < 28426:
            ampI1 *= inv

        if shn < 27687:
            ampI2 *= inv

        self.amp = dict(R=ampR, V=ampV, I1=ampI1, I2=ampI2, VF=ampVF)


if __name__ == "__main__":
    shn = 28469

    XPR = ProbeXPR(shn=shn)
    XPR.analyze(plunge=0)

    XPR.plot()
    plt.show()



