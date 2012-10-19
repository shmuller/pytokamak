import copy
import numpy as np

import scipy.interpolate as interp

if __name__ == "__main__":
    import matplotlib
    #matplotlib.use('TkAgg')
    matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as plt

import probe
reload(probe)

import config_xpr as config

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
        self.mdsfmt = '_s = augsignal(%d,"%s","%%s","AUGD",*,*,*,*,*,"raw")' % (self.shn, diag)

        self.datadeco = '%s; word_unsigned(data(_s))'
        self.timedeco = '%s; dim_of(_s)'
        self.sizedeco = '%s; size(_s)'


class IOFileAUG(IOFile):
    def __init__(self, shn, diag='XPR'):
        IOFile.__init__(self, shn=shn, suffix="_"+diag, subdir="AUG")


class DigitizerXPR(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock, name='XPR')

        self.IO_mds = IOMdsAUG(shn, sock, diag='XPR')
        self.IO_file = IOFileAUG(shn, diag='XPR')
        self.nodes = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 't')

        #f = 20./16283
        #offs = (8759, 8830, 8756, 8770, 8759, 8747, 8746, 8752, 0)
        #self.amp = {node: Amp(fact=f, offs=-f*o) for node, o in zip(self.nodes, offs)}
        self.amp = {node: amp14Bit.copy() for node in self.nodes}
        self.amp['t'] = ampUnity.copy()

    def update_window(self):
        """Detect failing data readout and cut window
        """
        t = self.x['t']
        dt = np.abs(np.diff(t))
        cnd = dt > 2*dt[0]
        iM = np.argmax(cnd)

        stop = None
        if cnd[iM] == True:
            stop = iM+1
        self.window = slice(None, stop)
        
    def calib(self):
        self.update_window()
        Digitizer.calib(self)


class DigitizerLPS(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock, name='LPS')

        self.IO_mds = IOMdsAUG(shn, sock, diag='LPS')
        self.IO_file = IOFileAUG(shn, diag='LPS')
        self.nodes = ('CUR1', 'VOL1', 'CUR2', 'VOL2', 'VOL3', 'VOL4', 't')

        self.window = slice(2048, None)

        #f = 10./4095
        #offs = (0,0,0,0,0,0,0)
        #self.amp = {node: Amp(fact=f, offs=-f*o) for node, o in zip(self.nodes, offs)}
        self.amp = {node: amp12Bit.copy() for node in self.nodes}
        self.amp['t'] = ampUnity.copy()

        for node in ('CUR1', 'VOL3'):
            self.amp[node] *= ampInv


class ProbeXPR(Probe):
    def __init__(self, shn, sock=None, dig=None):
        
        reload(config)
        
        self.config = config.campaign.find_shot(shn)
        if dig is None:
            dig = self.config.dig
        
        if dig == 'LPS':
            DigitizerClass = DigitizerLPS
        else:
            DigitizerClass = DigitizerXPR

        digitizer = DigitizerClass(shn, sock)
        Probe.__init__(self, digitizer)

    def mapsig(self):
        
        x = self.x[self.digitizer.window]

        self.S = self.config.mapsig(x, self.digitizer.name)

    def calib(self):        
        self.config.calib(self.S, self.digitizer.name)

        s = slice(5000)
        for S in self.get_type('Current'):
            S.norm_to_region(s)

        for S in self.get_type('Voltage'):
            S.norm_to_region(s)

        self.S['V'] = self.S['I1'].V
        self.S['It'] = self.S['I1'] + self.S['I2']

    def get_meas(self, Isat, Vf, Te, meas):
        head = self.config.head

        #i_up = head.get_tip_number_by_position('lower left')
        #i_dn = head.get_tip_number_by_position('lower right')

        tips = self.config.head.tips

        II = self.get_type('Current')
        for i in xrange(len(II)):
            num = II[i].number
            if num > 0:
                tip = tips[num - 1]
                j = Isat[i] / tip.area
                if tip.pos == 'lower left':
                    meas.jp = j
                elif tip.pos == 'lower right':
                    meas.jm = j
            else:
                meas.Vf = Vf[i]
                meas.Te = Te[i]

    def position_calib(self):
        R = self['R']
        sp = interp.UnivariateSpline(R.t, R.x, s=10)
        Rs = sp(R.t)
        return Rs

    def voltage_calib(self):
        self.load(trim=False, calib=False)
        return self['V'].x.mean()
   

def get_dwell_params():
    shots = config.campaign.shots_with_min_stars('*')
    for shot in shots:
        XPR = ProbeXPR(shn=shot)
        tM, RM = XPR.get_dwell_params()
        print "%d: tM = %s s, RM = %s m" % (shot, tM, RM)


if __name__ == "__main__":
    shn = 28469

    XPR = ProbeXPR(shn=shn)
    XPR.analyze(plunge=0)

    XPR.plot()
    plt.show()



