import copy
import os
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
Signal = probe.Signal
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
        raw  = kw.pop('raw', False)

        IOMds.__init__(self, *args, **kw)
        
        if os.uname()[1] == 'plaspc04':
            self.mdsserver, self.mdsport = "localhost", "8001"
        else:
            self.mdsserver, self.mdsport = "mdsplus.aug.ipp.mpg.de", "8000"

        if raw:
            mdsfmt = '_s = augsignal(%d,"%s","%%s","AUGD",*,*,*,*,*,"raw")'
            self.datadeco = '%s; word(data(_s))'
        else:
            mdsfmt = '_s = augsignal(%d,"%s","%%s","AUGD")'
            self.datadeco = '%s; data(_s)'

        self.mdsfmt =  mdsfmt % (self.shn, diag)

        self.timedeco = '%s; dim_of(_s)'
        self.sizedeco = '%s; size(_s)'


class IOFileAUG(IOFile):
    def __init__(self, shn, diag='XPR'):
        IOFile.__init__(self, shn=shn, suffix="_"+diag, subdir="AUG")


class DigitizerXPR(Digitizer):
    def __init__(self, shn, sock=None, raw=False):
        Digitizer.__init__(self, shn, sock, name='XPR')

        self.IO_mds = IOMdsAUG(shn, sock, diag='XPR', raw=raw)
        self.IO_file = IOFileAUG(shn, diag='XPR')
        self.nodes = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 't')

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


class DigitizerXPRRaw(DigitizerXPR):
    def __init__(self, shn, sock=None):
        DigitizerXPR.__init__(self, shn, sock, raw=True)

        self.amp = {node: amp14Bit.copy() for node in self.nodes[:-1]}


class DigitizerXPRRawPos(DigitizerXPRRaw):
    def __init__(self, shn, sock=None):
        DigitizerXPRRaw.__init__(self, shn, sock)
        self.nodes += ('PosL', 'PosH')
        self.more_nodes = ('Pos',)

    def _load_raw_factory(name):
        def load_raw(self):
            getattr(DigitizerXPRRaw, name)(self)

            PosL = self.x['PosL'].astype(np.int16).view(np.uint16)
            PosH = self.x['PosH'].astype(np.int16).view(np.uint16) & 0b0011111111111111
        
            Pos  = np.c_[PosL, np.r_[PosH[0], PosH[:-1]]]
        
            self.x['Pos'] = Pos.copy().view(np.uint32)[:,0]
            return self.x
        return load_raw

    load_raw_mds  = _load_raw_factory('load_raw_mds')
    load_raw_file = _load_raw_factory('load_raw_file')


class DigitizerLPS(Digitizer):
    def __init__(self, shn, sock=None, raw=False):
        Digitizer.__init__(self, shn, sock, name='LPS')

        self.IO_mds = IOMdsAUG(shn, sock, diag='LPS', raw=raw)
        self.IO_file = IOFileAUG(shn, diag='LPS')
        self.nodes = ('CUR1', 'VOL1', 'CUR2', 'VOL2', 'VOL3', 'VOL4', 't')

        self.window = slice(2048, None)


class DigitizerLPSRaw(DigitizerLPS):
    def __init__(self, shn, sock=None):
        DigitizerLPS.__init__(self, shn, sock, raw=True)

        self.amp = {node: amp12Bit.copy() for node in self.nodes[:-1]}

        for node in ('CUR1', 'VOL3'):
            self.amp[node] *= ampInv


class DigitizerXPOS(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock, name='XPOS')

        self.IO_mds = IOMdsAUG(shn, sock, diag='LPS')
        self.IO_file = IOFileAUG(shn, diag='LPS_XPOS')
        self.nodes = ('XPOS', 't')


class DigitizerLPSOld(DigitizerLPS):
    def __init__(self, shn, sock=None):
        DigitizerLPS.__init__(self, shn, sock)

        self.dig_xpos = DigitizerXPOS(self.shn, self.sock)

        self.more_nodes = ('XPOS',)

        for node in ('CUR1', 'CUR2', 'VOL3'):
            self.amp[node] = ampInv.copy()

    def _load_raw_factory(name):
        def load_raw(self):
            getattr(DigitizerLPS, name)(self)
            getattr(self.dig_xpos, name)()

            x = self.dig_xpos.x
            R = Signal(x['XPOS'].astype('d'), x['t'].astype('d'))
            R.despike()

            self.x['XPOS'] = R(self.x['t'])
            return self.x
        return load_raw

    load_raw_mds  = _load_raw_factory('load_raw_mds')
    load_raw_file = _load_raw_factory('load_raw_file')

    def save(self):
        DigitizerLPS.save(self)
        self.dig_xpos.save()


DigitizerClasses = dict(
        LPS = DigitizerLPSRaw,
        XPR = DigitizerXPRRaw,
        XPR_pos = DigitizerXPRRawPos,
        LPS_old = DigitizerLPSOld)


class ProbeXPR(Probe):
    def __init__(self, shn, sock=None, dig=None):
        
        reload(config)
        
        self.config = config.campaign.find_shot(shn)
        if dig is None:
            dig = self.config.dig

        digitizer = DigitizerClasses[dig](shn, sock)
        Probe.__init__(self, digitizer)

    def mapsig(self):
        
        x = self.x[self.digitizer.window]

        self.S = self.config.mapsig(x, self.digitizer.name)

    def calib(self):
        S = self.S
        self.config.calib(S, self.digitizer.name)

        s = slice(5000)
        for I in self.get_type('Current'):
            I.norm_to_region(s)

        for V in self.get_type('Voltage'):
            V.norm_to_region(s)

        S['Rs'] = S['R'].copy().mediansmooth(100)
        S['V']  = S['I1'].V
        S['It'] = S['I1'] + S['I2']

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



