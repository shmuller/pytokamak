import numpy as np
import os

from digitizer import TdiError, IOMds, IOFile, Digitizer

class IOMdsAUG(IOMds):
    def __init__(self, shn, diag='XPR', raw=False):
        # augdiag(_shot, _diag, _signame, _experiment, _edition, 
        #   _t1, _t2, _oshot, _oedition, _qual)

        IOMds.__init__(self, shn)
        
        if os.uname()[1] == 'plaspc04':
            self.mdsserver, self.mdsport = "localhost", "8001"
        else:
            self.mdsserver, self.mdsport = "mdsplus.aug.ipp.mpg.de", "8000"

        self.mdsplaceholder = 'f_float($)'

        fmtargs = '%d,"%s","%%s","AUGD",*,f_float($),f_float($)' % (self.shn, diag)
        if raw:
            self.mdsfmt = 'augdiag(%s,*,*,"raw")' % fmtargs
            self.datadeco = 'word(data(%s))'
        else:
            self.mdsfmt = 'augdiag(%s)' % fmtargs


class IOFileAUG(IOFile):
    def __init__(self, *args, **kw):
        kw.setdefault('subdir', "AUG")
        IOFile.__init__(self, *args, **kw)


class DigitizerAUG(Digitizer):
    def __init__(self, shn, diag, nodes, suffix='_AUG', group=None, raw=False, **kw):
        Digitizer.__init__(self, shn, name=diag, **kw)
        if group is None:
            group = diag

        self.IO_mds = IOMdsAUG(shn, diag=diag, raw=raw)
        self.IO_file = IOFileAUG(shn, suffix=suffix, group=group)

        if self.tnode not in nodes:
            nodes += (self.tnode,)
        self.nodes = nodes

    def load(self, **kw):
        try:
            return Digitizer.load(self, **kw)
        except TdiError:
            return {node: np.zeros(1) for node in self.nodes}

    def calib(self):
        for node, x in self.x.iteritems():
            if x.ndim > 1:
                self.x[node] = x.transpose(np.roll(np.arange(x.ndim), 1))
        Digitizer.calib(self)


class DigitizerAUGMAC(DigitizerAUG):
    def __init__(self, shn):
        DigitizerAUG.__init__(self, shn, diag='MAC', nodes=('Ipolsola', 'Ipolsoli'))
        self.dig_Tdiv = DigitizerAUG(shn, diag='MAC', group='MAC/Tdiv', nodes=('Tdiv',))

    def __getitem__(self, indx):
        try:
            return DigitizerAUG.__getitem__(self, indx)
        except KeyError:
            return self.dig_Tdiv[indx]


class DigitizerAUGEQI(DigitizerAUG):
    def __init__(self, shn):
        DigitizerAUG.__init__(self, shn, diag='EQI', 
                              nodes=('Ri', 'Zj', 'PFM', 'ikCAT', 'RPFx', 'zPFx', 'PFxx'))
        self.mapspec = dict(magnaxis=0, xpoint=1, innerlim=2, xpoint2=3, outerlim=4)
   
    def get_R_z_psi(self):
        return self.x['Ri'][0], self.x['Zj'][0], self['PFM']

    def get_R_z_psi_special(self, spec):
        i = self.mapspec[spec]
        return self['RPFx'][:,i], self['zPFx'][:,i], self['PFxx'][:,i]


