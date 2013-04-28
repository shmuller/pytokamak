import numpy as np
import os

from LP.sig import memoized_property

from digitizer import TdiError, IOMds, IOFile, Digitizer

from sm_pyplot.tight_figure import get_axes

from matplotlib.path import Path
from matplotlib.collections import PathCollection

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


class DigitizerAUGYGC(DigitizerAUG):
    def __init__(self, shn):
        DigitizerAUG.__init__(self, shn, diag='YGC', 
                nodes=('RGC2', 'zGC2', 'inxbeg', 'inxlen', 'inxlps', 'ixplin',
                       'RrGC', 'zzGC'))
    
    def calib(self):
        DigitizerAUG.calib(self)

        self.x['chGCnm'] = np.array(('PCup', 'PClow', 'PCled3', 'PCled4', 'TPRT', 'TPRB', 
                'TPLT', 'TPLB', 'LIM09', 'SBi', 'ICRHa', 'VESiR', 'LIaux13', 'LIaux14', 
                'LIaux15', 'D2c.i1', 'D2c.i2', 'D2c.TPib', 'D2c.TPi', 'D2c.TPic', 
                'D2c.domL', 'D2c.dome', 'D2c.domR', 'D2d.BG10', 'D2d.BG1', 'D2d.BG2', 
                'D2d.Bl3', 'D2d.Bl2', 'D2d.Bl1', 'D2.SBiu', 'TPLT_1', 'TPLT_2', 'TPLT_3',
                'TPLT_4', 'TPLT_5', 'TPRT_1', 'TPRT_2', 'TPRT_3', 'TPRT_4', 'TPRT_5',
                'D2d.Bu4', 'D2d.Bu3', 'D2d.Bu2', 'D2d.Bu1'))

    @memoized_property
    def label(self):
        x = self.x
        return x['chGCnm'][x['ixplin'] > 0]

    @memoized_property
    def xy(self):
        x = self.x
        i = x['inxbeg'] - 1
        ij = np.c_[i[:-1], i[1:]][x['ixplin'] > 0]
        xy = np.c_[x['RrGC'], x['zzGC']]
        return [xy[i:j] for i, j in ij]

    def plot(self, ax=None, unfilled=(4, 5), col='k'):
        ax = get_axes(ax, xlab="R (m)", ylab="z (m)")
        ax.set_aspect('equal')
        ax.set_xlim((0.5, 2.5))
        ax.set_ylim((-1.5, 1.5))

        paths = [Path(xy) for i, xy in enumerate(self.xy) if i not in unfilled]
        pc = PathCollection(paths, facecolors=col, edgecolors='none')
        ax.add_collection(pc)

        paths = [Path(self.xy[i]) for i in unfilled]
        pc = PathCollection(paths, facecolors='none', edgecolors=col)
        ax.add_collection(pc)
        return ax


dig_YGC = DigitizerAUGYGC(shn=25891)



