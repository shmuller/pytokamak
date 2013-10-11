import numpy as np
import os

from utils.utils import memoized_property
from utils.sig import Amp, get_axes

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
        self.timedeco = 'f_float(dim_of(%s))'

        self.sfhfmt = 'augparam(%d,"%s","%%s","%%s")' % (self.shn, diag)

    def get_param(self, *args):
        return self.mdsvalue(self.sfhfmt % args)


class IOFileAUG(IOFile):
    def __init__(self, *args, **kw):
        kw.setdefault('subdir', "AUG")
        IOFile.__init__(self, *args, **kw)


class DigitizerAUG(Digitizer):
    def __init__(self, shn, diag, suffix='_AUG', group=None, raw=False, 
                 IOMdsClass=IOMdsAUG, IOFileClass=IOFileAUG, **kw):
        if group is None:
            group = diag
        
        Digitizer.__init__(self, shn, name=diag, 
                IO_mds = IOMdsClass(shn, diag=diag, raw=raw),
                IO_file = IOFileClass(shn, suffix=suffix, group=group), **kw)

    def get_node(self, node, **kw):
        try:
            x = Digitizer.get_node(self, node, **kw)
            if x.ndim > 1:
                perm = np.roll(np.arange(x.ndim), 1)
                x = np.ascontiguousarray(x.transpose(perm))
            return x
        except TdiError:
            return np.zeros(1)


class DigitizerAUGMAC(DigitizerAUG):
    def __init__(self, shn):
        DigitizerAUG.__init__(self, shn, diag='MAC', nodes=('Ipolsola', 'Ipolsoli'))
        self.dig_Tdiv = DigitizerAUG(shn, diag='MAC', group='MAC/Tdiv', nodes=('Tdiv',))

    def __getitem__(self, indx):
        try:
            return DigitizerAUG.__getitem__(self, indx)
        except KeyError:
            return self.dig_Tdiv[indx]


class DigitizerAUGDCR(DigitizerAUG):
    alias = ('H-1', 'H-2', 'H-3', 'H-4', 'H-5')

    def __init__(self, shn):
        DigitizerAUG.__init__(self, shn, diag='DCR', nodes=('dbl_res',))

    def __getitem__(self, indx):
        try:
            return DigitizerAUG.__getitem__(self, indx)
        except KeyError:
            try:
                i = self.alias.index(indx)
            except:
                raise KeyError(indx)
            S = DigitizerAUG.__getitem__(self, 'dbl_res')[:, 5 + i]
            S.update(name=indx, label=indx)
            return S


class IOMdsAUGMIR(IOMdsAUG):
    numA = (1, 2, 3, 4, 5, 6, 26)
    numD = (7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 20, 21, 22)
    numE = (23, 24, 25, 27, 28, 29, 30, 31, 32)
    
    formatter = "C09-{:02d}".format
    nodesA = map(formatter, numA)
    nodesD = map(formatter, numD)
    nodesE = map(formatter, numE)

    node_groups = (nodesA, nodesD, nodesE)
    diags = ('MHA', 'MHD', 'MHE')

    node2diag = dict()
    for nodes, diag in zip(node_groups, diags):
        for node in nodes:
            node2diag[node] = diag

    def __init__(self, shn, **kw):
        kw['diag'] = '{diag}'
        IOMdsAUG.__init__(self, shn, **kw)

    def get_mdsbasestr(self, node):
        return self.mdsfmt.format(diag=self.node2diag[node]) % node
        
    def get_param(self, node, param):
        sfhfmt = self.sfhfmt.format(diag=self.node2diag[node])
        return self.mdsvalue(sfhfmt % ('C' + node, param))


class DigitizerAUGMIR(DigitizerAUG):
    def __init__(self, shn, **kw):
        kw.setdefault('t0', 0.)
        kw.setdefault('t1', 6.)
        kw.setdefault('s', slice(None, None, 4))

        nodes = sorted(IOMdsAUGMIR.nodesA + IOMdsAUGMIR.nodesD + IOMdsAUGMIR.nodesE)
        DigitizerAUG.__init__(self, shn, diag='MIR', suffix='_MIR', nodes=nodes, 
                              IOMdsClass=IOMdsAUGMIR, **kw)


amp_2pi     = Amp(fact=0.5/np.pi, offs=0)
amp_mu0_2pi = Amp(fact=2e-7, offs=0)

class DigitizerAUGFPP(DigitizerAUG):
    def __init__(self, shn, diag='FPP'):
        DigitizerAUG.__init__(self, shn, diag=diag, 
                nodes=('PFM', 'Ri', 'Zj', 'ikCAT', 'RPFx', 'zPFx', 'PFxx',
                       'Lpf', 'PFL', 'TFLx', 'Qpsi', 'Jpol', 'Pres', 'Vol', 'Area', 'CLE'))
        
        self.amp.update(PFM=amp_2pi, PFxx=amp_2pi, PFL=amp_2pi, Jpol=amp_mu0_2pi)
        self.alias = dict(psii='PFL', q='Qpsi')
        self.alias_primed = dict(f='Jpol', p='Pres', V='Vol', A='Area')
        
        self.mapspec = dict(magnaxis=0, xpoint=1, innerlim=2, xpoint2=3, outerlim=4)
    
    def __getitem__(self, indx):
        get = DigitizerAUG.__getitem__
        try:
            return get(self, self.alias[indx])
        except KeyError:
            pass

        indx, prime, tail = indx.partition('prime')
        try:
            return get(self, self.alias_primed[indx])[:,bool(prime)::2]
        except KeyError:
            return get(self, indx)

    def get_R_z_psi(self):
        return self.x['Ri'][0], self.x['Zj'][0], self['PFM']

    def get_R_z_psi_special(self, spec):
        i = self.mapspec[spec]
        return self['RPFx'][:,i], self['zPFx'][:,i], self['PFxx'][:,i]


class DigitizerAUGEQI(DigitizerAUGFPP):
    def __init__(self, shn, diag='EQI'):
        DigitizerAUGFPP.__init__(self, shn, diag=diag)

        self.nodes += ('FFP', 'CLD', 'Rinv', 'R2inv', 'Bave', 'B2ave', 'FTRA')
        self.amp.update(FFP=amp_mu0_2pi*amp_mu0_2pi)
        self.alias.update(ffprime='FFP', Bt='Bave')


class DigitizerAUGEQH(DigitizerAUGEQI):
    def __init__(self, shn, diag='EQH'):
        DigitizerAUGEQI.__init__(self, shn, diag=diag)


eqi_digitizers = dict(FPP=DigitizerAUGFPP, EQI=DigitizerAUGEQI, EQH=DigitizerAUGEQH)


class DigitizerAUGYGC(DigitizerAUG):
    def __init__(self, shn):
        DigitizerAUG.__init__(self, shn, diag='YGC', 
                nodes=('RGC2', 'zGC2', 'inxbeg', 'inxlen', 'inxlps', 'ixplin',
                       'RrGC', 'zzGC'))

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
    def _xy(self):
        xy = np.c_[self.x['RrGC'], self.x['zzGC']]

        # fix inner limiter
        perm = np.concatenate((np.arange(21, 42), np.arange(0, 22)))
        xy[74:117] = xy[74 + perm]
        return xy

    @memoized_property
    def _lbdry(self):
        lbdry = self.x['inxlps'].copy()
        
        # fix inner limiter and upper tiles
        lbdry[9] = 21
        lbdry[30:34] = 2
        lbdry[34] = 7
        lbdry[37:40] = 2
        return lbdry
        
    def _get_xy(self, lenname, mask=True):
        x = self.x
        i = x['inxbeg'][:-1] - 1
        j = i + x[lenname]
        ij = np.c_[i, j]
        if mask:
            ij = ij[x['ixplin'] > 0]
        return [self._xy[i:j] for i, j in ij]

    @memoized_property
    def xy(self):
        return self._get_xy(lenname='inxlen')
    
    @memoized_property
    def xy_bdry(self):
        i = self.x['inxbeg'][:-1] - 1
        ij = np.c_[i, i + self._lbdry]
        idx = np.r_[15:28, 10, 41:44, 36:40, 30:35, 9]
        return [self._xy[i:j] for i, j in ij[idx]]

    def plot(self, ax=None, unfilled=(4, 5), edgecolors='k', facecolors=(0.75, 0.75, 0.75)):
        from matplotlib.path import Path
        from matplotlib.collections import PathCollection

        ax = get_axes(ax, xlab="R (m)", ylab="z (m)")
        ax.set_aspect('equal')
        ax.set_xlim((0.5, 2.5))
        ax.set_ylim((-1.5, 1.5))

        paths = [Path(xy) for i, xy in enumerate(self.xy) if i not in unfilled]
        pc = PathCollection(paths, facecolors=facecolors, edgecolors=edgecolors)
        ax.add_collection(pc)

        paths = [Path(self.xy[i]) for i in unfilled]
        pc = PathCollection(paths, facecolors='none', edgecolors=edgecolors)
        ax.add_collection(pc)
        return ax


dig_YGC = DigitizerAUGYGC(shn=25891)



