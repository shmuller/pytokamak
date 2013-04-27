from digitizer import TdiError, IOMds, IOFile, Digitizer

from LP.sig import Amp

class IOMdsD3D(IOMds):
    def __init__(self, *args, **kw):
        diag = kw.pop('diag', '')
        suffix = kw.pop('suffix', '')

        IOMds.__init__(self, *args, **kw)
        self.mdsport = "8020"
        self.mdstree = diag
        self.mdsfmt = '\%s' + suffix


class IOFileD3D(IOFile):
    def __init__(self, *args, **kw):
        kw.setdefault('subdir', "D3D")
        IOFile.__init__(self, *args, **kw)


class DigitizerD3DEFIT(Digitizer):
    def __init__(self, shn):
        Digitizer.__init__(self, shn, name='EFIT')
        self.IO_mds = IOMdsD3D(shn, diag='EFIT01')
        self.IO_file = IOFileD3D(shn, suffix='_D3D', group='EFIT01')
        self.tnode, self.tunits = 'GTIME', 's'
        self.nodes = ('GTIME', 'R', 'Z', 'PSIRZ', 
                      'RMAXIS', 'ZMAXIS', 'SSIMAG', 'SSIBRY', 'BDRY')

        self.amp = dict(GTIME=Amp(1e-3))

    def get_R_z_psi(self):
        return self.x['R'], self.x['Z'], self['PSIRZ']

    def get_R_z_psi_special(self, spec):
        if spec == 'magnaxis':
            return self['RMAXIS'], self['ZMAXIS'], self['SSIMAG']
        elif spec == 'xpoint':
            return None, None, self['SSIBRY']


