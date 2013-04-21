import numpy as np

from matplotlib.pyplot import contour

from mdsclient import mdsopen, mdsvalue, mdsclose

from sm_pyplot.vtk_contour import VtkContour


class FluxSurf:
    def __init__(self, vtkCtr):
        self.vtkCtr = vtkCtr

    def plot(self):
        self.vtkCtr.plot()


class EFIT:
    def __init__(self, sock, shn, ti=None, dt=np.zeros(2)):
        self.sock, self.shn, self.ti, self.dt = sock, shn, ti, dt

    def load_node(self, node, idx=''):
        nodestr = "_x=data(%s); _x%s" % (node, idx)
        return mdsvalue(self.sock, nodestr)

    def load(self):
        mdsopen(self.sock,'EFIT01',self.shn)

        # constant quantities
        self.psin = self.load_node('\PSIN')
        self.r    = self.load_node('\R')
        self.z    = self.load_node('\Z')
        self.t    = self.load_node('\GTIME')

        if self.ti is not None:
            self.ind = np.abs(self.t-self.ti).argmin()
    
            end = str(self.ind) + ']'
            idx = ['['+end, '[*,'+end, '[*,*,'+end]
        else:
            idx = ['','','']

        # 1D quantities
        self.bcentr = self.load_node('\EFIT_G_EQDSK.BCENTR',idx[0])
        self.rmaxis = self.load_node('\RMAXIS',idx[0])
        self.zmaxis = self.load_node('\ZMAXIS',idx[0])
        self.ssimag = self.load_node('\SSIMAG',idx[0])
        self.ssibry = self.load_node('\SSIBRY',idx[0])
        self.cpasma = self.load_node('\CPASMA',idx[0])

        # 2D quantities
        self.ffprim = self.load_node('\FFPRIM',idx[1])
        self.pprime = self.load_node('\PPRIME',idx[1])
        self.pres   = self.load_node('\PRES'  ,idx[1])
        self.fpol   = self.load_node('\FPOL'  ,idx[1])
        self.qpsi   = self.load_node('\QPSI'  ,idx[1])
        self.epoten = self.load_node('\EPOTEN',idx[1])
        self.rhovn  = self.load_node('\RHOVN' ,idx[1])

        # 3D quantities
        self.bdry  = self.load_node('\BDRY' ,idx[2]).T
        self.psirz = self.load_node('\PSIRZ',idx[2])

        self.psi_n = (self.ssimag-self.psirz)/(self.ssimag-self.ssibry)

        mdsclose(self.sock)

    def getFluxSurf(self, Lvls=None):
        x = np.double(self.r)
        y = np.double(self.z)
        z = np.zeros(1,'d')
        f = np.double(self.psi_n)
        vtkCtr = VtkContour(x,y,z,f,Lvls)
        vtkCtr.contour()

        return FluxSurf(vtkCtr)

        #cs = contour(self.r,self.z,self.psi_n,[lvl,lvl])

        #p = cs.collections[0].get_paths()[0]
        #return FluxSurf(p.vertices)




