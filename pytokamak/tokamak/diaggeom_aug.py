import numpy as np
import os

from pytokamak.utils.utils import memoized_property, BoundingBox
from pytokamak.utils.sig import get_axes

from vtk_aug import VtkProxy, VtkWindow, VtkRotatingPolygon

class Vessel(VtkProxy):
    def __init__(self):
        from digitizer_aug import dig_YGC
        self.digitizer = dig_YGC

    @memoized_property
    def bdry(self):
        return np.concatenate(self.digitizer.xy_bdry)

    @memoized_property
    def vtk(self):
        kw = dict(phi1=3./2.*np.pi)
        return [VtkRotatingPolygon(*xy.T, **kw) for xy in self.digitizer.xy]

    def get_bbox(self):
        return BoundingBox(self.bdry.min(axis=0), self.bdry.max(axis=0))

    def plot(self, *args, **kw):
        return self.digitizer.plot(*args, **kw)

    def prerender(self, **kw):
        alpha = np.ones(len(self.vtk))
        alpha[[0, 1]] = 0.2
        alpha[[2, 4, 5]] = 0

        win = VtkWindow(axes='2d', **kw)
        for rpoly, a in zip(self.vtk, alpha):
            if a > 0:
                win = rpoly.prerender(win=win, alpha=a)
        return win


class DCNGeom:
    # from DCN_Strahlengaenge20071212.pdf
    DCN_all = np.array(
        [[(1.0060, 0.1447), (2.1662, 0.1351), (2.1667, 0.1715), (2.1664, 0.1533)],
         [(1.0067, 0.3154), (2.1697, 0.3054), (2.1698, 0.3418), (2.1695, 0.3236)],
         [(1.0065,-0.1465), (2.1611,-0.1778), (2.1617,-0.1429), (2.1614,-0.1603)],
         [(1.1287, 1.0566), (2.1671, 0.1754), (2.1662, 0.1299), (2.1668, 0.1528)],
         [(1.0946, 0.8032), (2.1712, 0.4246), (2.1718, 0.4606), (2.1715, 0.4426)]])

    names = np.array(('H-1', 'H-2', 'H-3', 'H-4', 'H-5'))
    DCN_rays = DCN_all[:,[0, 3]]

    def get_rays(self, names=names):
        idx = np.searchsorted(self.names, names)
        return self.DCN_rays[idx]

    def plot_rays(self, ax=None, names=names, **kw):
        kw.setdefault('linewidth', 1)
        rays = self.get_rays(names)
        ax = get_axes(ax)
        return ax.plot(rays[:,:,0].T, rays[:,:,1].T, **kw)

    def plot(self, ax=None, names=names, **kw):
        rays = self.get_rays(names)
        lines = self.plot_rays(ax, names, **kw)

        offs = np.array([-0.05, 0])
        for lab, xy, l in zip(names, rays[:,0,:], lines):
            ax.annotate(lab, xycoords='data', xy=xy + offs , fontsize=15,
                ha='right', va='center', color=l.get_color(), backgroundcolor='w')
        return ax


class DopplerGeom:
    antenna = np.array([(2.2717,-0.5087), (2.1998,-0.4665), 
                        (2.2528,-0.4033), (2.3066,-0.4667)])

    arrow = np.array([(2.2491,-0.4582), (1.8415,-0.0897)])

    facecolor = np.array((242, 110, 69)) / 255.
    edgecolor = np.array((235,  34, 44)) / 255.

    def plot(self, ax=None, facecolor=facecolor, edgecolor=edgecolor):
        from matplotlib.patches import Polygon
        ax = get_axes(ax)

        arrowprops = dict(edgecolor='none', facecolor=facecolor, shrink=0)
        ax.annotate('', xytext=self.arrow[0], xy=self.arrow[1], arrowprops=arrowprops)

        ax.add_patch(Polygon(self.antenna, facecolor=facecolor, edgecolor=edgecolor))

        ax.annotate("V-band\nX-mode\nfixed ant.", xycoords='data', xy=(2.25, -0.35), 
            fontsize=15, ha='left', va='bottom', color=facecolor, backgroundcolor='w')
        return ax

    @memoized_property
    def _transform_from_image(self):
        from matplotlib.transforms import Bbox, BboxTransform

        ref_im = np.array([(15.152, 57.079), (15.152, 65.091),
                           (12.949, 65.091), (12.949, 62.575),
                           ( 5.613, 62.575), ( 5.613, 60.587),
                           (12.949, 60.587), (12.949, 57.079)])

        ref = Vessel().digitizer.xy[2]

        bbox_im, bbox = Bbox.unit(), Bbox.unit()
        bbox_im.update_from_data_xy(ref_im)
        bbox.update_from_data_xy(ref)
        trans = BboxTransform(bbox_im, bbox)
        return trans

    def _from_image(self):
        antenna_im = np.array([(133.848, 82.243), (126.594, 86.477), 
                               (131.941, 92.802), (137.367, 86.454)])

        arrow_im = np.array([(131.566,87.309), (90.461,124.212)])

        trans = self._transform_from_image
        antenna = trans.transform(antenna_im)
        arrow = trans.transform(arrow_im)
        return dict(antenna=antenna, arrow=arrow)


class MIRGeom:
    def __init__(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data_aug')
        fname = os.path.join(data_dir, 'angle.dat')

        names = ('phi', 'theta', 'R', 'z', 'shift', 'dtheta', 'Ohm', 'Aeff')
        dtype = np.dtype([('name', np.str_, 8)] + zip(names, (np.float32,)*len(names)))
        self.record = np.loadtxt(fname, dtype, skiprows=2)

        self.lut = dict(zip(self.record['name'], self.record))

    def get_C09(self):
        import re
        matcher = re.compile('C09-\d\d').match
        cnd = np.array(map(matcher, self.record['name']), dtype=bool)
        rec = self.record[cnd]
        lut = dict(zip(rec['name'], rec))
        return lut

    def plot(self, ax=None):
        lut = self.get_C09()

        ax = get_axes(ax)
        for coil in lut.itervalues():
            ax.plot(coil['R'], coil['z'], 'rs')
        return ax
            


