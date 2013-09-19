import numpy as np

from utils.utils import memoized_property, BoundingBox

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



