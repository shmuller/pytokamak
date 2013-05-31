from sm_pyplot import vtk_plot

_VtkWindow = vtk_plot.VtkWindow
_VtkRotatingPolygon = vtk_plot.VtkRotatingPolygon
_VtkContour = vtk_plot.VtkContour

class VtkWindow(_VtkWindow):
    def __init__(self, **kw):
        _VtkWindow.__init__(self, campos=(0., -10., 5.), parscale=2., **kw)


class VtkRotatingPolygon(_VtkRotatingPolygon):
    def __init__(self, *args, **kw):
        _VtkRotatingPolygon.__init__(self, *args, VtkWindowClass=VtkWindow, **kw)


class VtkContour(_VtkContour):
    def __init__(self, *args, **kw):
        _VtkContour.__init__(self, *args, VtkWindowClass=VtkWindow, **kw)


