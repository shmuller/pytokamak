from matplotlib import get_backend

if get_backend().lower() == 'qt4agg':
    from sm_pyplot import vtk_plot

    _VtkCamera = vtk_plot.VtkCamera
    _VtkWindow = vtk_plot.VtkWindow
    _VtkRotatingPolygon = vtk_plot.VtkRotatingPolygon
    _VtkContour = vtk_plot.VtkContour
    _VtkPolyline = vtk_plot.VtkPolyline
else:
    class Dummy:
        def __init__(self, *args, **kw):
            pass

    _VtkCamera = _VtkWindow = _VtkRotatingPolygon = _VtkContour = _VtkPolyline = Dummy


class VtkProxy:
    def prerender(self, *args, **kw):
        return self.vtk.prerender(*args, **kw)

    def render(self, *args, **kw):
        return self.prerender(*args, **kw).render()


class VtkCamera(_VtkCamera):
    def __init__(self, **kw):
        _VtkCamera.__init__(self, campos=(0., -10., 5.), parscale=3., **kw)


class VtkWindow(_VtkWindow):
    def __init__(self, **kw):
        _VtkWindow.__init__(self, bgcolor=(1., 1., 1.), VtkCameraClass=VtkCamera, **kw)


class VtkRotatingPolygon(_VtkRotatingPolygon):
    def __init__(self, *args, **kw):
        _VtkRotatingPolygon.__init__(self, *args, VtkWindowClass=VtkWindow, **kw)


class VtkContour(_VtkContour):
    def __init__(self, *args, **kw):
        _VtkContour.__init__(self, *args, VtkWindowClass=VtkWindow, 
                             mode='contour', **kw)


class VtkPolyline(_VtkPolyline):
    def __init__(self, *args, **kw):
        _VtkPolyline.__init__(self, *args, VtkWindowClass=VtkWindow, **kw)


