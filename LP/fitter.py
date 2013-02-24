import numpy as np

import scipy.optimize as opt
import scipy.odr as odr

import minpack as mp

from sig import memoized_property, get_axes

class FitterError(Exception):
    pass

class Fitter:
    def __init__(self, x, y, args=(), engine='custom', 
            use_rms=True, use_diff=True, use_fast=True):
        self.x, self.y, self.args = x, y, args

        self.do_var = None

        self.buf = np.empty_like(x)

        if use_fast and hasattr(self, 'fitfun_fast'):
            self.fun = self.factory_buf(self.fitfun_fast)
        else:
            self.fun = self.fitfun

        if use_diff and hasattr(self, 'fitfun_diff'):
            self.fun_diff = self.factory_buf(self.fitfun_diff)
        else:
            self.fun_diff = self.factory_diff(self.fun)

        if use_rms and hasattr(self, 'fitfun_rms'):
            self.fun_rms = self.fitfun_rms
        else:
            self.fun_rms = self.factory_rms(self.fun_diff)

        self.set_engine(engine)

    # overload
    def set_OK(self):
        self.OK = True

    def set_norm(self):
        self.X, self.Y = self.x, self.y

    def set_unnorm(self):
        self.p, self.p0 = self.P, self.P0

    def set_guess(self):
        self.P0 = 0.

    @classmethod
    def fitfun(cls, P, X):
        pass

    # fitfunction factories
    def factory_buf(self, fun):
        def fun_buf(p, x, *args):
            return fun(p, x, self.buf[:x.size], *args)
        return fun_buf

    def factory_diff(self, fun_buf):
        def fun_diff(p, x, y, *args):
            return fun_buf(p, x, *args) - y
        return fun_diff

    def factory_rms(self, fun_diff):
        def fun_rms(*args):
            dy = fun_diff(*args)
            return dy.dot(dy)/dy.size
        return fun_rms

    def _factory_engine(wrap_engine, fun_name):
        def engine(self, p0, *args):
            fun = getattr(self, fun_name)

            if self.do_var is None:
                return wrap_engine(self, fun, p0, *args)
            else:
                do_var = self.do_var.astype(bool)

                p = p0.copy()
                def fun_do_var(p_var, *args):
                    p[do_var] = p_var
                    return fun(p, *args)
            
                p[do_var] = wrap_engine(self, fun_do_var, p0[do_var], *args)
                return p
        return engine

    # wrappers for fitting engines
    def wrap_fmin(self, fun, p0, *args):
        return opt.fmin(fun, p0, args=args, disp=False)

    def wrap_odr(self, fun, p0, x, y, *args):
        return odr.odr(fun, p0, y, x, extra_args=args)[0]

    def wrap_leastsq(self, fun, p0, *args):
        return opt.leastsq(fun, p0, args=args)[0]

    def wrap_leastsq2(self, fun, p0, x, *args):
        return mp.leastsq(fun, (p0.copy(), x, self.buf[:x.size]) + args)

    fit_fmin     = _factory_engine(wrap_fmin, 'fun_rms')
    fit_odr      = _factory_engine(wrap_odr, 'fun')
    fit_leastsq  = _factory_engine(wrap_leastsq, 'fun_diff')
    fit_leastsq2 = _factory_engine(wrap_leastsq2, 'fitfun_diff')

    def fit_custom(self, p0, x, *args, **kw):
        if self.do_var is not None:
            kw.setdefault('do_var', self.do_var)
        return self.custom_engine(p0.copy(), x, self.buf[:x.size], *args, **kw)

    def set_engine(self, engine):
        self.engine = getattr(self, 'fit_' + engine)

    # static
    @memoized_property
    def X(self):
        self.set_norm()
        return self.X

    @memoized_property
    def Y(self):
        self.set_norm()
        return self.Y

    @memoized_property
    def OK(self):
        self.set_OK()
        return self.OK

    @memoized_property
    def P0(self):
        self.set_guess()
        return self.P0

    @memoized_property
    def P(self):
        self.fit()
        return self.P

    @memoized_property
    def p0(self):
        self.set_unnorm()
        return self.p0

    @memoized_property
    def p(self):
        self.set_unnorm()
        return self.p

    def fit(self, P0=None):
        if not self.OK:
            raise FitterError("Cannot fit data that failed is_OK() check")
        if P0 is None:
            P0 = self.P0
        self.P = self.engine(P0, self.X, self.Y, *self.args)

    def eval_guess_norm(self, X):
        return self.fitfun(self.P0, X, *self.args)

    def eval_guess(self, x):
        return self.fitfun(self.p0, x, *self.args)

    def eval_norm(self, X):
        return self.fitfun(self.P, X, *self.args)

    def eval(self, x):
        return self.fitfun(self.p, x, *self.args)

    def get_XY(self):
        if self.OK: 
            return self.X, (self.Y, self.eval_norm(self.X))
        else:
            return self.X, (self.Y,)

    def get_xy(self):
        if self.OK:
            return self.x, (self.y, self.eval(self.x))
        else:
            return self.x, (self.y,)

    def plot(self, ax=None, fun='get_xy', lines=None):
        x, y = getattr(self, fun)()
        sty = (dict(linewidth=1.5),)
       
        ax = get_axes(ax)

        if lines is None: lines = []
        Nl, Ny = len(lines), len(y)

        for li, yi in zip(lines, y):
            li.set_data(x, yi)
            li.set_visible(True)

        for li in lines[Ny:]:
            li.set_visible(False)

        if Nl > 0:
            color = lines[0].get_color()
        else:
            color = ax._get_lines.color_cycle.next()
            sty = (dict(),) + sty

        for yi, styi in zip(y[Nl:], sty):
            li, = ax.plot(x, yi, color=color, **styi)
            lines.append(li)

        return lines

