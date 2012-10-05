import numpy as np

import scipy.interpolate as interp
import scipy.optimize as opt
import scipy.odr as odr

import logging
reload(logging)
logging.basicConfig(level=logging.WARN)
logger = logging

from pdb import set_trace

from sig import *

class ArrayView(np.ndarray):
    def __new__(subtype, x, fields):
        dtype = {f: x.dtype.fields[f] for f in fields}
        return np.ndarray.__new__(subtype, x.shape, dtype,
                                  buffer=x, strides=x.strides)


class MouseMotionObserverViewer:
    def __init__(self, observers, viewer, plotfun):
        self.observers, self.viewer = observers, viewer
        self.plotfun = plotfun

        self.observer_canvas = self.observers[0].figure.canvas
        self.viewer_canvas = self.viewer.figure.canvas
        
        self.cid = self.observer_canvas.mpl_connect('motion_notify_event', self.on_move)
        self.viewer_canvas.mpl_connect('resize_event', self.on_resize)

        self.on_resize(None)

    def on_resize(self, event):
        logger.debug("Resizing")
        save_vis = [line.get_visible() for line in self.viewer.lines]
        for line in self.viewer.lines:
            line.set_visible(False)
        
        self.viewer_canvas.draw()
        self.background = self.viewer_canvas.copy_from_bbox(self.viewer.bbox)
        
        for line, vis in zip(self.viewer.lines, save_vis):
            line.set_visible(vis)
            self.viewer.draw_artist(line)

    def on_move(self, event):
        logger.debug("Moving")
        if event.inaxes in self.observers:
            self.viewer_canvas.restore_region(self.background)
                
            self.plotfun(event)

            for line in self.viewer.lines:
                self.viewer.draw_artist(line)
            self.viewer_canvas.blit(self.viewer.bbox)

    def __del__(self):
        logger.debug("Destroying")
        self.observer_canvas.mpl_disconnect(self.cid)


class PiecewiseLinear:
    def __init__(self, y, x):
        self.y, self.x = y, x
        self.N = self.x.shape[0]
        self.shape = self.y.shape[2:]

    def __getitem__(self, index):
        if not isinstance(index, tuple): index = (index,)
        index = (slice(None), slice(None)) + index
        return PiecewiseLinear(self.y[index], self.x)

    def plot(self, ax=None):
        shape = (self.N*2,)
        x = self.x.reshape(shape)
        y = self.y.reshape(shape + self.shape)
        ax = get_axes(ax)
        return ax.plot(x, y)

    def plot2(self, ax=None):
        shape = (self.N, 1)
        nanx = np.empty(shape) + np.nan
        nany = np.empty(shape + self.shape) + np.nan

        shape = (self.N*3,)
        x = np.concatenate((self.x, nanx), 1).reshape(shape)
        y = np.concatenate((self.y, nany), 1).reshape(shape + self.shape)
        ax = get_axes(ax)
        return ax.plot(x, y)


class PiecewisePolynomial:
    def __init__(self, c, x, **kw):
        self.kw = {'fill': None, 'i0': np.arange(x.size), 'i1': None}
        self.kw.update(kw)
        for key, val in self.kw.iteritems():
            setattr(self, key, val)

        self.c, self.x = c, x
        self.N = self.i0.size
        self.shape = self.c.shape[2:]

    def __getitem__(self, index):
        if not isinstance(index, tuple): index = (index,)
        index = (slice(None), slice(None)) + index
        return PiecewisePolynomial(self.c[index], self.x, **self.kw)

    def __call__(self, X, side='right'):
        xi = self.x[self.i0]
        ind = np.searchsorted(xi, X, side) - 1
        outl, outr = ind < 0, ind > self.N-2
        ind[outl], ind[outr] = 0, self.N-2
        
        dX = X - xi[ind]

        #Y = reduce(lambda Y, c: Y*dX + c[ind], self.c)
        
        c = self.c
        Y = c[0,ind].copy()
        for a in c[1:]:
            Y = Y*dX + a[ind]

        if self.fill is not None:
            Y[outl | outr] = self.fill
        return Y

    @property
    def T(self):
        return PiecewisePolynomial(self.c.swapaxes(2,3), self.x, **self.kw)

    @property
    def savefields(self):
        return DictView(self.__dict__, ('c', 'x', 'i0'))

    def _mask(self, w):
        ind0, ind1 = np.searchsorted(self.i0, w)
        return np.concatenate(map(np.arange, ind0, ind1))

    @staticmethod
    def cat(a, axis=0):
        return a[0].__array_wrap__(np.concatenate(a, axis))

    def eval(self, x=None, w=None):
        i0 = self.i0
        if w is not None:
            i0 = i0[self._mask(w)]

        xi = self.x[i0]
        xl, xr = xi[:-1], xi[1:]
        yl, yr = self(xl, 'right'), self(xr, 'left')

        if x is not None:
            xi = x[i0]
            xl, xr = xi[:-1], xi[1:]

        shape = (xl.size + xr.size,)
        x = self.cat((xl[:,None], xr[:,None]), 1).reshape(shape)
        y = self.cat((yl[:,None], yr[:,None]), 1).reshape(shape + self.shape)
        return x, y

    def plot(self, ax=None, x=None, w=None):
        x, y = self.eval(x=x, w=w)
                
        ax = get_axes(ax)
        return ax.plot(x, y)


def wrap_fmin(fun, p0, x, y):
    def dy2(p):
        dy = fun(p, x) - y
        return dy.dot(dy)/dy.size

    return opt.fmin(dy2, p0, disp=False)

def wrap_odr(fun, p0, x, y):
    mod = odr.Model(fun)
    dat = odr.Data(x, y)
    o = odr.ODR(dat, mod, p0)
    out = o.run()
    return out.beta

engines = {'fmin': wrap_fmin, 'odr': wrap_odr}


class Fitter:
    def __init__(self, x, y, engine='fmin'):
        self.x, self.y = x, y

        self.OK = None
        self.X = self.Y = None
        self.P = self.P0 = None
        self.p = self.p0 = None

        self.set_engine(engine)

    # overload
    def set_OK(self):
        self.OK = True
        return self.OK

    def set_norm(self):
        self.X, self.Y = self.x, self.y
        return self.X, self.Y

    def set_unnorm(self):
        self.p, self.p0 = self.P, self.P0
        return self.p, self.p0

    def set_guess(self):
        self.P0 = 0.
        return self.P0

    def fitfun(self, P, X):
        pass

    # static
    def is_OK(self):
        if self.OK is None:
            self.set_OK()
        return self.OK

    def get_norm(self):
        if self.X is None:
            self.set_norm()
        return self.X, self.Y

    def get_unnorm(self):
        if self.p is None:
            self.set_unnorm()
        return self.p, self.p0

    def get_guess(self):
        if self.P0 is None:
            self.set_guess()
        return self.P0
    
    def set_engine(self, engine):
        self.engine = engines[engine]

    def fit(self, P0=None):
        if not self.is_OK():
            raise RuntimeError("Cannot fit data that failed is_OK() check")
        
        if P0 is None:
            P0 = self.get_guess()
        self.P = self.engine(self.fitfun, P0, self.X, self.Y)
        self.set_unnorm()

        return self.p

    def eval_guess_norm(self, X):
        return self.fitfun(self.P0, X)

    def eval_guess(self, x):
        return self.fitfun(self.p0, x)

    def eval_norm(self, X):
        return self.fitfun(self.P, X)

    def eval(self, x):
        return self.fitfun(self.p, x)

    def get_XY(self):
        if self.is_OK(): 
            if self.P is None:
                self.fit()
            return self.X, (self.Y, self.eval_norm(self.X))
        else:
            return self.X, (self.Y,)

    def get_xy(self):
        if self.is_OK(): 
            if self.p is None:
                self.fit()
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


def LP_unnormalize(P, Vm, VM, Im, IM):
    dV, dI = VM-Vm, IM-Im
    a = 2./dI
    b = -(Im*a+1)

    p = np.empty_like(P)
    p[0] = (P[0]-b)/a
    p[1] = (P[1]+P[2]*np.log(1.-b/P[0]))*dV+Vm
    p[2] = P[2]*dV
    return p

def LP_normalize(p, Vm, VM, Im, IM):
    dV, dI = VM-Vm, IM-Im
    a = 2./dI
    b = -(Im*a+1)

    P = np.empty_like(p)
    P[0] = a*p[0]+b
    P[1] = (p[1]-p[2]*np.log(1.-b/P[0])-Vm)/dV
    P[2] = p[2]/dV
    return P


class FitterIV(Fitter):
    def __init__(self, V, I, mask=None, cut_at_min=True, **kw):
        self.mask_ind = np.arange(V.size)
        if mask is not None:
            self.mask_ind = self.mask_ind[mask]

        self.sort_ind = V[self.mask_ind].argsort()
        self.ind = self.mask_ind[self.sort_ind]

        self.V, self.I = V[self.ind], I[self.ind]

        self.im = self.I.argmin()
        self.Vm, self.VM = self.V[0], self.V[-1]
        self.Im, self.IM = self.I[self.im], np.median(self.I[:self.I.size/2])
        self.dV = self.VM - self.Vm
        self.dI = self.IM - self.Im

        self.cut_at_min = cut_at_min
        self.M = self.V.size

        Fitter.__init__(self, self.V, self.I, **kw)

    def get_ind(self):
        return self.ind[:self.M]

    def set_OK(self):
        def medianstd(x):
            xm = np.median(x)
            dx = x-xm
            xs = np.sqrt(dx.dot(dx)/dx.size)
            return xm, xs

        N = self.I.size
        Ilm, Ils = medianstd(self.I[:N/10])
        Irm, Irs = medianstd(self.I[N*9/10:])

        cnd1 = Ilm > Ils
        cnd2 = Irm < -2*Ils
        self.OK = cnd1 & cnd2
        return self.OK

    def set_norm(self):
        if self.cut_at_min:
            self.M = self.im+1
        self.X = (self.V[:self.M].astype('d') - self.Vm)/self.dV
        self.Y = (self.I[:self.M].astype('d') - self.Im)/self.dI*2 - 1
        return self.X, self.Y

    def set_unnorm(self):
        self.p0 = LP_unnormalize(self.P0, self.Vm, self.VM, self.Im, self.IM)
        self.p  = LP_unnormalize(self.P , self.Vm, self.VM, self.Im, self.IM)
        return self.p, self.p0

    def set_guess(self):
        def find_i0(x):
            return np.flatnonzero(np.diff(np.sign(x)))[-1]

        V, I = self.get_norm()
        I0 = 1.
        i0 = find_i0(I)
        Vf = (V[i0]+V[i0+1])/2
        Te = (V[self.im]-Vf)/np.log(1-I[self.im]/I0)

        self.P0 = np.array([I0, Vf, Te])
        return self.P0

    def fitfun(self, P, X):
        return P[0]*(1.-np.exp((X-P[1])/P[2]))

    def fit(self):
        Fitter.fit(self)
        save = self.X, self.Y
        Y0 = 0.
        while True:
            self.M *= 0.95
            self.X, self.Y = self.X[:self.M], self.Y[:self.M]
            P_old = self.P
            Fitter.fit(self, P0=self.P)
            
            if np.any(self.P > P_old) or (self.eval_norm(self.X[-1]) > Y0):
                self.P = P_old
                break

        self.X, self.Y = save
        self.set_unnorm()
        return self.p


class IVChar:
    def __init__(self, V, I, mask=None, **kw):
        self.V, self.I = V, I
        self.fitter_IV = FitterIV(self.V.x, self.I.x, **kw)

    def fit(self, out=None):
        if out is None:
            out = np.empty(3)
        try:
            out[:] = self.fitter_IV.fit()
        except RuntimeError:
            out[:] = np.nan
        return out
    
    def mask(self, out=None):
        if out is None:
            out = np.empty(self.I.size, bool)
        out.fill(False)
        if self.fitter_IV.p is not None:
            out[self.fitter_IV.get_ind()] = True
        return out

    def plot(self, ax=None, fun='get_xy', lines=None):
        return self.fitter_IV.plot(ax, fun, lines)


class IVGroup:
    def __init__(self, V, II, s=slice(None), **kw):
        self.IV_char = np.empty(len(II), object)
        for j, I in enumerate(II):
            self.IV_char[j] = IVChar(V[s], I[s], **kw)

    def __getitem__(self, index):
        return self.IV_char[index]

    def fit(self, out=None):
        if out is None:
            out = np.empty((self.IV_char.size, 3))
        for p, IV_char in zip(out, self.IV_char):
            IV_char.fit(p)
        return out

    def mask(self, out=None):
        if out is None:
            out = np.empty((self.IV_char.size, self.IV_char[0].I.size), bool)
        for I, IV_char in zip(out, self.IV_char):
            IV_char.mask(I)
        return out

    def plot(self, ax=None, fun='get_xy'):
        ax = get_axes(ax)
        cache = getattr(ax, 'lines_cache', None)
        if cache is None:
            ax.lines_cache = [x.plot(ax, fun) for x in self.IV_char]
        else:
            ax.lines_cache = [x.plot(ax, fun, l) for x, l in zip(self.IV_char, cache)]


class IVSeries:
    def __init__(self, V, II, iE, **kw):       
        self.V, self.II, self.iE = V, II, iE
        N = iE.shape[0]
        self.siz = (N, len(II))
        self.ti = V.t[iE]

        self.V_range = V.plot_range()
        self.I_range = self._plot_range(II)

        self.IV_group = np.empty(N, object)
        for j in xrange(N):
            s = self._slice(j)
            self.IV_group[j] = IVGroup(V, II, s, **kw)

    def __getitem__(self, index):
        return self.IV_group[index]

    def _plot_range(self, II):
        I_range = np.array([I.plot_range() for I in II])
        return I_range[:,0].min(), I_range[:,1].max()

    def _slice(self, j):
        return slice(self.iE[j,0], self.iE[j,1]+1)

    def fit(self):
        out = np.empty(self.siz + (3,))
        for p, IV_group in zip(out, self.IV_group):
            IV_group.fit(p)

        i0 = np.r_[self.iE[:,0], self.iE[-1,1]]
        self.PP = PiecewisePolynomial(out[None], self.V.t, i0=i0)
        return self.PP

    def mask(self):
        out = np.zeros((len(self.II), self.II[0].size), bool)
        for j, IV_group in enumerate(self.IV_group):
            s = self._slice(j)
            IV_group.mask(out[:,s])
        return out

    def eval(self):
        V = self.V
        PP = self.PP(V.t)
        fitfun = self[0][0].fitter_IV.fitfun
        out = fitfun(PP.T, V.x)
        out[~self.mask()] = np.nan
        return out

    def plot(self):
        fig = figure()
        n = len(self.II)
        IIfit = self.eval()
        for i, I, Ifit in zip(xrange(n), self.II, IIfit):
            ax = fig.add_subplot(n, 1, 1+i)
            ax.plot(I.t, I.x, I.t, Ifit)
            ax.grid(True)
            ax.set_ylabel(I.name)
        ax.set_xlabel("t [s]")
        
    def animate(self, fun='get_xy'):
        ion()
        fig = figure()
        ax = fig.gca()
        if fun == 'get_xy':
            ax.set_xlim(self.V_range)
            ax.set_ylim(self.I_range)
        else:
            ax.set_xlim(( 0.0, 1.0))
            ax.set_ylim((-1.2, 1.2))
        
        canvas = fig.canvas
        canvas.draw()
        background = canvas.copy_from_bbox(ax.bbox)

        for IV_group in self.IV_group:
            canvas.restore_region(background)
            IV_group.plot(ax, fun)
            for line in ax.lines:
                ax.draw_artist(line)
            canvas.blit(ax.bbox)
            

class IVSeriesViewer:
    def __init__(self, IV_series):
        self.IV_series = IV_series
        self.MMOV = None

    def on_close(self, event):
        logger.debug("Closing")
        self.MMOV = None

    def plotfun(self, event):
        t_event = event.xdata
        ti = self.IV_series.ti[:,1]
        i = min(np.searchsorted(ti, t_event), ti.size-1)
        self.IV_series.IV_group[i].plot(self.ax, 'get_xy')

    def toggle(self, event):
        fig = event.inaxes.figure
        if self.MMOV is not None:
            self.MMOV.viewer_canvas.manager.destroy()
            return
            
        fig2 = figure(figsize=(6,5))
        self.ax = fig2.add_subplot(1, 1, 1)
        V_range, I_range = self.IV_series.V_range, self.IV_series.I_range
        #corners = (V_range[0], I_range[0]), (V_range[1], I_range[1])
        #ax.dataLim.update_from_data_xy(corners, ignore=True)
        self.ax.set_xlim(V_range)
        self.ax.set_ylim(I_range)
        self.ax.grid(True)

        fig2.canvas.mpl_connect('close_event', self.on_close)

        self.MMOV = MouseMotionObserverViewer(fig.axes, self.ax, self.plotfun)
        fig2.show()


class FitterIV2(Fitter):
    def __init__(self, V, I, mask, p0, **kw):
        self.V, self.I = V[mask], I[mask]
        self.Vm, self.VM = self.V.min(), self.V.max()
        self.Im, self.IM = self.I.min(), self.I.max()
        self.dV = self.VM - self.Vm
        self.dI = self.IM - self.Im

    def set_norm(self):
        self.X = (self.V.astype('d') - self.Vm)/self.dV
        self.Y = (self.I.astype('d') - self.Im)/self.dI*2 - 1
        return self.X, self.Y

    def set_unnorm(self):
        self.p = LP_unnormalize(self.P , self.Vm, self.VM, self.Im, self.IM)
        return self.p


class IVSeries2:
    def __init__(self, V, II, PP, mask, iE=None, **kw):
        self.V, self.II, self.PP, self.mask, self.iE = V, II, PP, mask, iE

    def plot(self, i, j):
        s = slice(self.iE[i], self.iE[j]+1)
        V = self.V[s]
        I = self.II[0][s]
        pp = self.PP[s, 0]
        m = self.mask[0, s]

        t = V.t[m]
        V = V.x[m]
        I = I.x[m]
        pm = pp[m]

        print pp.shape

        def fitfun(p, V): 
            return p[0]*(1.-np.exp((V-p[1])/p[2]))

        a = (t-t[0])/(t[-1]-t[0])

        def fitfun2(p, V):
            Is = p[0] + a*(p[3]-p[0])
            Vf = p[1] + a*(p[4]-p[1])
            Te = p[2] + a*(p[5]-p[2])
            return Is*(1.-np.exp((V-Vf)/Te))

        Ifit = fitfun(pm.T, V)

        """
        p0 = pm.mean(0)
        p = wrap_fmin(fitfun, p0, V, I)
        Ifit2 = fitfun(p, V)
        """

        p0 = np.r_[pm[0], pm[-1]]
        p = wrap_fmin(fitfun2, p0, V, I)
        Ifit2 = fitfun2(p, V)

        ax = get_axes()
        ax.plot(t, I, t, Ifit, t, Ifit2)
        return p


class PhysicalResults:
    def __init__(self, shn, R, i0, meas):
        self.shn, self.R = shn, R

        self.keys = ('n_cs', 'Mach', 'nv', 'j', 'Vf', 'Te', 'cs', 'n', 'v', 'pe')

        self.units = dict(n_cs='m$^{-2}$ s$^{-1}$', Mach=None, nv='m$^{-2}$ s$^{-1}$', 
                j='kA m$^{-2}$', Vf='V', Te='eV', cs='km s$^{-1}$', 
                n='m$^{-3}$', v='km s$^{-1}$', pe='Pa')

        self.fact = dict.fromkeys(self.keys, 1)
        self.fact['j'] = self.fact['cs'] = self.fact['v'] = 1e-3

        res = self.calc_res(meas)
        
        self.PP = PiecewisePolynomial(res[None], R.t, i0=i0)

    def calc_res(self, meas):
        qe = 1.6022e-19
        mi = 2*1.67e-27

        Gp = meas.jp/qe
        Gm = meas.jm/qe

        dtype = zip(self.keys, [np.double]*len(self.keys))
        res = np.empty(Gp.size, dtype).view(np.recarray)

        res.Mach = Mach = 0.5*np.log(Gp/Gm)
        res.n_cs = n_cs = np.e*np.sqrt(Gp*Gm)
        res.nv = nv = n_cs*Mach
        res.j  = qe*nv
        res.Vf = Vf = meas.Vf
        res.Te = Te = meas.Te

        Ti = Te
        res.cs = cs = np.sqrt(qe/mi*(Te+Ti))
        res.n  = n = n_cs/cs
        res.v  = v = Mach*cs
        res.pe = n*qe*Te
        return res

    def plot_key(self, key, x, y, xlab, ax=None):
        ax = get_axes(ax, figure=tfigure)

        ylab = key
        if self.units[key] is not None:
            ylab += ' [' + self.units[key] + ']'
        ax.set_ylabel(ylab)

        ax.set_xlabel(xlab)
        ax.grid(True)
        ax.plot(x, self.fact[key]*y[key], label=self.shn)

    def plot(self, fig=None, keys=None, x=None, plunge=None, inout=None):
        if x is None:
            xlab = 't [s]'
        elif x == 'R':
            x = self.R.x
            xlab = "R [mm]"

        w = self.R.plunges(plunge, inout)

        x, y = self.PP.eval(x=x, w=w)

        if keys is None:
            keys = ('n', 'Mach'), ('Te', 'j')
        keys = np.array(keys, ndmin=2)

        if fig is None:
            fig = tfigure()
            gs = gridspec.GridSpec(*keys.shape)
            for i in xrange(keys.size):
                fig.add_subplot(gs[i])

        ax = fig.axes
        for i in xrange(keys.size):
            self.plot_key(keys.flat[i], x, y, xlab, ax=ax[i])

        fig.axes[0].legend(loc='upper left')
        fig.canvas.draw()
        return fig

    def plot_R(self, **kw):
        return self.plot(x='R', **kw)


class LoadResultsError(Exception):
    pass


class Probe:
    def __init__(self, digitizer=None):
        self.digitizer = digitizer
        self.PP = self.IV_series = self.S = None

        self.xlab = "t [s]"
        self.ylab = ("Isat [A]", "Vf [V]", "Te [eV]")

    def __getitem__(self, index):
        if not isinstance(index, tuple):
            return self.S[index]
        else:
            return tuple([self.S[i] for i in index])

    def mapsig(self):
        pass

    def calib(self):
        pass
        
    def load_raw(self, loadfun='load', plunge=None, calib=True, corr_capa=False):
        self.x = getattr(self.digitizer, loadfun)()

        self.mapsig()
        if calib:
            self.calib()
        if corr_capa:
            self.corr_capa()
        if plunge is not None:
            self.trim(plunge)
        
    def get_type(self, type):
        istype = lambda x: x.type == type
        return filter(istype, self.S.itervalues())

    def plot_raw(self, fig=None, **kw):
        if self.S is None:
            self.load_raw(**kw)

        types = ['Current', 'Voltage', 'Position']
        if fig is None:
            fig = figure()
            for i, typ in enumerate(types):
                ax = fig.add_subplot(3, 1, 1+i)
                ax.grid(True)
                ax.set_ylabel(typ)
            ax.set_xlabel(self.xlab)
        
        for typ, ax in zip(types, fig.axes):
            for S in self.get_type(typ):
                S.plot(ax)

        fig.canvas.draw()
        return fig
        
    def trim(self, plunge='all'):
        S = self.get_type('Position')
        i0, iM, i1 = S[0].t_ind

        if plunge == 'all':
            s = np.concatenate(map(np.arange, i0, i1))
        else:
            s = slice(i0[plunge], i1[plunge])

        for S in self.S.itervalues():
            S.trim(s)

    def smooth_I(self, w=10):
        for I in self.get_type('Current'):
            I.smooth(w)

    def corr_capa(self):
        for I in self.get_type('Current'):
            I.x[:] -= I.I_capa()

    def calc_IV_series(self, n=1, **kw):
        II = self.get_type('Current')
        V = II[0].V
        if V.iE is None:
            V.chop_sweeps()
        iE = np.c_[V.iE[:-n], V.iE[n:]]
        return IVSeries(V, II, iE, **kw)

    def analyze(self, **kw):
        if self.S is None:
            self.load_raw(**kw)
        self.IV_series = self.calc_IV_series(engine='fmin')
        self.PP = self.IV_series.fit()

    def results(self, **kw):
        if self.PP is None:
            try:
                self.load(**kw)
            except LoadResultsError:
                self.analyze()
                self.save_res()

        tips = self.config.head.tips

        Isat = self.PP.c[0,:,:,0].T
        Vf   = self.PP.c[0,:,:,1].T
        Te   = self.PP.c[0,:,:,2].T

        keys = ('jp', 'jm', 'Vf', 'Te')
        dtype = zip(keys, [np.double]*len(keys))
        meas = np.empty(Isat.shape[1], dtype).view(np.recarray)

        meas.jp = Isat[0]/tips[0].proj_area
        meas.jm = Isat[1]/tips[0].proj_area
        meas.Vf = Vf[-1]
        meas.Te = Te[-1]

        shn = self.digitizer.shn
        self.res = PhysicalResults(shn, self['R'], self.PP.i0, meas)
        return self.res
        
    @property
    def h5name_res(self):
        return self.digitizer.IO_file.h5name[:-3] + "_res.h5"

    def save_res(self):
        IO = IOH5(self.h5name_res)
        d = self.PP.savefields
        IO.save(d)

    def load_res(self):
        IO = IOH5(self.h5name_res)
        try:
            d = IO.load()
        except IOError:
            raise LoadResultsError("no results available")
        self.PP = PiecewisePolynomial(**d)

    def load(self, **kw):
        self.load_raw(**kw)
        self.load_res()

    def analyze2(self, n=2, **kw):
        if self.IV_series is None:
            self.analyze(**kw)

        V = self.S['V']
        II = self.get_type('Current')
        self.IV_series2 = self.calc_IV_series(n=2, engine='fmin')
        self.PP2 = self.IV_series2.fit()

        """
        V = self.S['V']
        II = self.get_type('Current')
        PP = self.PP(V.t)
        mask = self.IV_series.mask()
        
        iE = self.IV_series.iE
        self.IV_series2 = IVSeries2(V, II, PP, mask, iE)
        """

    def plot(self, fig=None, PP='PP', x=None, plunge=None, inout=None):
        if self.PP is None:
            self.analyze()

        if PP == 'PP':
            ylab = self.ylab
        else:
            ylab = ("n [au], M",) + self.ylab[1:]
        PP = getattr(self, PP)

        if x is None:
            xlab = self.xlab
        elif x == 'R':
            x = self['R'].x
            xlab = "R [mm]"

        w = self['R'].plunges(plunge, inout)

        if fig is None:
            IV_series_viewer = IVSeriesViewer(self.IV_series)
            menu_entries_ax = (('IV viewer', IV_series_viewer.toggle),)

            fig = figure(figsize=(10,10), menu_entries_ax=menu_entries_ax)
            for i in xrange(3):
                ax = fig.add_subplot(3, 1, 1+i)
                ax.grid(True)
                ax.set_ylabel(ylab[i])
            ax.set_xlabel(xlab)
            fig.IV_series_viewer = IV_series_viewer

        for pp, ax in zip(PP.T, fig.axes):
            pp.plot(ax, x=x, w=w)

        fig.canvas.draw()
        return fig

    def plot_R(self, **kw):
        return self.plot(x='R', **kw)

    


