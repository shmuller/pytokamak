from sm_pyplot.tight_figure import get_fig, show

from mdsclient import *

shn = 30017
diag = 'EQI'
i = 20

sock = mdsconnect('localhost:8001')

mdsfmt = 'augsignal({shn}, "{diag}", "%s", "AUGD")'.format(shn=shn, diag=diag)

t = mdsvalue(sock, mdsfmt % 'time')
R = mdsvalue(sock, mdsfmt % 'Ri')[:t.shape[0]]
z = mdsvalue(sock, mdsfmt % 'Zj')[:t.shape[0]]
psi = mdsvalue(sock, mdsfmt % 'PFM')[:t.shape[0],:z.shape[1],:R.shape[1]]


labi = ('Magn. axis', 'X-point', 'Inner limiter', '2nd X-point', 'Outer limiter')
cati = mdsvalue(sock, mdsfmt % 'ikCAT')
psii = mdsvalue(sock, mdsfmt % 'PFxx')
Ri = mdsvalue(sock, mdsfmt % 'RPFx')
zi = mdsvalue(sock, mdsfmt % 'zPFx')

mdsdisconnect(sock)

fig = get_fig()
ax = fig.axes[0]
ax.axis('equal')

ax.contour(R[i], z[i], psi[i], 20)
ax.plot(Ri[i], zi[i], '*')

ax.contour(R[i], z[i], psi[i], psii[i])

#fig.colorbar(ax._gci())

show()

