import mdsclient
import EFIT

from sm_pyplot.tight_figure import get_fig, show

sock = mdsclient.mdsconnect('localhost:8020')

E = EFIT.EFIT(sock, 141451, ti=1650.)

E.load()

FS = E.getFluxSurf()

fig = get_fig()
FS.plot()
show()



