import config_xpr
import probe_xpr

for expt in config_xpr.campaign:
    print ""
    print "Experiment from: %s" % expt.date
    print "-------------------------"
    print ""

    for shot in expt:
        if len(shot.stars) == 0:
            continue

        print "=========================================================="
        print shot
        print ""
        print shot.descr
        print ""

        XPR = probe_xpr.ProbeXPR(shn=shot.shn)
        
        tM, RM = XPR.get_dwell_params()
        for par in zip(range(len(tM)), tM, RM):
            print "Plunge %d: tM = %.02f s, RM = %.02f m" % par

        print ""


