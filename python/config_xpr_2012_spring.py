######################
# 2012 SPRING CAMPAIGN
######################

tip1 = TipXPR(number=1, pos='lower left', V_keys='ampV', I_keys='ampI1')
tip2 = TipXPR(number=2, pos='lower right', V_keys='ampV', I_keys='ampI2')
tip3 = TipXPR(number=3, pos='upper', V_keys='ampVF', I_keys=None)

head = Head(tips=(tip1, tip2, tip3), R_keys='ampR')

tip3I = TipXPR(number=3, pos='upper', V_keys='ampV', I_keys='ampI3')
headI = Head(tips=(tip1, tip2, tip3I), R_keys='ampR')


############################################
E = campaign.add_experiment(date="20120405")

E.add(27684, "", 
        head = head,
        ampI1 = ampInv*CurrentProbe1[10]*Preamp1[5],
        ampI2 = ampInv*CurrentProbe2[10]*Preamp2[5], **def_LPS)

E.rep(27685, 27684)
E.rep(27686, 27684, "Both current probes on tip 1")

E.rep(27687, 27686, "Changed direction of 2nd current probe",
        ampI2 = CurrentProbe2[10]*Preamp2[5])

E.rep(27688, 27687)
E.rep(27689, 27687, "First shot of Leena's experiment")
E.rep(27690, 27687)

E.rep(27691, 27690, "Current measurement from 10 mA/div to 20 mA/div",
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
        ampI2 = CurrentProbe2[20]*Preamp2[5])

E.rep(27692, 27691, "All the way through, but signals saturated")

E.rep(27693, 27692, "Current probes 50 mA/div",
        ampI1 = ampInv*CurrentProbe1[50]*Preamp1[5],
        ampI2 = CurrentProbe2[50]*Preamp2[5],
        stars = '*****')

E.rep(27694, 27693, "HWM resumed experiment", 
        stars = '*****')

E.rep(27695, 27693, "Calibration after this shot", 
        stars = '*****')


############################################
E = campaign.add_experiment(date="20120621")

# Henrik Mayer
E.add(28232, "304 mm, 1.5 s",
        head = head,
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
        ampI2 = CurrentProbe2[20]*Preamp2[5], **def_LPS)

E.rep(28233, 28232, "440 mm, 5.0 s, disrupted before")
E.rep(28234, 28232, "204 mm, 4.8 s, disrupted before")


############################################
E = campaign.add_experiment(date="20120622")

# Rachael McDermott
E.add(28239, "304 mm, 3.8 s, 20 mA/div", 
        head = head,
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[5],
        ampI2 = CurrentProbe2[20]*Preamp2[5], **def_LPS)

E.rep(28240, 28239, "440 mm, 3.8 s, 50 mA/div",
        ampI1 = ampInv*CurrentProbe1[50]*Preamp1[5],
        ampI2 = CurrentProbe2[50]*Preamp2[5])

E.rep(28241, 28239, "440 mm, 3.8 s, 20 mA/div")
E.rep(28242, 28239, "440 mm, 3.2 s, 20 mA/div")
E.rep(28243, 28239, "440 mm, 3.2 s, 20 mA/div -> saturated I2")

# Tim Happel
E.rep(28244, 28243, "204 mm, 1.5 s, 20 mA/div -> 440 mm acc.")

E.rep(28245, 28244, "204 mm, 1.0 s, 20 mA/div, preamps 2x (was 5x)",
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
        ampI2 = CurrentProbe2[20]*Preamp2[2])

E.rep(28246, 28245, "254 mm, 0.85 s, 20 mA/div, preamps 2x -> L-H transition!",
        descr = "L-mode on way in, ELMing H-mode on way out",
        stars = '*****')

# Hendrik Meyer
E.rep(28250, 28245, "304 mm, 4.0 s, 20 mA/div, preamps 2x")
E.rep(28251, 28245, "304 mm, 1.8 s, 20 mA/div, preamps 2x, (bias voltage less positive)")
E.rep(28252, 28245, "354 mm, 1.75 s, 20 mA/div, preamps 2x -> 1.8 s acc.")
E.rep(28253, 28245, "304 mm, 1.75 s, 20 mA/div, preamps 2x (repeat 28250)")
E.rep(28254, 28245, "304 mm, 1.75 s, 20 mA/div, preamps 2x (repeat 28251)")


############################################
E = campaign.add_experiment(date="20120712")

E.add(28379, "Fixed probe @2564.05 mm", 
        head = head,
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
        ampI2 = CurrentProbe2[20]*Preamp2[2], **def_XPR_LPS)

E.rep(28380, 28379, "Fixed probe @2569.05 mm -> no data")
E.rep(28381, 28379, "-200 V bias -> Kepco breaks in at 0.5 s")
E.rep(28382, 28379, "@2569.10 mm, sweep -> data all the way to 6 s")
E.rep(28383, 28379, "@2589.08 mm, DC offset with small sweep -> worked",
        descr = "Early time base error")


############################################
E = campaign.add_experiment(date="20120713")

E.add(28389, "Fixed probe @2569.05 mm, -180 V with sweep", 
        head = head,
        ampI1 = ampInv*CurrentProbe1[20]*Preamp1[2],
        ampI2 = CurrentProbe2[20]*Preamp2[2], 
        stars = '', **def_XPR_LPS)

E.rep(28390, 28389, "-80 V / 150 V sweep at 100 Hz",
        stars = '')

E.rep(28394, 28389, "",
        stars = '')

E.rep(28395, 28394, "Turn 2nd current probe", 
        ampI1 = CurrentProbe1[20]*Preamp1[2],
        ampI2 = CurrentProbe2[20]*Preamp2[2])

E.rep(28403, 28395, "5 plunges, DC biasing with small rect sweeps")
E.rep(28404, 28395, "repeat")
E.rep(28405, 28395, "4 plunges with decreasing depth")
E.rep(28406, 28395, "4 plunges to 10 cm")
E.rep(28407, 28395, "5 plunges to 10 cm")


############################################
E = campaign.add_experiment(date="20120717")

E.add(28419, "Fcn gen. 20 Vpp, +8 VDC, 0.5 kHz (saturates in Isat regime)", 
        head = head,
        ampI1 = CurrentProbe1[5],
        ampI2 = CurrentProbe2[5], **def_XPR_LPS)

E.rep(28420, 28419, "Fcn gen. 12 Vpp, 4 VDC")
E.rep(28421, 28420)
E.rep(28422, 28421)
E.rep(28423, 28422)
E.rep(28424, 28423, "Asym. waveform, 7 Vpp, 1 VDC, 200 Hz")

E.rep(28425, 28424, "6.9 Vpp, 1 VDC, 100 Hz, 1 mA/div", 
        ampI1 = CurrentProbe1[1],
        ampI2 = CurrentProbe2[1])

E.rep(28426, 28425, "First data! Kepco breaks in in Isat")
E.rep(28427, 28426, "Change sweep pars -> Kepco still breaks")
E.rep(28428, 28427, "Back to other fcn gen. -> saturated at 1 mA/div")

E.rep(28429, 28428, "Back to 5 mA/div -> no data", 
        ampI1 = CurrentProbe1[5],
        ampI2 = CurrentProbe2[5])

E.rep(28434, 28429, "20 mA/div, 0.1 kHz, all 3 tips on bias voltage", 
        head = headI,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20], 
        ampI3 = CurrentProbe3[20], 
        descr = "Plasma died before first plunge")

E.rep(28435, 28434, "0.5 kHz, plunge at 1 s",
        descr = "Strong arc")

E.rep(28436, 28435, "1 plunge at 4 s",
        descr = "Plasma died before plunge")


############################################
E = campaign.add_experiment(date="20120719")

E.add(28442, "0.5 kHz, 3rd pin VF", 
        head = head,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20], **def_XPR_LPS)

E.rep(28444, 28442, "Max plunge at 1.75 s and 3.95 s")
E.rep(28445, 28442, "100 V Kepco")
E.rep(28446, 28442, "200 V Kepco, 12 Vpp, 5 VDC, 0.5 kHz")

# Francois Ryter
E.rep(28448, 28442, "1 kHz, 16 Vpp (digital fcn gen.), VDC from Kepco")
E.rep(28449, 28442, "0.5 kHz, reduced VDC slightly")
E.rep(28450, 28442, "2nd plunge to 1.6 s")
E.rep(28451, 28442, "Max penetration -> shot didn't run")
E.rep(28452, 28442, "Max penetration -> arcs")


############################################
E = campaign.add_experiment(date="20120720")

E.add(28455, "Acquisition with turned-off Kepco", 
        head = head,
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20], **def_XPR_LPS)

E.rep(28466, 28455, "0.5 kHz, 16 Vpp, Kepco offset just avoids saturation")
E.rep(28467, 28455)
E.rep(28468, 28455)

E.rep(28469, 28455, "H: 14 Vpp, max offset on Kepco")
E.rep(28472, 28455, "He: 3 plunges, max penetration", 
        descr = "Data on first plunge, arcs, good comparison between DAQs")

E.rep(28473, 28455, "He again: 1 plunges at 150 mm")


############################################
E = campaign.add_experiment(date="20120726")

E.add(28504, "Calibration, no signals attached", 
        dig='XPR', head=head, amp_default=amp_default_unity, 
        lines=dict(XPR=dict(amp={}, mapping=mapping_XPR), 
                   LPS=dict(amp={}, mapping=mapping_LPS)),
        ampI1 = CurrentProbe1[20],
        ampI2 = CurrentProbe2[20],
        stars = '')

E.rep(28507, 28504, "Calibration, 10 Vpp into 50 Ohm (+/-0.1 A)",
        stars = '')

E.rep(28508, 28504, "Signal also on bias voltage",
        stars = '')


