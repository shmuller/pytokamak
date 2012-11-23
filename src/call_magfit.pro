;	------------------------------------------------------------------
;
;				call_magfit
;				===========
;
;	------------------------------------------------------------------
;
;	INTERFACE:	ruft den C-Fit mag_fit auf
;
;
;	x	x-Vektor, hier Spannung
;	y	zugehöriger abhängiger y-Vektor, hier der Sondenstrom
;	w	Gewichtung der Abweichung, Standardeinstellung = 1.
;
;	yfit	Rückgabe-Vektor, enthält gefittete y-Werte
;
;	params	Parameter, die an Fitroutinen übergeben werden,
;		mit folgenden Bedeutungen:
;
;	params(0)	n_e
;	params(1)	T_e
;	params(2)	psi1
;	params(3)	alpha1
;	params(4)	beta
;	params(5)	d_V
;	params(6)	psi2
;	params(7)	alpha2
;	params(8)	optional, V_plasma - V_wall als Rueckgabewert
;	params(9)	unbenutzt
;	params(10)	optional, p_height einer free standing probe (Länge)
;	params(11)	optional, a_width einer flush mounted probe (Breite)
;	params(12)	optional, a_height einer flush mounted probe (Länge)
;	params(13)	optional, a_mass Massezahl des Füllgases
;	params(14)	optional, a_z Kernladungszahl des Füllgases
;	params(15)	optional, B_t Hauptfeld
;	params(16)	optional, Quotinen T_i / T_e
;
;		==> bei Rückkehr aus Fit, stehen in params die neuen
;		    angepaßten Werte, soweit die in fixed (s.u.) verlangt
;		    wurde
;
;	fixed		welche Parameter sollen variiert werden, welche
;			Parameter sollen nur Werte übergeben
;
;	fixed(i) = 0 	==> Fit in entsprechender Variablen
; 	fixed(i) = 1 	==> keine Fit, lediglich konstanter Parameter
;
;	nb_weights	Gewichte für Fit mit Nebenbedingungen 
;
;	nb_weights(i) = 0 	==> keine Nebenbedingung
;	nb_weights(i) = 1 	==> Nebenbedingung mit entsprechender 
;				    Gewichtung,
;				==> Sollwert wird als Startwert in params
;				    übergeben
;
;	func		welche Fitfunktion soll verwendet werden?
;			bis jetzt implementiert:
;
;	func = fast_doppel(C)	schneller Fit, enthält Näherung für
;				flush mounted probes
;	func = mag_doppel(C)	vollständiger Fit für flush mounted
;				probes, für pin-probes nicht geeignet
;
;	--------------------------------------------------------------------
;
;	spezielle Hinweise für flush mounted probes:
;
;	aus Kompatibilitätsgründen kann params als lediglich 8 elementiger
;	Vektor übergeben werden, die fehlenden Informationen werden dann
;	aus dem common-Block auswertung entnommen. Falls ein vollständiger
;	Vektor übergeben wird, ist darauf zu achten, daß
;		params(10) < 0 und fixed(10) = 1
;	gesetzt wird, um die pin-probes auszuschalten. Wird für params(10),
;	der Länge der pin-probe, ein positiver Wert angegeben, so soll
;			----> in einer späteren Version <-------
;	die Kombination einer pin-probe berechnet werden, deren abschließende
;	Oberfläche nicht parallel zum Magnetfeld angeordnet ist.
;
;	--------------------------------------------------------------------
;
;	spezielle Hinweise für pin probes:
;	
; 	Wahl der Startwerte, Belegung der Variablen
;	params(0) = Schätzung für n_e	fixed(0) = 0
;	params(1) = Schätzung für T_e	fixed(1) = 0
;	params(2) = 0.05		fixed(2) = 1
;	params(3) = 0.17		fixed(3) = 1
;	params(4) = 4.			fixed(4) = 0
;	params(5) = 0.			fixed(5) = 0
;	params(10) = l_Sonde		fixed(10) = 1
;	params(11) = b_sonde		fixed(11) = 1
;	params(12) = -1.		fixed(12) = 1
;	params(13) = Masse in AMU	fixed(13) = 1
;	params(14) = Kernladungszahl	fixed(14) = 1
;
;	alle anderen Prameter können unbelegt beleiben, fixed(*) muß
;	gleich 1 gesetzt werden.
;
;	--------------------------------------------------------------------
;
;	M. Weinlich, 	Max-Planck-Institut für Plasmaphysik
;			mnw@ipp-garching.mpg.de
;			089/3299-1746
;
;	Stand: 28.03.95
;
;	--------------------------------------------------------------------


FUNCTION call_magfit, x, y, w, yfit, params, fixed, nb_weights, func=func




if n_elements(func) ne 1 then func = 'mag_doppel(C)'
func_str = strupcase(strtrim(func,2))
which_func = 0L
if func_str eq 'MAG_DOPPEL(C)' then which_func = 0L $
else if func_str eq 'FAST_DOPPEL(C)' then which_func = 1L $
else which_func = 10L


;
; alle Parameter muessen als double oder long an C uebergeben werden
;

n_x    = long(n_elements(x))

if (size(x))(2) ne 5 then x = double(x)
if (size(y))(2) ne 5 then y = double(y)
if (size(w))(2) ne 5 then w = double(w)
if n_elements(nb_weights) ne n_elements(params) then  $
			nb_weights=dblarr(n_elements(params))
if (size(nb_weights))(2) ne 5 then nb_weights = double(nb_weights)
yfit = dblarr(n_x)

if n_elements(y) ne n_x then stop
if n_elements(w) ne n_x then stop

; Startwerte fuer Parameter in C
c_params = dblarr(20)
c_nb_val = dblarr(20)
c_nb_cst = dblarr(20)
do_var   = lonarr(20)

c_params(0)  = params(0)	; n_e
c_params(1)  = params(1)	; T_e
c_params(2)  = params(3)	; alpha1
c_params(3)  = params(4)	; beta
c_params(4)  = params(5)	; d_V
c_params(5)  = params(7)	; alpha2
c_params(6)  = params(2)	; psi1
c_params(7)  = params(6)	; psi2
if n_elements(a_width) gt 0 then $
	c_params(10) = double(a_width)	; Breite der Sonden
if n_elements(a_height) gt 0 then $
	c_params(11) = double(a_height) ; Laenge der Sonden
if n_elements(a_mass) gt 0 then $
	c_params(12) = double(a_mass)   ; Massezahl der Ionenen in AMU


if n_elements(a_z) lt 1 then begin
   c_params(13) = 1.d
endif else begin
   c_params(13) = double(a_z)      ; Kernladungszahl der Ionenen in AMU
   endelse

if n_elements(a_b_loc) ne 1 then begin
   c_params(14) = 2.d
endif else begin
   c_params(14) = double(a_b_loc)	; Betrag des lokalen Magnetfeldes
   endelse

if n_elements(a_tau) lt 1 then begin
   c_params(16) = 1.d
endif else begin
   c_params(16) = double(a_tau)	; T_i / T_e
   endelse

; optionale Parameter
if (n_elements(params) ge 11 ) then if (abs(params(10)) gt 1.e-8) $
	then c_params(17) = params(10)  ; p_height
if (n_elements(params) ge 12 ) then if (abs(params(11)) gt 1.e-8) $
	then c_params(10) = params(11)  ; a_width
if (n_elements(params) ge 13 ) then if (abs(params(12)) gt 1.e-8) $
	then c_params(11) = params(12)  ; a_height
if (n_elements(params) ge 14 ) then if (abs(params(13)) gt 1.e-8) $
	then c_params(12) = params(13)  ; a_mass
if (n_elements(params) ge 15 ) then if (abs(params(14)) gt 1.e-8) $
	then c_params(13) = params(14)  ; a_z
if (n_elements(params) ge 16 ) then if (abs(params(15)) ge 1.e-8) $
	then c_params(14) = params(15)  ; B_t


; Gewichte fuer Nebenbedingung
; --> eigentliche Gewichte = Gewichtsfaktor / tolerabler Abweichung
;
c_nb_cst(0) = nb_weights(0) / ( alog(params(0)) - alog(0.95*params(0)) )^2
c_nb_cst(1) = nb_weights(1) / ( alog(params(1)) - alog(0.95*params(1)) )^2
c_nb_cst(2) = nb_weights(3) / ( alog(1.0) - alog(1.4) )^2
c_nb_cst(3) = nb_weights(4) / ( alog(params(4)) - alog(params(4)+1.) )^2
c_nb_cst(4) = nb_weights(5) / ( 1.5 )^2
c_nb_cst(5) = nb_weights(7) / ( alog(1.0) - alog(1.4) )^2
c_nb_cst(6) = nb_weights(2) / ( cos(params(2)) - cos(params(2)-0.5*!dtor) )^2
c_nb_cst(7) = nb_weights(6) / ( cos(params(6)) - cos(params(6)-0.5*!dtor) )^2

; Nebenbedingungen: moeglichst nahe am Startwert
for i=0, 9 do begin
   if c_nb_cst(i) gt 0.0 then c_nb_val(i) = c_params(i)
   endfor

; welche Parameter sollen festgehalten werden
do_var(*) = 1L
do_var(0) = fixed(0)
do_var(1) = fixed(1)
do_var(2) = fixed(3)
do_var(3) = fixed(4)
do_var(4) = fixed(5)
do_var(5) = fixed(7)
do_var(6) = fixed(2)
do_var(7) = fixed(6)

; in IDL: 0 --> variiere Parameter,    1 --> Parameter = const.
; in C  : 0 --> Parameter = const.,    1 --> variiere Parameter
ind0 = where( do_var eq 0, c0)
ind1 = where( do_var ne 0, c1)
if c0 gt 0 then do_var(ind0) = 1L
if c1 gt 0 then do_var(ind1) = 0L

n_pars = long(n_elements(c_params))

chi_sqr    = double(0.0)

iter_max = 200L
eps_abs  = 0.d
eps_rel  = 1.d-4


lib = 'mag_fit.so'
routine = 'run_magfit'

; keine Ueberpruefung mathematischer Fehler
junk = check_math(0,1,trap=1)

c_err = call_external( lib, routine, x, y, w, yfit, n_x, c_params, do_var, $
			n_pars, c_nb_val, c_nb_cst, chi_sqr, $
			iter_max, eps_abs, eps_rel, which_func, $
			value=replicate(0B,15) )

; irgendwelche Fehler aufgetreten?
error = check_math(0,0,trap=1) 
if (error ne 0) and (c_err eq 0) then begin
   print, 'math-Error waehrend C-Fit detektiert '
   if error ge 128 then begin
      print, '   ... illegal floating operand'
      error = error - 128
      endif
   if error ge 64 then begin
      print, '   ... illegal floating operand'
      error = error - 64
      endif
   if error ge 32 then begin
      print, '   ... floating overflow'
      error = error - 32
      endif
   if error ge 16 then begin
      print, '   ... floating underflow'
      error = error - 16
      endif
   if error ge 8 then begin
      print, '   ... unknown floating error'
      error = error - 8
      endif
   if error ge 4 then begin
      print, '   ... unknown floating error'
      error = error - 4
      endif
   if error ge 2 then begin
      print, '   ... integer overflow'
      error = error - 2
      endif
   if error ge 1 then begin
      print, '   ... integer divided by zero'
      error = error - 1
      endif
   print,''
   endif

; Rueckschreiben der Parameter
params(0) = c_params(0)
params(1) = c_params(1)
params(3) = c_params(2)
params(4) = c_params(3)
params(5) = c_params(4)
params(7) = c_params(5)
params(2) = c_params(6)
params(6) = c_params(7)

; wenn genuegend Platz vorhanden, dann soll auch V_wand - V_plasma
; abgespeichert werden
if n_elements(params) gt 8 then begin
   params(8) = c_params(15)
   endif


return, c_err

end
