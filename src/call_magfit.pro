PRO check_lib


common		mag_so_aktuell,		was_checked, lib_file, interface


; Initialisiserung wenn noetig
if n_elements(was_checked) lt 1 then was_checked = 0 eq 1 else return

print
print, 'mag_fit.so: verfuegbar unter Solaris, SunOS und AIX'
print
 
; keine Unterscheidungen nach Betriebssystem mehr noetig, wird hier nur
; belassen, um evtl. Testversionen leichter einspielen zu koennen  
hv = strupcase(!version.os)
lib_file = '/afs/ipp/home/l/lang/bin/mag_fit.so'
if strpos( hv, 'AIX', 0) ge 0 then begin
   ;lib_file = '/afs/ipp/u/mnw/bin/@sys/mag_fit.so'
endif else if strpos( hv, 'SUNOS', 0) ge 0 then begin
   ;lib_file = '/afs/ipp/u/mnw/bin/@sys/mag_fit.so'
endif else begin
   print
   print, 'unknown os-version ....'
   print
   stop
   endelse

interface = 'idl2mag'      

rc = findfile( lib_file, count=count1)
if count1 lt 1 then begin
   print, 'check_lib: kann ' + lib_file + ' nicht finden!'
   stop
   endif

real_lib = follow_link( lib_file)

was_checked = 0 eq 0
print, 'using: ' + interface 
print, '       aus ' + lib_file
if strlen(real_lib) gt 1 then print, '       --> ' + real_lib
print

if strlen(real_lib) gt 1 then lib_file = real_lib



end	; check_lib





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
;	params(10)	p_height einer free standing probe (Länge)
;	params(11)	a_width einer flush mounted probe (Breite)
;	params(12)	a_height einer flush mounted probe (Länge)
;	params(13)	a_mass Massezahl des Füllgases
;	params(14)	a_z Kernladungszahl des Füllgases
;	params(15)	B_t Hauptfeld
;	params(16)      Quotient T_i / T_e
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
;			weinlich@ipp-garching.mpg.de
;			089/3299-1813
;
;	Stand: 28.03.95
;
;	28.07.95	Anpassung an plattform-uebergreifenden Betrieb 
;			auf SunOs4.x und AIX3.2x
;
;	--------------------------------------------------------------------


FUNCTION call_magfit, x, y, yfit, params, fixed


;@auswertung.common

common		mag_so_aktuell,		was_checked, lib_file, interface
common		mag_fit_errors,		n_err_report, last_err


; Initialisiserung wenn noetig
check_lib


routine = 'mag_doppel' 

   

; ------------------------------------------------------------------------
;
; alle Parameter muessen als double oder long an C uebergeben werden
;

if n_elements(x) ne n_elements(y) then begin
   hnx = n_elements(x)
   hny = n_elements(y)
   hn  = hnx < hny
   message, 'x- und y-Vektoren haben unterschiedliche Laenge!', /info
   message, 'n(x) = ' + strtrim(hnx,2) + ',  n(y) = ' + strtrim(hny,2), /info
   message, '--> Reduktion beider Vektoren auf ' + strtrim(hn,2) + ' Punkte', /info
   x = x(0:hn-1)
   y = y(0:hn-1)
   endif
   
n_x  = long(n_elements(x))
yfit = dblarr(n_x, /nozero)

if (size(x))(2) ne 5 then x = double(x)
if (size(y))(2) ne 5 then y = double(y)
w = y*0+1.

nb_weights=dblarr(n_elements(params))


; Startwerte fuer Parameter in C
c_params = dblarr(20)
c_nb_val = dblarr(20)		; Nebenbedingung, ausgeschaltet
c_nb_cst = dblarr(20)		; Nebenbedingung, ausgeschaltet
do_var   = lonarr(20, /nozero)

c_params(0)  = params(0)	; n_e
c_params(1)  = params(1)	; T_e
c_params(2)  = params(3)	; alpha1
c_params(3)  = params(4)	; beta
c_params(4)  = params(5)	; d_V
c_params(5)  = params(7)	; alpha2
c_params(6)  = params(2)	; psi1
c_params(7)  = params(6)	; psi2

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



; optionale Parameter
if n_elements(params) gt 10 then begin

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

   endif


; in IDL: 0 --> variiere Parameter,    1 --> Parameter = const.
; in C  : 0 --> Parameter = const.,    1 --> variiere Parameter
ind0 = where( do_var eq 0, c0)
ind1 = where( do_var ne 0, c1)
if c0 gt 0 then do_var(ind0) = 1L
if c1 gt 0 then do_var(ind1) = 0L

n_pars = n_elements(c_params)


; --------------------------------------------------------------------------
;
; fest vorgegebene Einstellung fuer Fitgenauigkeit, Iterationstiefe etc.
;
chi_sqr    = double(0.0)
iter_max = 200L
eps_abs  = 0.d
eps_rel  = 1.d-4

; keine Ueberpruefung mathematischer Fehler
junk = check_math(0,1,trap=1)
			
c_err = call_external( lib_file, interface, routine, x, y, w, yfit, n_x,  $
			c_params, do_var, n_pars, c_nb_val, c_nb_cst, $
			chi_sqr, iter_max, eps_abs, eps_rel,  $
			value=replicate(0B,15) )

; irgendwelche Fehler aufgetreten?
error = check_math(0,0,trap=1) 
if n_elements(n_err_report) ne 1 then n_err_report = 0L
n_err_max = 10L
if (n_err_report gt n_err_max) then begin
   if (error ne last_err) then n_err_report = n_err_max-1
   endif
if (error ne 0) and (c_err eq 0) and (n_err_report lt n_err_max) then begin
   n_err_report = n_err_report+1L
   last_err = error
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





;	------------------------------------------------------------------
;
;				call_mag_solver
;				==================
;
;	------------------------------------------------------------------
;
;	INTERFACE:	ruft den nichtlinearen solver fuer Kennlinien aus
;			der Bibliothek mag_fit.so auf
;
;	--> abgespeckte Version, um Triple-, Quadrupol- oder Penta- Daten 
;	    auszuwerten
;
;
;	x	x-Vektor, hier Spannung
;	y	zugehöriger abhängiger y-Vektor, hier der Sondenstrom
;
;	yfit	Rückgabe-Vektor, enthält summierte Abweichung von vorgegebenen 
;		Stromwerten
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
;	params(8)	V_plasma - V_wall als Rueckgabewert
;	params(9)	unbenutzt
;	params(10)	p_height einer free standing probe (Länge)
;	params(11)	a_width einer flush mounted probe (Breite)
;	params(12)	a_height einer flush mounted probe (Länge)
;	params(13)	a_mass Massezahl des Füllgases
;	params(14)	a_z Kernladungszahl des Füllgases
;	params(15)	B_t Hauptfeld
;	params(16)      Quotinen T_i / T_e
;
;		==> bei Rückkehr aus Fit, stehen in params die neuen
;		    angepaßten Werte, soweit die in fixed (s.u.) verlangt
;		    wurde
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
;	spezielle Hinweise für PIN PROBES:
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
;			weinlich@ipp-garching.mpg.de
;			089/3299-1813
;
;	Stand: 19.06.95
;
;	--------------------------------------------------------------------


FUNCTION call_mag_solver, params, ne_res, te_res, $
			beta=beta, delta_i=delta_i, psi=psi, $
			v1=v1, i1=i1, v2=v2, i2=i2, v3=v3, i3=i3, v4=v4, i4=i4


;@auswertung.common

common		mag_so_aktuell,		was_checked, lib_file, interface


; Initialisierung wenn noetig
check_lib

;
; alle Parameter muessen als double oder long an C uebergeben werden
;

if not keyword_set(v1) then begin
   print, 'call_mag_solver: keine Spannungswerte uebergeben, returning ...'
   return, -1
   endif

n_data = n_elements(v1)
n_points = 1

if not keyword_set(v2) then begin
   print, 'call_mag_solver: nicht genuegend Spannungswerte uebergeben, ' + $
				'returning ...'
   return, -1
   endif

n_points = 2
n_pars   = 2L

if keyword_set(v3) then begin
   n_points = 3
   n_pars   = 3
   endif

if (n_pars gt 3) and (not keyword_set(v4)) then begin
   print, 'call_mag_solver: nicht genuegend Spannungswerte uebergeben, ' + $
				'returning ...'
   return, -1
endif else if keyword_set(v4) then begin
   n_points = 4
   endif

if n_pars gt n_points then begin
   print, 'call_mag_solver: mehr freie Parameter als Stuetzstellen, ' + $
				'returning ...'
   return, -1
   endif 


; erzeuge Datenfelder
data_x = dblarr( n_points, n_data, /nozero)
data_y = dblarr( n_points, n_data, /nozero)

data_x(0,*) = v1
data_y(0,*) = i1

if n_points ge 2 then begin
   data_x(1,*) = v2
   data_y(1,*) = i2
   endif

if n_points ge 3 then begin
   data_x(2,*) = v3
   data_y(2,*) = i3
   endif

if n_points ge 4 then begin
   data_x(3,*) = v4
   data_y(3,*) = i4
   endif

; Datenfelder fuer Rueckgabewerte
ne_res = dblarr(n_data, /nozero )
te_res = dblarr(n_data, /nozero )
beta   = dblarr(n_data, /nozero )
delta_i = dblarr(n_data, /nozero )

; Startwerte fuer Parameter in C
c_params = dblarr(20)
c_params(0)  = params(0)	; n_e
c_params(1)  = params(1)	; T_e
c_params(2)  = params(3)	; alpha1
c_params(3)  = params(4)	; beta
c_params(4)  = params(5)	; d_V
c_params(5)  = params(7)	; alpha2
c_params(6)  = (params(2))	; psi1
c_params(7)  = (params(6))	; psi2


; wenn kein Winkel uebergeben, dann entnehme ihn aus params
if n_elements(psi) lt n_data then dpsi = ne_res*0 + c_params(6) $
			     else dpsi = double(psi)   

; optionale Parameter
c_params(17) = params(10)  ; p_height
c_params(10) = params(11)  ; a_width
c_params(11) = params(12)  ; a_height
c_params(12) = params(13)  ; a_mass
c_params(13) = params(14)  ; a_z
c_params(14) = params(15)  ; B_t


n_pars = long(n_pars)
n_points = long(n_points)
n_data = long(n_data)

; keine Ueberpruefung mathematischer Fehler
junk = check_math(0,1,trap=1)

t1 = systime(1)

c_err = call_external( lib_file, interface, 'mag_solver', $
			data_x, data_y, dpsi, $
			n_points, n_data, $
			c_params, n_pars, $
			ne_res, te_res, beta, $
			delta_i,  $
			value=replicate(0B,20) )

t2 = systime(1)

error = check_math(0,0,trap=1) 

if n_elements(ne_res) gt 100 then begin
   print, ''
   print, 'Zeitaufwand fuer eigentlichen C-Aufruf = ' + strtrim(t2-t1,2) + ' sec'
   print, ''
   endif
   
; alle internen Werte: single-precision
ne_res   = float(temporary(ne_res))
te_res   = float(temporary(te_res))
beta     = float(temporary(beta))
delta_i  = float(temporary(delta_i))

return, c_err

end		; call_mag_solver


