;@init_all
;@call_magfit

; ------------------------------------------------------------------------
; 
; Auswertung von Demo-Kennlinien von Divertorsonden und Mittelebenenmanip.
;
; Referenzwerte:
;
;			   Divertor		   Mittelebene
; mag_doppel		1.69e18 / 7.10		not implemented
;
;
; -------------------------------------------------------------------------

pro main


;@auswertung.common
@physconst.dat


!p.background = !d.n_colors-1
!p.color=0

temp = 'xxxx'
for i=0,5 do print
print, 'welche Testdaten sollen geladen werden?'
print, '1) nichts"attigende LSF-Daten (#1869, 13ua51)'
print, '2) s"attigende LSM-Daten (#5288, 3_5288.dat)'
print
read, 'Datensatz > ', temp

if temp eq '2' then begin

   openr, unit, 'lsm_sample.dat', /get_lun

   vpr = fltarr(500)
   strom = fltarr(500)
   arr2 = [0.,0.]
   i = 0L

   while not eof(unit) do begin
      readf, unit, arr2
      vpr(i) = arr2(0)
      strom(i) = arr2(1)
      i = i+1
      endwhile

   free_lun, unit

   vpr = vpr(0:i-1)
   strom = strom(0:i-1)

   ind = where(vpr lt 85.)
   vpr = vpr(ind)
   strom = strom(ind)
   n = n_elements(vpr)

   l_sonde = 0.002
   b_sonde = 0.0009

   is_lsf = 0 eq 1

endif else begin

   openr, unit, '/afs/ipp/u/mnw/idl/single/demo.dat', /get_lun
   readf, unit, n
   n = fix(n)
   vpr = fltarr(n,/nozero)
   strom = fltarr(n,/nozero)
   readf, unit, vpr
   readf, unit, strom
   close, unit
   free_lun, unit

   l_sonde = 0.040
   b_sonde = 0.005

   is_lsf = 0 eq 0

   endelse

last_num = -1
window, last_num+1, title='Fit mit a_sonde_proj, nichtsaettigende flush probe'
plot, vpr, strom, psy=4


x = double(vpr)
y = double(strom)
yfit = double(strom)

; startwerte fuer Parameter in IDL
idl_params = dblarr(17)
idl_params(0) = 1.e18
idl_params(1) = 10.
idl_params(2) = 89.4 * !dtor
idl_params(3) = 1.
idl_params(4) = 1.
idl_params(5) = 0.
idl_params(6) = 89.6 * !dtor
idl_params(7) = 1.

if is_lsf then begin
   idl_params(10) = -l_sonde
   idl_params(11) = b_sonde
   idl_params(12) = l_sonde
endif else begin
   idl_params(2) = 10. * !dtor
   idl_params(3) = 0.05
   idl_params(10) = -l_sonde
   idl_params(11) = b_sonde
   idl_params(12) = l_sonde / cos(idl_params(2))
   endelse

idl_params(13) = 2.	; Massezahl
idl_params(14) = 1.	; Kernladungszahl
idl_params(15) = 2.	; Magnetfeld, B_t
idl_params(16) = 1.	; T_i / T_e

sp = idl_params


idl_fixed = intarr(n_elements(idl_params))
idl_fixed(*) = 1	; 1 --> kein Fit in entsprechendem Parameter 
idl_fixed(0) = 0	; 0 --> Fit in entsprechendem Parameter
idl_fixed(1) = 0
; kein Fit im Winkel: idl_fixed(2) = 0
idl_fixed(3) = 0
idl_fixed(4) = 0
idl_fixed(5) = 0
idl_fixed(7) = 0


;a_width = b_sonde
;a_height = l_sonde

start = systime(1)
rc = call_magfit( x, y, c_fit, idl_params, idl_fixed)
ende = systime(1)

oplot, vpr, c_fit, line = 2

xyouts, 0.20, 0.8, /normal, 'n_e = '
xyouts, 0.20, 0.77, /normal, 'T_e = '
xyouts, 0.20, 0.74, /normal, 'd_V = '
xyouts, 0.20, 0.71, /normal, 'beta = '

xyouts, 0.32, 0.8, /normal, strtrim(idl_params(0),2)
xyouts, 0.32, 0.77, /normal, strtrim(idl_params(1),2)
xyouts, 0.32, 0.74, /normal, strtrim(idl_params(5),2)
xyouts, 0.32, 0.71, /normal, strtrim(idl_params(4),2)


for i=1, 5 do print
print, 'Zeit in C-Routine '+ strtrim(ende-start) + ' sec.'
for i=1,2 do print
print , 'Parameter:', idl_params
for i=1,2 do print



end		; main



; kurzes Hauptprogramm
!p.background = !d.n_colors-1
!p.color = 0

main

end
