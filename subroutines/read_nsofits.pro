pro read_nsofits, file, wnum, spec, head, resol=resol $
                , units=units, silent=silent
;
;Read an FTS spectrum from a National Solar Observatory (NSO) FITS file.
;
;Inputs:
; file (string) name of NSO FITS file to read
; [/silent] (switch) supresses printing of diagnostics
;
;Outputs:
; wnum (vector) uniform wavenumber vector constructed from header keywords
; spec (vector) spectrum read from primary data unit in FITS file
; [head] (string vector) FITS header
; [units=] (integer vector) list of units where diagnostics will be written
; [resol=] (scalar) nominal resolution of spectrum (in 1/cm)

  if n_params() lt 3 then begin
    print, 'syntax: read_nsofits, file, wnum, spec' $
         + ' [,head ,resol= ,units= ,/silent]'
    return
  endif

;Default values for optional keywords.
  if n_elements(units) eq 0 then units = -1	;write only to terminal

;Read data and header from FITS file.
  fits_read, file, spec, head			;astronomy library routine
  n = n_elements(spec)

;Extract wavenumber scale information from FITS header.
  wstart = sxpar(head, 'WSTART', count=count)
  if count eq 0 then begin
    diag, units, 'no WSTART keyword in header of ' + file, /stop
  endif
  wstop = sxpar(head, 'WSTOP', count=count)
  if count eq 0 then begin
    diag, units, 'no WSTOP keyword in header of ' + file, /stop
  endif
  delw = sxpar(head, 'DELW', count=count)
  if count eq 0 then begin
    diag, units, 'no DELW keyword in header of ' + file, /stop
  endif

;Construct wavenumber scale.
  wnum = wstart + delw * dindgen(n)
  if abs(wnum[n-1] - wstop) gt 1d-4 then begin
    diag, units, 'error generating wavenumber scale'
  endif

;Extract spectral resolution from FITS header.
  resol = sxpar(head, 'RESOLUTN', count=count)
  if keyword_set(silent) then return

;Calculate diagnostics.
  wair = 1d8 / [wstart, wstop]
  vactoair, wair

;Print diagnostics.
  fmt = '(f19.5)'
  diag, units, file
  diag, units, 'observed ' + strtrim(sxpar(head, 'DAY'), 2) $
             + ' for ' + strtrim(sxpar(head, 'USER'), 2)
  diag, units, '"' + strtrim(sxpar(head, 'ID'), 2) + '"'
  diag, units, strtrim(n, 2) + ' points'
  diag, units, 'wavenumbers: ' $
             + strtrim(string(wstart, form=fmt), 2) + '-' $
             + strtrim(string(wstop, form=fmt), 2) + ' 1/cm'
  diag, units, 'dispersion: ' $
             + strtrim(string(delw, form=fmt), 2) + ' 1/cm'
  diag, units, 'resolution: ' $
             + strtrim(string(resol, form=fmt), 2) + ' 1/cm'
  diag, units, 'resolving power: ' $
             + strtrim(round(wstart/resol), 2) + '-' $
             + strtrim(round(wstop/resol), 2)
  diag, units, 'air wavelengths: ' $
             + strtrim(string(wair[1], form=fmt), 2) + '-' $
             + strtrim(string(wair[0], form=fmt), 2) + ' Angstroms'

end
