Function dop_transpar,par,toggle, npix

;Transform from original free parameters (par) to new parameters (newpar).
;The new parameters are linear combinations of the original.
;The NEWPAR are designed to enhance CHI-SQ minimization ---
; they lie parallel or perpendicular to the Chi-Sq valleys.
; from Marcy & Butler, modified by Fischer nov09 for dop_starsolve test

;INPUT:
;      PAR     fltarr(15)    input parameters
;      TOGGLE  1 or -1       1 (-1) for forward (backward) transformation

;OUTPUT: 
;     The function returns the new parameters, NEWPAR
;     The new DSTEP array is also passed (toggle=1), via the common block
;     
;CAUTION:  DSTEP array is automatically transformed and passed in common
;
;Aug-95 GWM

  c_light = 2.99792458d8
  newpar = par                  ;Initialize newpar

  dvel = 20.                    ;20 m/s  = standard velocty unit
  ; dvel = 10.                  ;10 m/s  = standard velocty unit ; dave 4/28/06
  dstep[2] = dvel/c_light       ;20 m/s achieves lowest Chi-sq

  IF toggle ne 1 and toggle ne -1 then begin ;Verify toggle = 1 or -1
     print,'TRANSPAR:  Error in toggle (+1 or -1)'
  END

  IF toggle eq 1 then begin     ;Forward Transformation
     newpar[0] = par[0] + (npix/2.)*par[1]  ;wav at pixel = 40
     newpar[1] = par[0] - (npix/2.)*par[1]  ;wav at pixel = -40 (!)
     dstep[0] = fix(par[0]) * dvel/c_light  ;new dstep (~ 10 m/s)
     dstep[1] = fix(par[0]) * dvel/c_light  ;new dstep (~ 10 m/s)
  ENDIF

  IF toggle eq -1 then begin                 ;Transform Back
  ; De-Transform pars 11 and 13 (par[0] and par[1])
     newpar[0] = 0.5*(par[0] + par[1])
     newpar[1] = (par[0] - par[0])/(1.0*npix)
  ENDIF

  return, newpar
end
