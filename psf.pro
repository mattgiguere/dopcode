;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PSF
;
; PURPOSE: 
;   This function constructs the PSF as the 
;   with the width of a single Gaussian as a free parameter
;
; PROCEDURE
;   1) using the input psfpar, construct a central Gaussian
;   2) normalize the psf
;
; OUTPUTS:
;   IP: instrumental profile 
;
; initial coding: Debra Fischer, Jan 2008
; based on Valenti, Butler & Marcy, 1995, PASP, 107, 966
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION PSF, par

  psfpar=par

  ; GRAB PSF VARIABLES FROM THE COMMON DOPENV STRUCTURE

  ; SET UP SOME VARIABLES AND ARRAYS 
    osamp=4
    n_oversamp = 15*osamp 
    xarr = (dindgen(2*n_oversamp + 1) - n_oversamp) / double(osamp) ; 121 elements
    o_pix = n_elements(xarr)           ; 121 subpixels 
    tmp_ip=fltarr(o_pix)               ; fltarr(121) 

  ; CONSTRUCT THE CENTRAL GAUSSIAN
    if ~keyword_set(cntr) then cntr = 0.0
    cen = 0.0 - cntr
    a0 = abs(psfpar)
    cent_wid = a0*4. > 2.
    xx = where(xarr ge cen - cent_wid and xarr le cen + cent_wid)
    tmp_ip[xx] = exp(-0.5 * ((xarr[xx] - cen) / abs(psfpar))^2.) 


  ; SHIFT TO CENTER - CHECK THIS - MAY WANT TO ELIMINATE!
    fwhm=0.5*max(tmp_ip)
    x2=where(tmp_ip ge fwhm and abs(xarr) lt 5., nx2) ; peak points
    if nx2 ge 3 then begin                            ; shift to center
       dd = where(tmp_ip[x2] eq max(tmp_ip[x2]))  &  dd=dd[0]
       cntr=xarr[x2[dd]]  
       if abs(cntr) ge 0.1 and abs(cntr) lt 1.2 then $
       tmp_ip[x2] = exp(-0.5 * ((xarr[x2] - cntr) / abs(psfpar))^2.) 
    endif                       ; shift to center

  ; NORMALIZE 
    dx= (xarr[o_pix-1] - xarr[0]) / (o_pix-1) 
    ip = tmp_ip / (dx * total(tmp_ip))

  return, ip


end

