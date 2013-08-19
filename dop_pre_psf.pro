;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DOP_PRE_PSF
;
; CALLED FROM: 
;   DOP_FIT
;   ip = dop_pre_psf(osamp,par) 
;
; PURPOSE: 
;   This function constructs the PSF as the 
;   with the width of a single Gaussian as a free parameter
;        psfpar[4] is the only par; the width of the central Gaussian for the PSF
;
; PROCEDURE
;   1) using the input psfpar, construct a central Gaussian
;   2) normalize the psf
;
; INPUTS:
;   OSAMP: pass in the pixel oversample
;   PAR: width of central Gaussian, heights of little Gaussians
;
; OUTPUTS:
;   IP: instrumental profile 
;
; initial coding: Debra Fischer, Jan 2008
; based on Valenti, Butler & Marcy, 1995, PASP, 107, 966
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION DOP_PRE_PSF, psfpar, dopenv=dopenv, cntr=cntr, plot=plot

  ; GRAB PSF VARIABLES FROM THE COMMON DOPENV STRUCTURE

  ; SET UP SOME VARIABLES AND ARRAYS 
    osamp=dopenv.osamp
    n_oversamp = 15*osamp 
    xarr = (dindgen(2*n_oversamp + 1) - n_oversamp) / double(osamp) ; 121 elements
    o_pix = n_elements(xarr)           ; 121 subpixels 
    tmp_ip=fltarr(o_pix)               ; fltarr(121) 

  ; CONSTRUCT THE CENTRAL GAUSSIAN
    if ~keyword_set(cntr) then cntr = 0.0
    cen = 0.0 - cntr
    a0 = abs(psfpar)
    cent_wid = a0*4. > 5.
    xx = where(xarr ge cen - cent_wid and xarr le cen + cent_wid)
    tmp_ip[xx] = exp(-0.5 * ((xarr[xx] - cen) / abs(psfpar))^2.) 


;  ; SHIFT TO CENTER - CHECK THIS - MAY WANT TO ELIMINATE!
    fwhm=0.5*max(tmp_ip)
    x2=where(tmp_ip ge fwhm and abs(xarr) lt 12., nx2) ; peak points
    if nx2 ge 3 then begin                            ; shift to center
       dd = where(tmp_ip[x2] eq max(tmp_ip[x2]))  
       ndd=n_elements(dd)
       if ndd le 2  then dd=dd[0]
       if ndd gt 2 then dd=dd[fix(ndd/2.)]
       cntr=xarr[x2[dd]]  
       if abs(cntr) ge 0.1 and abs(cntr) lt 1.2 then $
       tmp_ip[x2] = exp(-0.5 * ((xarr[x2] - cntr) / abs(psfpar))^2.) 
    endif                       ; shift to center

  ; NORMALIZE 
    dx= (xarr[o_pix-1] - xarr[0]) / (o_pix-1) 
    ip = tmp_ip / (dx * total(tmp_ip))

   if keyword_set(plot) then begin
     plot,xarr,ip
  endif
 

  return, ip


end

