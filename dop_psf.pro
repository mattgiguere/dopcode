;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DOP_PSF
;
; CALLED FROM: 
;   DOP_FIT
;   ip = dop_psf(osamp,par) 
;
; PURPOSE: 
;   This function constructs the PSF as the 
;   sum of several Gaussians
;        psfpar[0] is the width of the central Gaussian for the PSF
;        psfpar[1] - psfpix[16] are amplitudes of little Gaussians for the PSF
;
;   Other free paramameters are not passed down to this program, but
;   solved in dop_fit.pro 
;        w0, the wavelength of the first pixel in the chunk
;        the linear dispersion across a chunk
;        z, the Doppler shift, v/c or dlam/lam
;        normalization factor
;
; PROCEDURE
;   1) using the input psfpar's, construct a central Gaussian
;   2) using the input psfpar's, construct flanking Gaussians
;   3) normalize the psf and recenter 
;
; INPUTS:
;   OSAMP: pass in the pixel oversample
;   PAR: width of central Gaussian, heights of little Gaussians
;
; OUTPUTS:
;   LSF: Line spread function (or instrumental profile)
;
; initial coding: Debra Fischer, Jan 2008
; based on Valenti, Butler & Marcy, 1995, PASP, 107, 966
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION DOP_PSF, psfpar, dopenv=dopenv, xarr=xarr, cntr=cntr, plot=plot

  ; GRAB PSF VARIABLES FROM THE COMMON DOPENV STRUCTURE
    osamp=dopenv.osamp
    ppix=dopenv.psfpix 
    psig=dopenv.psfsig 
    npsfpix=dopenv.npsfpix
    
  ; SET UP SOME VARIABLES AND ARRAYS 
    n_oversamp = 15.*osamp             ; n_elements(ppix) * osamp        
    if ~keyword_set(xarr) then $
       xarr = (dindgen(2*n_oversamp + 1) - n_oversamp) / double(osamp) ; 121 elements

    ip=dop_gauss(psfpar, dopenv, xarr=xarr, cntr=cntr)

  shft_style = dopenv.shft_style 

  if shft_style eq 1 then begin
  ; SHIFT TO CENTER - CHECK THIS - MAY WANT TO ELIMINATE!
    fwhm=0.3*max(ip)  ;pad it a bit so 0.4 instead of 0.5
    x2=where(ip ge fwhm and abs(xarr) lt 6., nx2) ; peak points +/- from center
	if nx2 lt 3 then print,'not enough pixels to center the psf'
	del=60.-x2[0]
    if nx2 ge 3 then begin                            ; shift to center
       dd = where(ip[x2] eq max(ip[x2]))  &  dd=dd[0]
       ndd=n_elements(dd)
       if ndd le 2  then dd=dd[0]
       if ndd gt 2 then dd=dd[fix(ndd/2.)]
       cntr=xarr[x2[dd]]  
       if abs(cntr) ge 0.1 and abs(cntr) lt 1.2 then $
          ip = dop_gauss(psfpar, dopenv, xarr=xarr, cntr=cntr)
    endif                       ; shift to center
 endif

 if shft_style eq 2 then begin 
    dx= (xarr[npsfpix-1] - xarr[0]) / (npsfpix-1) 
    c1 = dx * total(xarr*ip)
    c2 = dx * total(ip)
    com = c1 / c2 
    cntr = -com/dx
    if abs(cntr) gt 0.25 then begin  
        ;print, 'cntr: ',cntr
        ip = shift_interp(ip, cntr)
    endif
 endif ; shft_style=2

 
   if keyword_set(plot) then begin
     plot,xarr,ip
  endif
 
  return, ip


end
