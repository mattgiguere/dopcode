FUNCTION SUM_GAUSS, xarr, par, _EXTRA=functargs 

; called by my_psf to fit sum of Gaussians to the median yrmo SLSF
; fischer Jan 5, 2013

o_pix = n_elements(xarr)           ; 121 subpixels = 30pix*osamp + 1
    extra=functargs
    osamp=extra.osamp
    ppix=extra.psfpix              ; locations of Gaussians offset from central Gaussian
    psig=extra.psfsig              ; width of offset Gaussians 
    tmp_ip=fltarr(o_pix)               ; oversampled PSF 
       
  ; CENTRAL GAUSSIAN
    if ~keyword_set(cntr) then cntr = 0.0
    cen = 0.0 - cntr
    a0 = abs(par[0])
    ;;;fischer jan21,2013
    ;cent_wid = a0*7. > 2.  ;define the central gaussian over this restricted pixel interval
    cent_wid = a0*10. > 2.  ;define the central gaussian over this restricted pixel interval
    lo = cen - cent_wid
    hi = cen + cent_wid
    xx = where(xarr ge lo[0] and xarr le hi[0],nxx)
    if nxx eq 0 then stop
    tmp_ip[xx] = exp(-0.5 * ((xarr[xx] - cen) / a0)^2.) 

  ; CONSTRUCT AND ADD SURROUNDING GAUSSIANS
    for i=1,n_elements(psig)-1 do begin 
      cen=ppix[i]-cntr
      a_i=psig[i]
      if a_i gt 0 then begin        ; gaussian component exists 
         gd_range=5.*a_i            ; side gaussisan spans 5 times the width
         xx=where(xarr ge (cen - gd_range) and xarr le (cen + gd_range))
         tmp_ip[xx]=tmp_ip[xx] + par[i]*exp(-0.5 * ((xarr[xx]-cen) / a_i)^2)
      endif
    endfor 

  ; NORMALIZE 
    dx= (xarr[o_pix-1] - xarr[0]) / (o_pix-1) 
    ip = tmp_ip / (dx * total(tmp_ip))
        
    return, ip
 end
