FUNCTION dop_psf_fit, ip_old, dopenv  


  yobs=ip_old 
  n_psf=n_elements(y)
  psfpix=dopenv.psfpix
  psfsig=dopenv.psfsig
  osamp=dopenv.osamp

   names=['AMP0','AMP1','AMP2','AMP3',$
          'AMP4','AMP5','AMP6','AMP7','AMP8','AMP9','AMP10','AMP11','AMP12',$
          'AMP13','AMP14','AMP15','AMP16'] 
          
  par = [0.8,  0.0,0.0,0.001,0.003,0.006, 0.03, 0.11, 0.21, 0.21, 0.11,0.03,0.007,0.002,0.001,0.0, 0.0]
  npar = n_elements(par)

; SETUP THE PARINFO STRUCTURE FOR THE PSF DESCRIPTION
   parinfo = {value: 0.0d,        $    ; double precision  
              fixed: 0,          $
              limited: [0,0],    $ ; use with caution 
              limits: fltarr(2), $
              parname: '?',      $
              step: 0.01d,        $
              relstep: 0.00,     $
              mpside: 2}       ; 0, 1, -1, 2   
   parinfo=replicate(parinfo, npar)
   
   for i=0,npar-1 do begin
      parinfo[i].parname = names[i]
      parinfo[i].value = par[i]
	endfor   

   wts=sqrt(abs(par+0.001))
   n_osamp=15.*osamp
   xarr = (dindgen(2*n_osamp + 1) - n_osamp) / double(osamp) ; 121 elements
   functargs = {psfpix:psfpix, psfsig:psfsig, osamp:osamp, n_psf:n_psf} 

   newpar = mpfitfun('sum_gauss', xarr, yobs, parinfo=parinfo, functargs=functargs,$
                     yfit=ip_new, /quiet)


  return, ip_new
  
end ;pro
