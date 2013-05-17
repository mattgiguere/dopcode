pro fitpsf, gauss=gauss, bspline=bspline, plot=plot, allgauss=allgauss

  restore,'/tous/mir7/files_df/vdadgAVG_1207'   ;for example...
  indx=[5,10,23,66,126,458,500]  ;chunks to fit
  n_indx=n_elements(indx)

if keyword_set(gauss) then begin

; L-M fit the medpsf for each chunk with a sum of gaussians 
; store the amplitudes in a generic yrmo cdbfile 
; CRITICAL: must use the same psfpix and psfsig in ctio4k_init 
  osamp = 4
  psfpix = [0.0, -4.0, -3.6, -3.0, -2.4, -1.8, -1.2, -0.8, -0.4, 0.4, 0.8,  1.2,  1.8,  2.4, 3.0, 3.6, 4.0]
  psfsig = [1.05, 0.0,  0.8,  0.7,  0.7,  0.7,  0.7,  0.6,  0.5, 0.5, 0.6,  0.7,  0.7,  0.7, 0.7, 0.8, 0.0]
  par = [0.8, 0.0, 0.0,0.001,0.003,0.006, 0.03, 0.11, 0.21, 0.21, 0.11,0.03,0.007,0.002,0.001,0.0, 0.0]

  npar=n_elements(par) 
  n_psf=121

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
   par1info=parinfo  ; just the central gaussian width
   
   names=['AMP0','AMP1','AMP2','AMP3',$
          'AMP4','AMP5','AMP6','AMP7','AMP8','AMP9','AMP10','AMP11','AMP12',$
          'AMP13','AMP14','AMP15','AMP16'] 
          
   for i=0,npar-1 do begin
      parinfo[i].parname = names[i]
          parinfo[i].value = par[i]
       endfor   
;   parinfo[0].limited=[1,1] 
;   parinfo[0].limits=[0.68,1.2]

   n_osamp=15.*osamp
   xarr = (dindgen(2*n_osamp + 1) - n_osamp) / double(osamp) ; 121 elements
   functargs = {psfpix:psfpix, psfsig:psfsig, osamp:osamp, n_psf:n_psf} 
   title1='!6 Composite Gaussians for Obs (black) and Model (red-dashed)'

for i=0,n_indx-1 do begin 
   ip=vdav[indx[i]].psf
   yobs=ip
   newpar = mpfitfun('sum_gauss', xarr, yobs, parinfo=parinfo, functargs=functargs,$
                     yfit=ip_model,/quiet)

   if keyword_set(plot) and keyword_set(gauss) then begin ; overplot the SLSF and the fitted model
      ny=sum_gauss(xarr,newpar,_extra=functargs) 

      if keyword_set(hdcopy) then ps_open,'psf_components',/color
      plot,4.*(xarr+15),yobs,col=0,/ysty,yra=[-0.1,0.7],xra=[40,80],thick=2  
      oplot,4.*(xarr+15),ny,col=222,linesty=2, thick=2
      for j=1,npar-1 do plots,[(60+psfpix[j]*4), (60+psfpix[j]*4)],[0.0, 0.6],linesty=2,thick=0.5,col=155
      plots,[60,60],[0.0, 0.6],linesty=2,thick=0.5,col=155

      if keyword_set(allgauss) then begin
         tmp_ip=fltarr(n_psf)   ; oversampled PSF 
         cen=0
         a0 = abs(newpar[0])
         cent_wid = a0*5. > 2.  ;define the central gaussian over this restricted pixel interval
         lo = cen - cent_wid
         hi = cen + cent_wid
         xx = where(xarr ge lo[0] and xarr le hi[0],nxx)
         tmp_ip[xx] = exp(-0.5 * ((xarr[xx] - cen) / a0)^2.) 
         oplot,xx,tmp_ip[xx]*max(ny),col=90,linesty=4,thick=4
         xyouts,42,0.14,/data,'!6 Arrows point to psfpix positions',size=1.8
         for ii=1,npar-1 do begin ;ii is the psfpar number
            cen=psfpix[ii]
            a_i=psfsig[ii]
            if a_i gt 0 then begin    ; gaussian component exists 
               gd_range=5.*a_i        ; side gaussisan spans 5 times the width
               xx=where(xarr ge (cen - gd_range) and xarr le (cen + gd_range))
               tmp_ip[xx]= newpar[ii]*exp(-0.5 * ((xarr[xx]-cen) / a_i)^2)
               oplot,xx,tmp_ip[xx]*max(ny),col=ii*15,linesty=4,thick=4
            endif
            arrow, (60+(psfpix[ii]*4)), 0.1, (60+(psfpix[ii]*4)), 0.,/data,thick=2
         endfor                 ;ii
      endif ;allgauss
      if keyword_set(hdcopy) then ps_close
   endif   ; plot gauss
   stop
endfor
endif   ;gauss


if keyword_set(bspline) then begin
   psfbsplnplaces=[0,30,45,50,53,56,59,61,64,67,70,75,90,120] ;coefficients for the bspline model
   npar=n_elements(psfbsplnplaces)

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
;   par1info=parinfo  ; just the central gaussian width

;   psfbsinvvar=dblarr(n_elements(psfpix)) ;inverse variance for bspline weighting
   psfbsord=4  ; the order to use for the bspline
   
   for i=0,npar-1 do begin
       parinfo[i].value = psfbsplnplaces[i]
   endfor   

for i=0,n_indx-1 do begin 
   osamp=4
   n_osamp=15.*osamp
   xarr = (dindgen(2*n_osamp + 1) - n_osamp) / double(osamp) ; 121 elements
   ip=vdav[indx[i]].psf
   yobs=ip
   newpar = mpfitfun('my_bspline', xarr, yobs, parinfo=parinfo, functargs=functargs,$
                     yfit=ip_model,/quiet)
   ny=my_bspline(xarr,newpar) 
   plot,4.*(xarr+15),ip
   oplot,4.*(xarr+15),ny
stop
endfor  ;indx

endif  ; bspline

end   ;pro
