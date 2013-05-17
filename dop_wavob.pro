FUNCTION DOP_WAVOB, subch, ncol, poly_ord=poly_ord, demo=demo


   if ~keyword_set(poly_ord) then poly_ord = 3
   wc = dblarr(poly_ord+1)  ; wavelength coefficients
   zp = 0   ; zeroth Doppler pass par's
;   fp = 1   ; first Doppler pass par's
;   lp = 2   ; last Doppler pass par's

   lpix = subch.pixt     ; the iodine wavelength soln
   nind=n_elements(lpix) ; number of chunks
   lwav=dblarr(nind)
   lnwav=dblarr(nind)

   for jj=0, nind-1 do lwav[jj] = double((*subch[jj].free_par)[zp].wcof[0])

 ; FIRST REMOVE THE DOMINANT LINEAR TREND
   lcoef=robust_linefit(lpix,lwav,sav_lin)   ;returns sav_lin, the linear baseline
   nlwav = lwav - sav_lin   ;remove the large linear trend 
 ; FIT THE RESIDUALS 
   resid_lcoef=robust_poly_fit(lpix,nlwav,poly_ord,wmodel_nonlin)

 ;;;  resid_lcoef=poly_fit(lpix,nlwav,poly_ord,/double,yfit=lwavfit)

 ; REFILL THE ZERO-PASS WAVELENGTH COEFFICIENTS FROM THE FIT
  for jj=0, nind-1 do (*subch[jj].free_par)[zp].wcof[0]=double(wmodel_nonlin[jj]+sav_lin[jj])

 ; ONE WAVELENGTH FOR EVERY PIXEL
   xarr=dindgen(ncol)
   wavlinear=poly(xarr,lcoef)
   wavresid=poly(xarr,resid_lcoef)   ; a wavelength for every pixel 
   wav=double(wavlinear+wavresid)

   if keyword_set(demo) then begin 
         !x.margin=[12,5]
         !y.margin=[5,3]
         plot, lpix, nlwav, yra=[-1.5,0.4],/ysty, col=1, xtit='!6 Pixel', $
               ytit='!6 Residual to polynomial fit', ps=4, $
               tit='!6 Smoothed, continuous wavelength soln across each order'
         oplot,lpix,lwavfit,col=222
         oplot,xarr, wavresid, col=155
      endif

return, wav

end

