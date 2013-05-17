PRO DOP_FITWAV, chunk_arr, npix=npix, poly_ord=poly_ord, pass=pass, $
                plot=plot, demo=demo, new_chunk_arr

 ; START WITH VD0 AND UPDATE
   thresh=1.5                  ; rejection threshold
   c_light = 2.99792458d8

 ; SET UP: TO FIT A CONTINUOUS WAVELENGTH SOLN ACROSS EACH ORDER, 
 ; STARTING WITH THE I2 WAVELENGTH SOLN AND REPLACE W0 AND WCOF 
   ordr = chunk_arr.ordob              
   gg=uniq(ordr) 
   ordr = ordr(gg)               ; array of orders in the VD
   nord = n_elements(ordr)       ; number of orders

 ; SET UP: ORDER OF POLYNOMIAL FIT
   if ~keyword_set(poly_ord) then poly_ord = 3
   wc = dblarr(poly_ord+1, nord)  ; wavelength coefficients for each order 
   velrms = fltarr(nord) 
   vs = fltarr(nord) 
 
 ; ONE ORDER AT A TIME...
   for i=0,nord - 1 do begin 
     ind = where(chunk_arr.ordob eq ordr[i],nind)  
     if nind eq 0 then stop, 'no orders'
     subchunk_arr = chunk_arr[ind]          ; one order at a time
     fit=fltarr(nind)

     lpix=fltarr(nind)
     lwav=fltarr(nind)
     errl=fltarr(nind)          ; measurement error for poly_fit 
     cpix=fltarr(nind)
     cwav=fltarr(nind)
     errc=fltarr(nind)          ; measurement error for poly_fit 
     lnwav=fltarr(nind)
     cnwav=fltarr(nind)

     err=fltarr(nind*2)
     pix=fltarr(nind*2)
     wav=fltarr(nind*2)

     for j=0,nind-1 do fit[j] = subchunk_arr[j].fit[2]
     medfit=median(fit)
     good = where(fit le thresh*medfit, ngood)
     bad = where(fit gt thresh*medfit, nbad) 

;yuck - I don't like this
     for ii=0,nind-1 do errl[ii] = .125*fit[ii]  ; adhoc scaling for errors from chisq 
     if nbad gt 0 then errl[bad] = 1.5*errl[bad]   ; above, so average chi is about 1.
     errc=2.*errl   ; less weight to these points

     ; IDENTIFY PIXEL AND WAVELENGTHS AT LEFT EDGE AND CENTER OF EACH CHUNK
     lpix = subchunk_arr.pixob
     if pass eq 0 then fitp=2 else fitp=pass-1 
     for jj=0, nind-1 do lwav[jj] = double((*subchunk_arr[jj].free_par)[fitp].wcof[0])

     cpix = subchunk_arr.pixob + npix/2
     for jj=0, nind-1 do cwav[jj] = double((*subchunk_arr[jj].free_par)[fitp].wcof[0]) $
        + (*subchunk_arr[jj].free_par)[fitp].wcof[1]*(npix/2)


     err[0:1]=[errl[0],errc[0]]
     pix[0:1]=[lpix[0],cpix[0]]
     wav[0:1]=[lwav[0],cwav[0]]

     for l=1,nind - 1 do begin 
        m=l*2 
        err[m] = errl[l]
        err[m+1] = errc[l]
        pix[m] = lpix[l]
        pix[m+1] = cpix[l]
        wav[m] = lwav[l]
        wav[m+1] = cwav[l]
     endfor


     ; REMOVE DOMINANT LINEAR WAVELENGTH SOLN TO (NUMERICALLY)
     ; ASSIST POLY_FIT (ADD BACK IN LATER)
     lcoef=poly_fit(lpix,lwav,1, /double,yfit=lin_lpix)         ; same slope for left and center
     nlwav = lwav - lin_lpix

     coef=poly_fit(pix,wav,1,/double,yfit=lin_pix)
     nwav = wav - lin_pix

     ; FIT POLYNOMIAL TO RESIDUALS ACROSS THE ORDER
     ; LWAVFIT AND CWAVFIT ARE THE LEFT AND CENTRAL WAVELENGTH VALUES
     ; FOR EACH CHUNK, FOUND WITH THE NEW POLYNOMIAL WAVELENGTH SOLUTION

     ; COEFFICIENTS, FITTING LEFT EDGE PIXELS
     resid_lcoef=poly_fit(lpix,nlwav,poly_ord,/double,$
          measure_errors=err, yfit=lwavfit,yerror=lwv_rms, chisq=chil)
     resid_lcoef=reform(resid_lcoef)

     resid_coef=poly_fit(pix,nwav,poly_ord,/double,$
          measure_errors=err, yfit=wavfit,yerror=wv_rms, chisq=chi)
     resid_coef=reform(resid_coef)


     for l=0,nind-1 do begin
        m=l*2
        lnwav[l] = wavfit[m] + lin_pix[m]
        cnwav[l] = wavfit[m+1] + lin_pix[m+1]
     endfor

   ; update the wavelength par's wiht smoothed soln
     for jj = 0, nind-1 do (*subchunk_arr[jj].free_par)[pass].wcof[0] = lwavfit[jj] + lin_pix[jj]
     for jj = 0, nind-1 do (*subchunk_arr[jj].free_par)[pass].wcof[1] = (cnwav[jj] - lnwav[jj]) / (npix/2)

     if keyword_set(demo) then begin 
        if i eq 0 then begin 
           !x.margin=[12,5]
           !y.margin=[5,3]
           plot, lpix, lnwav,col=1, xtit='!6 Pixel', $
                 ytit='!6 Residual to polynomial fit', ps=4, $
                 tit='!6 Smoothed, continuous wavelength soln across each order'
           oplot,lpix,lwavfit,col=222
        endif
     endif
     chunk_arr[ind]=subchunk_arr          ; one order at a time

endfor

end

