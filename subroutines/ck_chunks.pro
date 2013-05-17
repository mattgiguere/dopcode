pro ck_chunks, subch_in=subch_in, subch_out=subch_out, disp=disp, wvlen=wvlen

;fischer Jan 27, 2013
; INPUT
; subch_in  ; set of chunk_arr with same ordt

; OUTPUT
; subch_out (if corrections)

flag=0           ; flag for replacement with subch_out
thresh_psf=0.1     ; threshold for changing the central PSF width
thresh_disp=0.002  ; threshold for changing the dispersion 
thresh_wav=0.01    ; threshold for changing the wavelength zero point
ord=subch_in[0].ordt

 ; CONTINUITY OF CENTRAL PSF PARS
   nch_ord = n_elements(subch_in.ordt) 
   ord_psf=fltarr(nch_ord)
   for jj=0,nch_ord-1 do ord_psf[jj]=(*subch_in[jj].free_par)[0].amp[0] 
   xar_tmp=indgen(nch_ord)
   coef=poly_fit(xar_tmp,ord_psf,2)  ; slowly changing
   mod_psf=poly(xar_tmp,coef)        ; model values 

 ; remove any big outliers and fit again
   ord_psf_init=ord_psf              ; save it - you might be curious
   diff=ord_psf-mod_psf
   xbp=where(abs(diff) gt thresh_psf,nxbp)  
   if nxbp gt 0 then begin           ; there be bad points
   	  plot,xar_tmp,ord_psf,title='!6 Order: '+strcompress(string(ord),/rem),ps=8,yra=[0.6,1.2],/ysty
   	  oplot,xar_tmp,mod_psf, linesty=2,thick=1,col=222
      print,'tweaking the psf for order: ',ord
      repeat begin
      	 ibad=where(diff eq max(diff) and abs(diff) gt thresh_psf, nibad)
      	 if nibad gt 0 then begin    ; knock out the worst points first
         	ord_psf[ibad] = mod_psf[ibad]
         	plots,[xar_tmp[ibad],xar_tmp[ibad]],[mod_psf[ibad],mod_psf[ibad]],ps=8,symsize=1.5,col=222
         	(*subch_in[ibad].free_par)[0].amp[0] = mod_psf[ibad] 
         	coef = poly_fit(xar_tmp, ord_psf,2)
         	mod_psf = poly(xar_tmp,coef)
         	diff = ord_psf - mod_psf
         endif
      endrep until nibad eq 0
      ;stop
    endif


if keyword_set(disp) then begin
  ; DISPERSION CONTINUITY 
    disp=fltarr(nch_ord) 		
  ; get the dispersion from the previous pass, across each order
  	for ii=0, nch_ord-1 do disp[ii]=(*subch_in[ii].free_par)[0].wcof[1]
  	xarr=findgen(nch_ord)
  	coef=poly_fit(xarr, disp, 2)
  	model_disp=poly(xarr,coef)  ; smoothed dispersion
  	
  ; remove any big outliers and fit again
    disp_init = disp                   ; save it	
  	diffdisp=disp-model_disp
	sdv_disp=stddev(diffdisp)
  	xbd=where(abs(diffdisp) gt thresh_disp,nxbd)
  	if nxbd gt 0 then begin
    	plot, xarr, disp, /xsty, /ysty, ps=8, col=1, symsize=1.2, xra=[-1,38],$
    		title='!6 Order: '+strcompress(string(ord),/rem),$
     		yra=[min(disp)-0.001, max(disp)+0.001],xtitl='!6Chunk',ytitle='!6 Dispersion'
    	oplot, xarr, model_disp, col=222, thick=1, linesty=2
  		flag=1
  		print,'tweaking the disp for order: ',ord
;  		repeat begin
  			ibadd=where(diffdisp eq max(diffdisp) and abs(diffdisp) gt thresh_disp, nibadd) 
  			if nibadd gt 0 then begin   ;knock out the worst first
  				plots,[xarr[ibadd],xarr[ibadd]],[mod_disp[ibadd],mod_disp[ibadd]],ps=8,symsize=1.5,col=222
;  				ans='n'
;  				read,'is this change in the dispersion OK? (y/n) ',ans
 ; 				if ans eq 'y' then begin
  					disp[ibadd] = mod_disp[ibadd]
  					(*subch_in[ibadd].free_par)[0].wcof[1] = mod_disp[ibadd] 
 ; 				endif
;  				if ans eq 'n' then begin
;  					repeat begin 
;  						print,'click on the desired dispersion' 
;  						cursor,x,ynew
;  						plots,[xarr[ibadd],xarr[ibadd]],[ynew,ynew],ps=8,symsize=1.5,col=fix(randomu(seed)*225)
;						read,'is this change in the dispersion OK? (y/n) ',ans
;						if ans eq 'y' then begin
;							disp[ibadd]=ynew
;							(*subch_in[ibadd].free_par)[0].wcof[1] = ynew
;						endif
;					endrep until ans eq 'y'
;				endif  ; didn't accept initial model dispersion as replacement
			endif ; nibadd 
			coef=poly_fit(xarr, disp, 2)
  			model_disp=poly(xarr,coef)  ; smoothed dispersion
  			diffdisp=disp-model_disp
;		endrep until nibadd eq 0
	endif					
endif ;dispersion

if keyword_set(wvlen) then begin
  ; WAVELENGTH CONTINUITY 
    wav=dblarr(nch_ord) 		
  
  ;get the wavelength soln from the previous pass
 	for ii=0, nch_ord-1 do wav[ii]=double( (*subch_in[ii].free_par)[0].wcof[0] )
  	xarr=dindgen(nch_ord)
	lcoef=poly_fit(xarr,wav,1)
	sav_lin=poly(xarr,lcoef)
	wav_nonlin=wav-sav_lin   ;remove the large linear trend 
	nlcoef=poly_fit(xarr,wav_nonlin,3)
	model_nonlin=poly(xarr,nlcoef)

  ; remove any big outliers and fit again
    wav_nonlin_init = wav_nonlin               ; save it	
  	diffwav=wav_nonlin-model_nonlin
	sdv=stddev(model_nonlin-wav_nonlin)
  	xbw=where(abs(diffwav) gt thresh_wav,nxbw)
  	if nxbw gt 0 then begin
    	plot, xarr, wav_nonlin, /xsty, /ysty, ps=8, col=1, symsize=1.2, xra=[-1,38],$
    		title='!6 Order: '+strcompress(string(ord),/rem),$
     		xtitl='!6Chunk',ytitle='!6 Wavelength'
    	oplot, xarr, model_nonlin, col=222, thick=1, linesty=2
  		flag=1
  		print,'tweaking the wavelength for order: ',ord
;  		ans='n'
;  		repeat begin
  			ibadw=where(diffwav eq max(diffwav) and abs(diffwav) gt thresh_wav, nibadw) 
  			if nibadw gt 0 then begin   ;knock out the worst first
  				plots,[xarr[ibadw],xarr[ibadw]],[model_nonlin[ibadw],model_nonlin[ibadw]],ps=8,symsize=1.5,col=222
;  				read,'is this change in the wavelength OK? (y/n) ',ans
;  				if ans eq 'y' then begin
  					wav_nonlin[ibadw] = model_nonlin[ibadw]
  					(*subch_in[ibadw].free_par)[0].wcof[0] = model_nonlin[ibadw]+sav_lin[ibadw]
;  				endif 
;  				if ans eq 'n' then begin
;  					repeat begin 
;  						print,'click on the desired wavelength' 
;  						cursor,x,yw
;  						plots,[xarr[ibadw],xarr[ibadw]],[yw,yw],ps=8,symsize=1.5,col=fix(randomu(seed)*225)
;						read,'is this change in the wavelength OK? (y/n) ',ans
;						if ans eq 'y' then begin
;							wav_nonlin[ibadw]=yw
;							(*subch_in[ibadw].free_par)[0].wcof[0] = yw + sav_lin[ibadw] 
;						endif
;					endrep until ans eq 'y'
				endif  ; didn't accept initial model wavelength as replacement
			endif ; nibadw 
			nlcoef=poly_fit(xarr,wav_nonlin,3)
			model_nonlin=poly(xarr,nlcoef)
			diffwav=wav_nonlin-model_nonlin
;		endrep until nibadw eq 0
;		wait,1
	endif					
; end ;wav
			
; if no changes then subch_out = subch_in 
	if flag eq 0 then subch_out = subch_in

end ;pro
