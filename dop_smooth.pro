FUNCTION DOP_SMOOTH, order, sub_chunk, PASS=PASS, DISP=DISP, WAV=WAV, PLOT=PLOT

;fischer july 2011
;fischer apr 2013
; improved the smoothing algorithm for dispersion and wavelength continuity

if keyword_set(disp) and keyword_set(wav) then stop,'choose either disp or wav' 

  ; DISPERSION 
if keyword_set(disp) then begin 
  nx=n_elements(sub_chunk.ordt)
  disp=fltarr(nx)
 		
  ; get the dispersion across the order
  	for ii=0, nx-1 do disp[ii]=(*sub_chunk[ii].free_par)[pass].wcof[1]
  	xarr=findgen(nx)
  	coef=svdfit(xarr, disp, 2, yfit=model_disp)
	sdv_disp=stddev(model_disp - disp)>0.00001

  ; fix the > sigthresh sigma outliers 
	sigthresh=4.
	diffd=abs(model_disp-disp)
	xbadd=where(diffd gt sigthresh*sdv_disp,nxbadd)
	if nxbadd gt 0 then begin 
		print,'adjusting dispersion'
     	plot, xarr, model_disp-disp, /xsty, /ysty, ps=8, col=1, symsize=1.2, xra=[-1,38],$
     			yra=[min(model_disp-disp)-0.002, max(model_disp-disp)+0.002],xtitl='!6Chunk',ytitle='!6 Dispersion',$
     			title='!6 Order: '+strcompress(order,/rem)
     	oplot, xarr, fltarr(nx), col=222, thick=3, linesty=2
		if nxbadd gt 0 then oplot, xarr[xbadd], model_disp[xbadd]-disp[xbadd], ps=8,col=222
			     		
      ; now iterate, fitting only points with stddev less than sigthresh * stddev
		xgd=where(diffd lt sigthresh*sdv_disp,nxgd) 
		newcoef=svdfit(xarr[xgd],disp[xgd],2,yfit=new_model_disp)
		diffd2 = abs(new_model_disp - disp)
		sdv_disp2 = stddev(new_model_disp - disp)
		xbad2=where(diffd2 gt sigthresh*sdv_disp2,nxbad2)
			
	  ; replace any bad values of dispersion
		if nxbad2 gt 0 then begin
			disp[xbad2] = new_model_disp[xbad2]
  			!p.charsize=2
			oplot,xarr[xbad2],model_disp[xbad2]-disp[xbad2], ps=8,col=155,symsize=1.2
		endif ; nxbad2 gt 0
	endif  ;nxbadd gt 0

	smdisp = disp
	return, smdisp
endif ;disp 

  ; WAVELENGTH 
if keyword_set(wav) then begin 
  nx=n_elements(sub_chunk.ordt)
  wav=dblarr(nx)
  polyord=4
  	;get the wavelength soln for an order
 		for ii=0, nx-1 do wav[ii]=double( (*sub_chunk[ii].free_par)[pass].wcof[0] )
  		xarr=dindgen(nx)
		nlcoef=svdfit(xarr, wav, polyord, yfit=wmodel_nonlin)
		diff=abs(wmodel_nonlin-wav)
		sdv=stddev(diff)
	
		; fix the > 3 sigma outliers 
		sigthresh=5.
		xbadw=where(diff gt sigthresh*sdv,nxbadw) 
		if nxbadw gt 0 then begin
			print,'adjusting wavelength zero pt'
			plot,xarr,wmodel_nonlin-wav,/xsty,/ysty, ps=8, col=1, symsize=1.2, $
				xtitl='!6Chunk',ytitle='!6Wavelength ', $
     			title='!6 Order: '+strcompress(order,/rem),$
     			yra=[min(wmodel_nonlin-wav)-0.05, max(wmodel_nonlin-wav)+0.05]
     		oplot,xarr,fltarr(nx),col=222,linesty=2,thick=2  ;model
			
     		; now iterate, fitting only points with stddev less than sigthresh * stddev			
			xgdw=where(diff lt sigthresh*sdv,nxgdw)
			newwnlcoef = svdfit(xarr[xgdw],wav[xgdw],polyord,yfit=nwmodel_nonlin)
			tmp=poly(xarr,newwnlcoef)	
			diff2 = abs(tmp - wav)
			sdv_diff2 = stddev(diff2) 
			xbadw2 = where(diff2 gt sigthresh*sdv_diff2,nxbadw2) 
			
			; replace any significant outliers for the wavelength
			if nxbadw2 gt 0 then begin
				oplot,xarr[xbadw2],tmp[xbadw2]-wav[xbadw2],ps=8,symsize=1.2,col=222  ;outlier
				wav[xbadw2]=tmp[xbadw2]
				!p.charsize=2
			oplot, xarr[xbadw2], tmp[xbadw2]-wav[xbadw2], ps=8, col=155, symsize=1.2 ;replacement point
			endif
		endif 	; nxbadw gt 0 

		wavsm = wav
		return, wavsm
endif ;wav 	

end
