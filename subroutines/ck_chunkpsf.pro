pro ck_chunkpsf, subch_in=subch_in, subch_out=subch_out

;fischer Jan 27, 2013
; INPUT
; subch_in  ; set of chunk_arr with same ordt

; OUTPUT
; subch_out (if corrections)

flag=0           ; flag for replacement with subch_out
thresh_psf=0.1     ; threshold for changing the central PSF width
ord=subch_in[0].ordt

 ; CONTINUITY OF CENTRAL PSF PARS
   nch_ord = n_elements(subch_in.ordt) 
   ord_psf=fltarr(nch_ord)  ;central gaussian
   for jj=0,nch_ord-1 do ord_psf[jj]=(*subch_in[jj].free_par)[0].amp[0] 
   xar_tmp=indgen(nch_ord)
   coef=poly_fit(xar_tmp,ord_psf,2)  ; slowly changing
   mod_psf=poly(xar_tmp,coef)        ; model values 

 ; remove any big outliers and fit again
   diff=ord_psf-mod_psf
   xbp=where(abs(diff) gt thresh_psf,nxbp) 
   xgp=where(abs(diff) lt thresh_psf,nxgp)
    
   if nxbp gt 0 then begin           ; there be bad points
   	  plot,xar_tmp,ord_psf,title='!6 Order: '+strcompress(string(ord),/rem),ps=8,yra=[0.6,1.2],/ysty
	  newcoef=poly_fit(xar_tmp[xgp],ord_psf[xgp],2)
	  nmod_psf=poly(xar_tmp,newcoef)
   	  oplot,xar_tmp,nmod_psf, linesty=2,thick=1,col=222
	  ord_psf[xbp]=nmod_psf[xbp]
	  oplot,xar_tmp[xbp],ord_psf[xbp],ps=8,symsize=1.2, col=155 
	  flag=1
	  for ii =0,nxbp-1 do (*subch_out[xbp[ii]].free_par)[0].amp[0]=ord_psf[xbp[ii]]
    endif



; if no changes then subch_out = subch_in 
	if flag eq 0 then subch_out = subch_in

end ;pro
