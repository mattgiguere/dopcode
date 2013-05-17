pro dop_star_psf, yrmo=yrmo, tag=tag, mode=mode, smooth=smooth, plot=plot, hdcopy=hdcopy

; PURPOSE
;  	calculate a median PSF and disp from the Bstar observations and use this 
;	as input for the stellar Dop models
;
; INPUT
; 	YRMO: year and month for extracting the median Bstar psf
; 	TAG: tag used in the standard Dop code analysis of Bstars 
;
; OPTIONAL
;  	PLOT: a 9-chunk montage of the individ SLSF models with the median overplotted
;	HDCOPY
;
; fischer 2013 jan

if ~keyword_set(tag) then tag='adg'
if ~keyword_set(yrmo) then yrmo='1207' 
if ~keyword_set(mode) then mode='narrow_slit'
year='20'+strmid(yrmo,0,2)
cdbfile_path='/tous/mir7/files_df/'

ff=file_search('/tous/mir7/logstructs/'+year+'/'+yrmo+'*log.dat',count=count)
cdb=' '  ;seed the string array

;collect all Bstar cdbfiles for the yrmo time interval
for j=0,count-1 do begin
	restore,ff[j]
	x=where(strcompress(log.decker, /rem) eq mode and $
        strcompress(log.iodcell, /rem) eq 'IN',nx)
	x1=strpos(ff[j], yrmo)
	date=strmid(ff[j],x1,6) 
	print,'DATE: ',date, nx
	; build array of cdb files
	for i=0,nx-1 do begin
   		obnm='a'+strcompress(log[x[i]].prefix,/rem)+'.'+strcompress(log[x[i]].seqnum,/rem)
   		;print,obnm
   		cdtmp='cd'+tag+'b'+strcompress(log[x[i]].object,/rem)+'_'+obnm
   		filecheck=file_search('/tous/mir7/files_df/'+cdtmp,count=fcount)
   		print,cdtmp, fcount
   		if fcount eq 1 then cdb=[cdb,cdtmp] 
   	endfor
endfor
n_cdb=n_elements(cdb)
psfstack=cdb[1:n_cdb-1]  ;trim the bogus seed

restore,cdbfile_path+psfstack[0]                 ;chunk_arr
ds=size(chunk_arr[0].ip)
n_chunks=n_elements(chunk_arr.ordt) ;760 chunks 
n_psf=ds[1]                         ;121 element array
nstack=n_elements(psfstack)         ;number of Bstars in the yrmo interval


; calculate the median and stddev SLSF and dispersion for each chunk 
; optionally, make a montage /plot for 9 chunks
	xarrpsf=indgen(n_psf)
	medpsf=fltarr(n_psf,n_chunks)
	errpsf=fltarr(n_psf,n_chunks)
	meddisp=fltarr(n_chunks)
	disp_rms=fltarr(n_chunks)
	medwav=fltarr(n_chunks)
	wav_rms=fltarr(n_chunks)
	amp0=fltarr(n_chunks) 
	limit_amp0=fltarr(n_chunks) 
	psf2darr=fltarr(n_psf, n_chunks, nstack) 
	disp2darr=fltarr(n_chunks, nstack) 
	wav2darr=fltarr(n_chunks, nstack)
	amp0_arr=fltarr(n_chunks, nstack)
	if keyword_set(hdcopy) then ps_open,'psf_arr',/color
	for k=0,nstack-1 do begin 
		restore, cdbfile_path+psfstack[k]
		for l=0,n_chunks-1 do psf2darr[*,l,k]=chunk_arr[l].ip[*,2]
		for l=0,n_chunks-1 do amp0_arr[l,k]=(*chunk_arr[l].free_par)[1].amp[0]
		for l=0,n_chunks-1 do disp2darr[l,k]=(*chunk_arr[l].free_par)[2].wcof[1]
		for l=0,n_chunks-1 do wav2darr[l,k]=(*chunk_arr[l].free_par)[2].wcof[0]
	endfor
	for l=0,n_chunks-1 do meddisp[l]=median(disp2darr[l,*])
	for l=0,n_chunks-1 do medwav[l]=median(wav2darr[l,*])
	for l=0,n_chunks-1 do disp_rms[l]=stddev(disp2darr[l,*])
	for l=0,n_chunks-1 do wav_rms[l]=stddev(wav2darr[l,*])
	for l=0,n_chunks-1 do amp0[l]=median(amp0_arr[l,*])
	for l=0,n_chunks-1 do limit_amp0[l]=stddev(amp0_arr[l,*])
	for m=0,n_psf-1 do begin  
		for l=0,n_chunks-1 do medpsf[m,l]=median(psf2darr[m,l,*])
		for l=0,n_chunks-1 do errpsf[m,l]=stddev(psf2darr[m,l,*])
	endfor ;m
	!p.multi=[0,3,3]
	ckarr=[81, 92, 109, 309, 320, 337, 537, 548, 565]   
	for l=0,8 do begin
		plot,psf2darr[*,ckarr[l],0],/xsty,/ysty, yra=[0,0.7],xra=[40,80]
		for n=0,nstack-1 do oplot,psf2darr[*,ckarr[l],n],col=fix(200/n)
		oplot,medpsf[*,ckarr[l]],col=222,linesty=2
		!p.color=222
		oploterr,xarrpsf,medpsf[*,ckarr[l]],errpsf[*,ckarr[l]], 3
		!p.color=0
		plots,[60, 60], [0,0.7]
	endfor
	if keyword_set(hdcopy) then ps_close
	!p.multi=[0,1,1]

; medpsf[*,chunk] gives the median psf for a chunk
; err[*,chunk] gives the stddev uncertainties for each oversampled psf point 
	savfile='SLSF_'+yrmo+'.dat'	
	save, medpsf,errpsf,f=savfile

; now L-M fit the medpsf for each chunk with a sum of gaussians 
; store the amplitudes in a generic yrmo cdbfile 
; CRITICAL: must use the same psfpix and psfsig in ctio4k_init 
  osamp = 4
  psfpix = [0.0, -4.0, -3.6, -3.0, -2.4, -1.8, -1.2, -0.8, -0.4, 0.4, 0.8,  1.2,  1.8,  2.4, 3.0, 3.6, 4.0]
  psfsig = [1.05, 0.0,  0.8,  0.7,  0.7,  0.7,  0.7,  0.6,  0.5, 0.5, 0.6,  0.7,  0.7,  0.7, 0.7, 0.8, 0.0]
  par = [0.8,  0.0,0.001,0.003,0.006, 0.03, 0.11, 0.21, 0.21, 0.11,0.03,0.007,0.002,0.001,0.0, 0.0, 0.0]

  npar=n_elements(par) 

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
   parinfo[0].limited=[1,1] 
   parinfo[0].limits=[0.68,1.2]

	n_osamp=15.*osamp
    xarr = (dindgen(2*n_osamp + 1) - n_osamp) / double(osamp) ; 121 elements
	functargs = {psfpix:psfpix, psfsig:psfsig, osamp:osamp, n_psf:n_psf} 
	
; Levenberg-Marquardt fit with mpfitfun	
; store in a template cdbfile
cdbf_mid=psfstack(fix(nstack)/2)
restore,cdbfile_path+cdbf_mid  ;chunk_arr template 
xh=strpos(cdbf_mid,'argb') 
vdnm_mid='vd'+tag+strmid(cdbf_mid,xh+4,strlen(cdbf_mid)-xh)

; RESTORE A TEMPLATE VD
restore,cdbfile_path+vdnm_mid 

for i =0,n_chunks-1 do begin   ;i is chunk number
	err_y=errpsf[*,i] ;stddev of median PSF for a chunk
	yobs=medpsf[*,i]  ;median PSF for a chunk
	; fitted Gaussian amplitudes 
	parinfo[0].value=amp0[i]
	lim=limit_amp0[0]>0.01
	parinfo[0].limits=[amp0[i]-lim, amp0[i]+lim]
	;FUNCTION DOP_PRE_PSF, psfpar, dopenv=dopenv, cntr=cntr
	newpar = mpfitfun('sum_gauss', xarr, yobs, parinfo=parinfo, functargs=functargs,$
					err=err_y, yfit=ip,/quiet)
	print,newpar[0]				
	; store median SLSF and disp in template files
	vd[i].psf=medpsf[*,i]  
	vd[i].par=newpar
	vd[i].w0=fix(medwav[i])
	vd[i].wcof[0]=medwav[i]-vd[i].w0
	vd[i].wcof[1]=meddisp[i]
	vd[i].vel=disp_rms[i]   ;IMPORTANT: this holds the limits in chunk disp for parinfo
	chunk_arr[i].ip[*,1]=medpsf[*,i] 
	chunk_arr[i].ip[*,2]=medpsf[*,i]
	for j=0,npar-1 do (*chunk_arr[i].free_par)[1].amp[j]=newpar[j]
	(*chunk_arr[i].free_par)[1].wcof[1]=meddisp[i] ;median dispersion
	(*chunk_arr[i].free_par)[2].wcof[1]=meddisp[i] ;median dispersion
	(*chunk_arr[i].free_par)[1].wcof[0]=medwav[i] ;median wav
	(*chunk_arr[i].free_par)[2].wcof[0]=medwav[i] ;median wav
		

    if keyword_set(plot) then begin ; overplot the median SLSF and the fitted model
        ny=sum_gauss(xarr,newpar,_extra=functargs) 
    	if keyword_set(hdcopy) then ps_open,'psf_components',/color
    	title1='!6 Composite Gaussians for Obs (black) and Model (red-dashed)'
    	plot,4.*(xarr+15),yobs,col=0,/ysty,yra=[-0.1,0.7],xra=[40,80],thick=2  
    	oplot,4.*(xarr+15),ny,col=222,linesty=2, thick=2
        tmp_ip=fltarr(n_psf)               ; oversampled PSF 
  		cen=0
        a0 = abs(newpar[0])
    	cent_wid = a0*5. > 2.  ;define the central gaussian over this restricted pixel interval
    	lo = cen - cent_wid
    	hi = cen + cent_wid
    	xx = where(xarr ge lo[0] and xarr le hi[0],nxx)
    	tmp_ip[xx] = exp(-0.5 * ((xarr[xx] - cen) / a0)^2.) 
    	oplot,xx,tmp_ip[xx]*max(ny),col=90,linesty=4,thick=4
    	xyouts,42,0.14,/data,'!6 Arrows point to psfpix positions',size=1.8
    	for ii=1,npar-1 do begin  ;ii is the psfpar number
    	    cen=psfpix[ii]
      		a_i=psfsig[ii]
      		if a_i gt 0 then begin        ; gaussian component exists 
         		gd_range=5.*a_i            ; side gaussisan spans 5 times the width
         		xx=where(xarr ge (cen - gd_range) and xarr le (cen + gd_range))
         		tmp_ip[xx]= newpar[ii]*exp(-0.5 * ((xarr[xx]-cen) / a_i)^2)
         		oplot,xx,tmp_ip[xx]*max(ny),col=ii*15,linesty=4,thick=4
         	endif
         	arrow, (60+(psfpix[ii]*4)), 0.1, (60+(psfpix[ii]*4)), 0.,/data,thick=2
        endfor  ;ii
        if keyword_set(hdcopy) then ps_close
        wait,1  ; take a look
    endif ; plot   
endfor ;n_chunks
 
;SMOOTH THE WAVELENGTH AND DISPERSION
if keyword_set(smooth) then begin
	ord=vd.ordt
	g=uniq(ord)
	ord=ord[g]
	nord=n_elements(ord)
	for i=0,nord-1 do begin	
		xx=where(vd.ordt eq ord[i],nxx)
		wav=dblarr(nxx)  & lwavresid=dblarr(nxx) & nwav=fltarr(nxx)
		disp=dblarr(nxx)
		wav=vd[xx].w0 + vd[xx].wcof[0]
		coef1=poly_fit(indgen(nxx),wav,1)
		lwav_resid=wav - poly(indgen(nxx),coef1)
		coef2=poly_fit(indgen(nxx),lwav_resid,3)
		nwav=poly(indgen(nxx),coef1)+poly(indgen(nxx),coef2)
		vd[xx].w0=fix(nwav)
		vd[xx].wcof[0]=nwav-fix(nwav)
		disp=vd[xx].wcof[1]
		coef1d=poly_fit(indgen(nxx),disp,1)
		ldisp_resid= disp - poly(indgen(nxx),coef1d)
		coef2d=poly_fit(indgen(nxx),ldisp_resid,3)
		ndisp=poly(indgen(nxx),coef1d)+poly(indgen(nxx),coef2d)
		vd[xx].wcof[1]=ndisp
;		for j=0,nxx-1 do (*chunk_arr[j].free_par)[1].wcof[0]=nwav[j]
;		for j=0,nxx-1 do (*chunk_arr[j].free_par)[2].wcof[0]=nwav[j]
;		for j=0,nxx-1 do (*chunk_arr[j].free_par)[1].wcof[1]=ndisp[j]
;		for j=0,nxx-1 do (*chunk_arr[j].free_par)[2].wcof[1]=ndisp[j]
	endfor
 
  for i=0,50 do print,vd[i].w0+vd[i].wcof[0]
  for i=0,50 do print,(*chunk_arr[i].free_par)[1].wcof[0]
endif
;stop

vdav=vd
chunk_avg=chunk_arr
 
vdoutfile='vd'+tag+'AVGs_'+yrmo
cdboutfile='cd'+tag+'bAVGs_'+yrmo
save,vdav,f=cdbfile_path+vdoutfile
save,chunk_avg,f=cdbfile_path+cdboutfile


stop
end
