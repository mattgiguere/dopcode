pro dop_diagnostics, cdbfile=cdbfile, wav_cont=wav_cont, disp_cont=disp_cont, $
                     rms_map=rms_map, psfstack=psfstack, $
                     hdcopy=hdcopy, ctio=ctio,lick=lick,keck=keck

; Fischer Dec 26, 2009 
; Look at some attributes of the output structure for Doppler analysis

; EXAMPLE
;   dop_diagnostics,cdbfile='cdamb10700_achi121015.1130',/wav_cont,/ctio
;   pfs='cdafg_achi120710.'+['1148','1149','1150','1151']
;   dop_diagnostics,cdbfile='cdadg10700_achi120710.1151',/ctio,psfstack=pfs


; INPUT
; cdbfile: name of the complete database (cdb) file holding the 
; second pass chunk_arr 
; diagnostic flags: wav_cont, disp_cont, psf

; 1. wavelength continuity? 
;    expect the wavelength at the end of one chunk to be different by
;    (only) one dispersion unit from the wavelength at the beginning 
;    of the next chunk. 
;    /wav_cont will check this and plot delta(wavelength)/dispersion

; 2. Dispersion continuity? 
;    expect a smoothly changing dispersion across an order 
;    /disp_cont will check this and plot delta(disp) 
;    after removing a linear trend 

; 3. Goodness of model 
;    /rms_map will create: 
;    a) 2-d map (order,chunk) of stddev(model-obs): rms_map.ps
;    b) counts map: counts_map.ps
;    c) chisq map: chisq_map.ps 
;    d) psf map:  psf-fwhm_map.ps

; 4. PSF stacks
;    psfstack=[cdbfile1,cdbfile2,cdbfile3...] takes a set of cdbfiles
;    


if ~keyword_set(cdbfile) then stop, 'Missing the cdbfile' ;cdbfile = 'cdbb128620_rqa10.4079'
if keyword_set(ctio) then cdbfile_path='/tous/mir7/files_df/'
if keyword_set(lick) then cdbfile_path='/tous/mir1/files/'
if keyword_set(keck) then cdbfile_path='/tous/mir3/files/'
restore,cdbfile_path+cdbfile ;chunk_arr and obs_info structures 
ind1=strpos(cdbfile,'_')
obsnm=strmid(cdbfile,ind1+1,strlen(cdbfile)-(ind1-1))

free_par_names=['WCOF0','WCOF1','Z','OFFSET','AMP0','AMP1','AMP2','AMP3',$
          'AMP4','AMP5','AMP6','AMP7','AMP8','AMP9','AMP10','AMP11','AMP12',$
          'AMP13','AMP14','AMP15','AMP16'] 

ord=chunk_arr.ordt
num_chk=n_elements(chunk_arr) 

ind2=fix( ( (max(ord)-min(ord))/2. + min(ord)))
ss=where(chunk_arr.ordt eq ind2,nchunk_max) 
if nchunk_max eq 0 then stop, 'not enough chunks?'

i=uniq(ord)
ord=ord[i]
nord=n_elements(ord)
npix=obs_info.npix

; PLOTTING
!p.charsize=1.0
!x.charsize=1.2
!y.charsize=1.2
;!x.margin=[10,5]
;!y.margin=[6,3]

; check the wavelength continuity chunk-to-chunk for each order
if keyword_set(wav_cont) then begin
   if keyword_set(hdcopy) then ps_open,'wav_cont',/color
   for i=0,nord-1 do begin      ;order-by-order 
      x=where(chunk_arr.ordt eq ord[i],nx)
      ch=chunk_arr[x]
      nchunks=n_elements(ch)
      w0=fltarr(nchunks)
      disp=fltarr(nchunks)
      wav_ordr=fltarr(npix,nchunks)  ;[80 pix, nchunks]
      w_start=fltarr(nchunks)
      w_end=fltarr(nchunks)
      delta_wav=fltarr(nchunks-2)
      for j=0,nchunks-1 do begin ;chunk-by-chunk in ord
         w0[j]=(*ch[j].free_par)[2].wcof[0]
         disp[j]=(*ch[j].free_par)[2].wcof[1]
         wav_ordr[*,j]=w0[j]+disp[j]*findgen(npix)
         w_start[j]=wav_ordr[0,j]
         w_end[j]=wav_ordr[79,j]
         if j gt 0 and j lt nchunks-1 then delta_wav[j-1]=(w_start[j]-w_end[j-1])/disp[j]
      endfor
      xch=indgen(n_elements(delta_wav))+1
      if i eq 0 then begin
         str1='!6 Final pass solution for orders in'
         plot,xch,delta_wav,title='!6 Wavelength change across chunk gap in Disp units',$
              xtitl='!6Chunk number',ytitl='!6Delta Wav / disp',$
              yra=[0,2],/xsty,/ysty,ps=-8
         xyouts,nchunks/2.,1.7,/data,str1,size=1.8
         xyouts,nchunks/2.,1.2,/data,cdbfile,size=1.8
      endif
      if i gt 0 then oplot,xch,delta_wav,ps=-8,col=i*12
   endfor
   if keyword_set(hdcopy) then ps_close
endif

; check the delta(dispersion)/dispersion across each order 
if keyword_set(disp_cont) then begin 
   if keyword_set(hdcopy) then ps_open,'disp_cont',/color
   for i=0,nord-1 do begin      ;order-by-order 
      x=where(chunk_arr.ordt eq ord[i],nx)
      ch=chunk_arr[x]
      nchunks=n_elements(ch)
      disp=fltarr(nchunks)
      delta_disp=fltarr(nchunks-2)
      for j=0,nchunks-1 do begin ;chunk-by-chunk in ord
         disp[j]=(*ch[j].free_par)[2].wcof[1]
         if j gt 0 and j lt nchunks-1 then delta_disp[j-1]=(disp[j]-disp[j-1])/disp[j]
      endfor
      xch=indgen(n_elements(delta_disp))+1
      if i eq 0 then begin
         str1='!6 Final pass solution for orders in'
         plot,xch,delta_disp,title='!6 Chunk-to-Chunk Dispersion Change [in Disp units]',$
              xtitl='!6Chunk number',ytitl='!6Delta Disp / Disp',$
              yra=[-0.05,0.06],/xsty,/ysty,ps=-8
         xyouts,nchunks/2.,.047,/data,str1,size=1.8
         xyouts,nchunks/2.,.04,/data,cdbfile,size=1.8
      endif
      if i gt 0 then oplot,xch,delta_disp,ps=-8,col=i*12
   endfor
   if keyword_set(hdcopy) then ps_close
endif

; plot up characteristics in 2-d (order,chunk)
rms=fltarr(nord,nchunk_max)
cts=fltarr(nord,nchunk_max)
fit=fltarr(nord,nchunk_max)
psf_fwhm=fltarr(nord,nchunk_max)
if keyword_set(rms_map) then begin
   for i=0,nord-1 do begin      ;order-by-order 
      x=where(chunk_arr.ordt eq ord[i],nx)
      ch=chunk_arr[x]
      nchunk=n_elements(ch.ordt)
      for j=0,nchunk-1 do begin  ;chunk-by-chunk in ord
         rms[i,j]=stddev( (*ch[j].smod) -(*ch[j].sobs) )
         cts[i,j]=ch[j].cts > 0.
         fit[i,j]=ch[j].fit[2]
         psf_fwhm[i,j]=(*ch[j].free_par)[1].amp[0]
      endfor
   endfor
   !p.charsize=2.5
   if keyword_set(hdcopy) then ps_open,'rms_map',/color
   surface,rms/100., xtitl='!6 ORDER', ytitl='!6 CHUNK',$
           ztitl='!6 RMS/100',$
           zcharsize=1.8,/lego,shades=findgen(n_elements(rms));,ax=85
   if keyword_set(hdcopy) then ps_close
   wait,2
   if keyword_set(hdcopy) then ps_open,'counts_map',/color
   surface,sqrt(cts),xtitl='!6 ORDER', ytitl='!6 CHUNK',$
           ztitl='!6 SNR',$
           zcharsize=1.8,/lego,shades=findgen(n_elements(rms)),ax=65
   if keyword_set(hdcopy) then ps_close
   wait,2
   if keyword_set(hdcopy) then ps_open,'chisq_map',/color
   surface,fit,xtitl='!6 ORDER', ytitl='!6 CHUNK',$
           ztitl='!6 CHISQ', $
           zcharsize=1.8,/lego,shades=findgen(n_elements(rms)),ax=65
   if keyword_set(hdcopy) then ps_close
   wait,2
   if keyword_set(hdcopy) then ps_open,'psf-fwhm_map',/color
   surface,psf_fwhm, xtitl='!6 ORDER', ytitl='!6 CHUNK',$
           ztitl='!6 PSF FWHM', $
           zcharsize=1.8,/lego,shades=findgen(n_elements(rms)),az=65,ax=60
   if keyword_set(hdcopy) then ps_close
endif

; check the PSF stability chunk-to-chunk for each order
if keyword_set(psfstack) then begin
	restore,cdbfile_path+psfstack[0]                 ;chunk_arr
	ds=size(chunk_arr[0].ip)
	n_chunks=n_elements(chunk_arr.ordt) ;760 chunks 
	n_psf=ds[1]                         ;121 element array
	nstack=n_elements(psfstack)
	xarrpsf=indgen(n_psf)
	medpsf=fltarr(n_psf,n_chunks)
	errpsf=fltarr(n_psf,n_chunks)
	psf2darr=fltarr(n_psf, n_chunks, nstack) 
	if keyword_set(hdcopy) then ps_open,'psf_arr',/color
	for k=0,nstack-1 do begin 
		restore, cdbfile_path+psfstack[k]
		for l=0,n_chunks-1 do psf2darr[*,l,k]=chunk_arr[l].ip[*,2]
	endfor
	for m=0,n_psf-1 do begin
		for l=0,n_chunks-1 do medpsf[m,l]=median(psf2darr[m,l,*])
		for l=0,n_chunks-1 do errpsf[m,l]=stddev(psf2darr[m,l,*])
	endfor ;m
	!p.multi=[0,3,3]
	ckarr=[81, 92, 109, 309, 320, 337, 537, 548, 565] 
	restore,'SLSF_1207.dat'  
	for l=0,8 do begin
		plot,psf2darr[*,ckarr[l],0],/xsty,/ysty, yra=[0,0.7],xra=[40,80]
		for n=0,nstack-1 do oplot,psf2darr[*,ckarr[l],n],col=fix(200/n)
		oplot,medpsf[*,ckarr[l]],col=222,linesty=2
		!p.color=222
		oploterr,xarrpsf,medpsf[*,ckarr[l]],errpsf[*,ckarr[l]], 3
		!p.color=0
		plots,[60, 60], [0,0.7]
		oplot, mdpsf[*,ckarr[l]],col=90,thick=3
		!p.thick=3 & !p.color=90
		oploterr,xarrpsf,mdpsf[*,ckarr[l]],err[*,ckarr[l]],3 
		!p.thick=1 & !p.color=0
	endfor
	if keyword_set(hdcopy) then ps_close
	!p.multi=[0,1,1]
 endif
	
stop
end
