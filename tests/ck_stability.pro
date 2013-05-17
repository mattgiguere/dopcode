pro ck_stability, wav=wav, disp=disp, psf=psf, tag=tag, run=run

loadct, 39
; fischer 24Mar2011
; ck_stability, /psf, tag='js',run='rj162' 
;                     tag='ad', run='rj101'
;                     tag='ad', run='rj102'

;; motion in spectrograph
;; 1. how do the iodine wavelength solutions change with time?
;; 2. what is the variation in PSF? 

;; 1. wavelength shifts: restore series of iodines with narrow slit, 
;;    sample chunks over the CCD  (ord 12 - 29, pix 80 - 3100)
;;    plot dW0 vs time (slopes?), print stddev in plot window 

restore,'/tous/mir3/bary/kbcvel.dat' ;mir7/bary/qbcvel.dat'  ;bcat
if ~keyword_set(tag) then tag='js' ;'kd'
;find all narrow slit iodines
;log_arr='/mir7/logsheets/'+['110309','110310','110311','110312','110313','110314','110315','110316',$
;		'110317','110318','110319','110320','110321','110322']+'.log'
;num=n_elements(log_arr) 
;obsnm=['?']
;ob_name=['?']

;obsnm1='rj162.'+['63','64','65','81','82','83','84','85','86','92','93','94',$
;	'139','140','141','147','148','149','150','176','177','178','184','185','186',$
;	'217','218','219','225','226','227','243','244','245','251','252','253',$
;	'269','270','271','277','278','279','285','286','287']

if ~keyword_set(run) then run = 'rj162'
x=where(strmid(bcat.obsnm,0,5) eq run and strmid(bcat.objnm, 0,2) eq 'HR') 
obsnm=bcat[x].obsnm 
if run eq 'rj162' then begin
	xx=where( fix(strmid(obsnm,6,4)) lt 900)  
	obsnm=obsnm[xx]
endif

stop
;for i=0,num-1 do begin  ;num of logsheets
;	ff_test=file_search(log_arr[i],count=ffcount)
;	if ffcount gt 0 then begin
;		readcol,log_arr[i], skip=9, obnm, objnm, i2, mdpt, exptm, bin, slit, f='(a5, a13, a4, a14, a8, a3, a6)'
;			x1=where(bin eq '3x1' and slit eq 'narrow' and strmid(objnm,0,2) eq 'HR',nx1)
;			if nx1 gt 0 then begin
;				obsnm=[ obsnm, obnm[x1] ]
;				ob_name=[ob_name, objnm[x1]]
;			endif
;	endif
;endfor
;g=where(obsnm ne '?')   &   obsnm='rqa31.'+obsnm[g]  &  ob_name=ob_name[g]

nobs=n_elements(obsnm)
chunk=[40, 56, 73, 344, 360, 377, 610, 626, 643]+1  &  nchunk=n_elements(chunk) 
          ; sample the chunks top, bottom, middle, left, right 
          ;   40     56     73 
          ;  344    360    377 
          ;  610    626    643   
wav_tmp=fltarr(1,nchunk)   &  wav_arr=wav_tmp
disp_tmp=fltarr(1,nchunk)  &  disp_arr=disp_tmp
p_tmp=fltarr(1,nchunk,41)  &  p_arr=p_tmp
jd_arr=[0.0d]
chi_tmp=fltarr(1,nchunk)   &  chi_arr=chi_tmp
counter=1
c_tmp=fltarr(1,702)   ;fltarr(1,684)  ; 684 chunks in a vd  
c_arr=c_tmp

;find corresponding vd's 
for i=0,nobs-1 do begin
	vdnm='/tous/mir3/files_df/vd'+tag+'*HR*_'+obsnm[i]
	ff_vd=file_search(vdnm,count=vdcount)
	if vdcount gt 0 then begin  ;obs found
		counter+=1
		restore,vdnm
		for k=0,n_elements(vd.fit)-1 do c_tmp[0,k]=vd[k].fit
		c_arr=[c_arr,c_tmp]
		for j=0,nchunk-1 do begin
			wav_tmp[0,j]=vd[chunk[j]].w0+vd[chunk[j]].wcof[0]
			disp_tmp[0,j]=vd[chunk[j]].wcof[1]
			chi_tmp[0,j]=vd[chunk[j]].fit
			p_tmp[0,j,*]=vd[chunk[j]].psf[40:80]
		endfor
		xb=where(obsnm[i] eq bcat.obsnm,nxb)  
		if nxb eq 0 then stop, 'no entry in qbcvel.dat?'
		if nxb eq 1 then begin
			jd_arr=[jd_arr,bcat[xb].jd]
			wav_arr=[wav_arr,wav_tmp]
			disp_arr=[disp_arr,disp_tmp]
			chi_arr=[chi_arr,chi_tmp]
			p_arr=[p_arr,p_tmp]
		endif
	endif ;reading vd
endfor 
	jd_arr=jd_arr[1:*]   &   wav_arr=wav_arr[1:*,*]    &   disp_arr=disp_arr[1:*,*]
    jd_arr=jd_arr-fix(min(jd_arr))
    
!p.multi=[0,3,3]
!p.charsize=2.0
!p.thick=2
if keyword_set(wav) then begin
	ps_open,'Keck_DoubleScrambler-wav'
	for i=0,8 do begin  ;i is the chunck
	    wav_arr[*,i]=wav_arr[*,i]-fix(min(wav_arr[*,i]))
		ymin=min(wav_arr[*,i])
		ymax=max(wav_arr[*,i])
		print,ymin,ymax
		plot,jd_arr,wav_arr[*,i],ps=8,/xsty,/ysty,$
			xra=[min(jd_arr)-0.1,max(jd_arr)+0.1],$
			yra=[ymin-0.02,ymax+0.02];,title='!6 Wavelength Shift [A]'
	endfor
	ps_close
endif ; keyword wav

if keyword_set(disp) then begin 
	ps_open,'Keck_DoubleScrambler-disp',/color
	for i=0,8 do begin  ;i is the chunck
	    disp_arr[*,i]=disp_arr[*,i]-min(disp_arr[*,i])
	    chi_arr[*,i]=chi_arr[*,i]-min(chi_arr[*,i])
		ymin=min(disp_arr[*,i])
		ymax=max(disp_arr[*,i])
		print,'Chisq = ',median(chi_arr[*,i]),format='(a8,f5.2)'
		plot,jd_arr,disp_arr[*,i],ps=8,/xsty,/ysty,$
			xra=[min(jd_arr)-0.1,max(jd_arr)+0.1],$
			yra=[ymin-0.0001,ymax+0.0001];,title='!6 Wavelength Shift [A]'
	endfor
	ps_close
endif ; keyword wav

if keyword_set(psf) then begin 
;	ps_open,'Keck_DoubleScrambler-psf', /color
	loadct,39
	for i=0,8 do begin  ;i is the chunk
		sz=size(p_arr)
		plot,indgen(41)-20,p_arr[0,i,*],ps=-8,/xsty,/ysty, yra=[0.0,0.6], col=0
		for j=1,sz[1]-1 do oplot,indgen(41)-20,p_arr[j,i,*], col=(j+1)*20
	endfor
;	ps_close
endif ; keyword wav

mchi=fltarr(702)
for k=0, 702-1 do mchi[k] = median(c_arr[1:12,k])
xx=where(mchi gt 1.3,nxx) 
for i=0,nxx-1 do print,xx[i], mchi[xx[i]]

stop

end
