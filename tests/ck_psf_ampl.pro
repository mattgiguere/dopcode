pro ck_psf_ampl, tag=tag, yrmo=yrmo, bstar=bstar, star=star, hdcopy=hdcopy


;fischer 2013 jan
; calculate a median PSF from the Bstar observations and use this for the stellar obs

if ~keyword_set(tag) then tag='adg'
if ~keyword_set(star) and ~keyword_set(bstar) then bstar='HR8425'
if ~keyword_set(yrmo) then yrmo='1207' 
year='20'+strmid(yrmo,0,2)
vdfile_path='/tous/mir7/files_df/'

ff=file_search('/tous/mir7/logstructs/'+year+'/'+yrmo+'*log.dat',count=count)

vdb=' '
ha_arr=' '
for j=0,count-1 do begin
	restore,ff[j]
	if keyword_set(bstar) then x=where(strmid(strcompress(log.object, /rem),0,2) eq 'HR' and $
        strcompress(log.decker, /rem) eq 'narrow_slit' and $
        strcompress(log.iodcell, /rem) eq 'IN',nx)
	if keyword_set(star) then x=where(strcompress(log.object, /rem) eq star and $
        strcompress(log.decker, /rem) eq 'narrow_slit' and $
        strcompress(log.iodcell, /rem) eq 'IN',nx)
	x1=strpos(ff[j], yrmo)
	date=strmid(ff[j],x1,6) 
	print,'DATE: ',date
	; build array of cdb files
	for i=0,nx-1 do begin
   		obnm='a'+strcompress(log[x[i]].prefix,/rem)+'.'+strcompress(log[x[i]].seqnum,/rem)
   		print,obnm
   		vdtmp='vd'+tag+strcompress(log[x[i]].object,/rem)+'_'+obnm
   		filecheck=file_search('/tous/mir7/files_df/'+vdtmp,count=fcount)
   		print,vdtmp, fcount
   		if fcount eq 1 then begin
   			vdb=[vdb,vdtmp]
   			ha=log[x[i]].ha
   			print,ha
   			ha_arr=[ha_arr,ha]
   		endif
   	endfor
endfor
num=n_elements(vdb)
vdb=vdb[1:num-1]  ; obs names for files extracted from logstructs
ha_arr=ha_arr[1:num-1]


nn=n_elements(vdb)/2
mid_bst=vdb[nn]
restore,'/tous/mir7/files_df/'+mid_bst  ;vd array
n_chunks=n_elements(vd.ordt)
n_vdb=n_elements(vdb)
ampl_arr=fltarr(17,n_chunks,n_vdb)
wav_arr=fltarr(n_chunks,n_vdb) 
disp_arr=fltarr(n_chunks,n_vdb) 
jd_arr=fltarr(n_vdb)
ha_arr=fltarr(n_vdb)


;restore,'/tous/mir7/logstruct/2012/date

restore,'/tous/mir7/bary/qbcvel.dat'
;find the JD, w0 and disp for each observation, each chunk
for i=0,n_vdb-1 do begin
	restore,vdfile_path+vdb[i]
	xx=where(strmid(vd.obnm,1,strlen(vd.obnm)-1) eq bcat.obsnm,nxx)
	if nxx ne 1 then stop
	jd_arr[i]=bcat[xx].jd
	for j=0,n_chunks-1 do begin	
		ampl_arr[*,j,i]=vd[j].par[*]
		wav_arr[j,i] = vd[j].w0+vd[j].wcof[0]
		disp_arr[j,i] = vd[j].wcof[1]
	endfor
endfor 

if keyword_set(bstar) and ~keyword_set(outfile) then outfile='ampl_arr_bst_'+yrmo+'.dat'
if keyword_set(star) and ~keyword_set(outfile) then outfile='ampl_arr_'+star+'_'+yrmo+'.dat'

restore,'/tous/mir7/files_df/'+mid_bst  ; template vd array

del_w=fltarr(n_chunks,n_vdb)
del_d=fltarr(n_chunks,n_vdb)
for j=0,n_chunks-1 do begin
	for i=0,n_vdb-1 do begin
		del_w[j,i]=wav_arr[j,0]-wav_arr[j,i]
		del_d[j,i]=disp_arr[j,0]-disp_arr[j,i]
	endfor
endfor

if keyword_set(hdcopy) then ps_open,'bstr_wav_shft_time',/color
plot,fix(jd_arr)-16100, del_w[0,*],yra=[-0.015,0.015],/ysty,$
	title='!6 Bstar wavelength shifts (color-coded chunks), through July 2012',$
	xtitle='!6 JD-16100', ytitle='!6 Day 1 Chunk wav - Chunk wav',ps=3
for i=0,n_chunks-1 do oplot,fix(jd_arr)-16100, del_w[i,*],col=i*25,ps=8,symsize=0.5	 
if keyword_set(hdcopy) then ps_close


if keyword_set(hdcopy) then ps_open,'bstr_disp_shft_time',/color
plot,fix(jd_arr)-16100, del_d[0,*],yra=[-0.001,0.001],/ysty,$
	title='!6 Bstar dispersion drifts (color-coded chunks), through July 2012',$
	xtitle='!6 JD-16100', ytitle='!6 Day1 Chunk disp - Chunk disp',ps=3
for i=0,n_chunks-1 do oplot,fix(jd_arr)-16100, del_d[i,*],col=i*25,ps=8,symsize=0.5	 
if keyword_set(hdcopy) then ps_close



save,ampl_arr,f=outfile
;surface of PSF pars and chunks (z is the median or stddev) 
;junk=fltarr(17,760)
;for i=0,16 do for j=0,759 do junk[i,j]=mean(ampl_arr[i,j,*])
;for i=0,16 do for j=0,759 do junk[i,j]=stddev(ampl_arr[i,j,*])
stop

end  ;pro