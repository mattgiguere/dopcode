pro dr_obs, star=star, tag=tag, date=date, seq_arr=seq_arr, mode=mode, dsst_nm=dsst_nm,demo=demo, $
	psfmod=psfmod, avg=avg, vdavgnm=vdavgnm, cdnear_name=cdnear_name, shft_style=shft_style, verbose=verbose

; Fischer 2013 Jan
; finds all observations for star on a given date with narrow slit and runs them
; OPTIONAL INPUT
;	seq_arr: 2-d array of logsheet observation numbers on date.log


if ~keyword_set(mode) then mode='slit' 
;if ~keyword_set(psfmod) then psfmod='gaussian'
year = '20'+strmid(date,0,2)   ;e.g., 2012
yrmo = strmid(date,0,4)        ;e.g., 1207

print,'Year: ',year, '  YRMO: ',yrmo, '  Mode: ', mode
if year ne '2011' then begin
	restore,'/tous/mir7/logstructs/'+year+'/'+date+'log.dat'
	if ~keyword_set(seq_arr) then $
		x=where(strcompress(log.object, /rem) eq star and $
        	strcompress(log.decker, /rem) eq mode and $
        	strcompress(log.ccdsum, /rem) eq '31' and $
        	strcompress(log.iodcell, /rem) eq 'IN',nx)
    if keyword_set(seq_arr) then begin
    	ob1=seq_arr[0]  &   ob2=seq_arr[1]
    	x=where(strcompress(log.object, /rem) eq star and $
        	strcompress(log.decker, /rem) eq mode and $
			log.seqnum ge ob1 and log.seqnum le ob2 and $
        	strcompress(log.ccdsum, /rem) eq '31' and $
        	strcompress(log.iodcell, /rem) eq 'IN',nx)
	endif
endif
if year eq '2011' then begin
	restore,'/tous/mir7/bary/qbcvel.dat'
	x=where(bcat.objnm eq star and bcat.obtype eq 'o' and strmid(bcat.obsnm,0,5) eq 'rqa39',nx) 
;	sdir='rqa39'
endif 

if star eq '128620' then dsst_nm = 'dsst128620adg_achi120606.1124.dat'
if star eq '128621' then dsst_nm = 'dsst128621adg_achi120606.1135.dat'
if star eq '10700' then dsst_nm = 'dsst10700adg_achi120710.1141.dat' 
if ~keyword_set(dsst_nm) then begin
	ans=''
	read,'enter the dsst name: ',dsst_nm
endif

if ~keyword_set(avg) then avg=0

for i=0,nx-1 do begin
   if year ne '2011' then begin
   		obnm='a'+strcompress(log[x[i]].prefix,/rem)+'.'+strcompress(log[x[i]].seqnum,/rem)
   		print,obnm
   		bad=['achi120612.1122', 'achi120913.1122', 'achi120921.1117', 'achi120811.1173',$
   			'achi120807.1167']
		xx=where(obnm eq bad,nbad) 
		if nbad eq 0 then $
   		dr_star, star=star,obsnm=obnm, tag=tag, date=date, yrmo=yrmo, dsst_nm=dsst_nm, $
   				demo=demo, avg=avg, vdavgnm=vdavgnm, cdnear_name=cdnear_name, $
   				shft_style=shft_style, verbose=verbose, psfmod=psfmod
   	endif
   	if year eq '2011' then begin
   		obnm='a'+strmid(bcat[x[i]].obsnm,1,strlen(bcat[x[i]].obsnm)-1)
   		print,obnm
   		dr_star, star=star,obsnm=obnm, tag=tag, date='', yrmo=yrmo, dsst_nm=dsst_nm, $
   				 demo=demo, avg=avg
	endif 
endfor

end


