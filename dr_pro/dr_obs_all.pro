pro dr_obs_all, star=star, tag=tag, yrmo=yrmo, mode=mode, dsst_nm=dsst_nm,$
	 avg=avg, vdavgnm=vdavgnm, psfmod=psfmod, demo=demo, cdnear_name=cdnear_name, $
	 shft_style=shft_style

; Fischer 2013 Jan
; finds all observations for star on a given yrmo with narrow slit and runs them

year = '20'+strmid(yrmo,0,2) 
if ~keyword_set(tag) then tag='adg'
if ~keyword_set(yrmo) then yrmo='1207' 
if ~keyword_set(mode) then mode='narrow_slit'
ff=file_search('/tous/mir7/logstructs/'+year+'/'+yrmo+'*log.dat',count=count)

if star eq '128620' then dsst_nm = 'dsst128620adg_achi120606.1124.dat'
if star eq '128621' then dsst_nm = 'dsst128621adg_achi120606.1135.dat'
if star eq '20794' then dsst_nm = 'dsst20794adg_achi121104.1129.dat'
if star eq '22049' then dsst_nm = 'dsst22049adg_achi121103.1143.dat'
if star eq '156274' then dsst_nm = 'dsst156274adg_achi120819.1141.dat'
if star eq '94683' then dsst_nm = 'dsst94683adg_achi120719.1118.dat'
if star eq '10700' then dsst_nm = 'dsst10700adg_achi120710.1141.dat' 
if keyword_set(vdavgnm) then avg = 1 else avg = 0

for j=0,count-1 do begin
	restore,ff[j]
	x=where(strcompress(log.object, /rem) eq star and $
        strcompress(log.decker, /rem) eq mode and $
        strcompress(log.ccdsum, /rem) eq '31' and $
        strcompress(log.iodcell, /rem) eq 'IN',nx)
	x1=strpos(ff[j], yrmo)
		
	date=strmid(ff[j],x1,6) 
	print,'DATE: ',date, nx
    if nx gt 0 then begin
        for i=0,nx-1 do begin
   		   obnm='a'+strcompress(log[x[i]].prefix,/rem)+'.'+strcompress(log[x[i]].seqnum,/rem)
   		   bb='achi130514.'+strcompress(string(indgen(50)+1195),/rem)
   		   bb2='achi121115.'+strcompress(string(indgen(303)+1000),/rem)
   		   bad=[bb, bb2, 'achi120612.1122', 'achi120615.1162',  $
   		   'achi120626.1165', 'achi120811.1173','achi120807.1167', $;'achi120911.1164',$
   		   'achi120913.1122', 'achi120921.1117', 'achi121201.1132', $
   		   'achi130111.1114', 'achi130111.1121', $   ; bad midpoint times.
   		   'achi130111.1124', 'achi130111.1125', 'achi130205.1114',$ 
   		   'achi130111.1126','achi130508.1169','achi130508.1170','achi130508.1171']
		   xx=where(obnm eq bad,nbad) 
		   if nbad eq 0 then $
           	dr_star, star=star,obsnm=obnm, tag=tag, date=date, yrmo=yrmo, $
           		cdnear_name=cdnear_name, mode=mode, psfmod=psfmod, dsst_nm=dsst_nm,avg=avg, $
           		vdavgnm=vdavgnm, demo=demo, shft_style=shft_style
        endfor             ; each obs
    endif  ;obs exist
endfor ; date 

end


