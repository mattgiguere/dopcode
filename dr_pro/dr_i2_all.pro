pro dr_i2_all, tag=tag, yrmo=yrmo, mode=mode, demo=demo, psfmod=psfmod, cdnear_name=cdnear_name,$ 
	avg=avg, vdavgnm=vdavgnm, shft_style=shft_style

year = '20'+strmid(yrmo,0,2) 
if ~keyword_set(tag) then tag='adg'
if ~keyword_set(yrmo) then yrmo='1206' 
if ~keyword_set(mode) then mode='narrow_slit'
if tag eq 'abg' then psfmod='bspline'

ff=file_search('/tous/mir7/logstructs/'+year+'/'+yrmo+'*log.dat',count=count)

for j=0,count-1 do begin
	restore,ff[j]
	x=where(strmid(strcompress(log.object, /rem),0,2) eq 'HR' and $
        strcompress(log.decker, /rem) eq mode and $
        strcompress(log.iodcell, /rem) eq 'IN',nx)

	x1=strpos(ff[j], yrmo)
	date=strmid(ff[j],x1,6) 
	print,'DATE: ',date
	for i=0,nx-1 do begin
   		obnm='a'+strcompress(log[x[i]].prefix,/rem)+'.'+strcompress(log[x[i]].seqnum,/rem)
   		bad=['achi130110.1117']
		xx=where(obnm eq bad,nbad) 
		if nbad eq 0 then $
   		dr_iod_soln, obnm,tag=tag,/ctio4k, date=date, yrmo=yrmo, demo=demo, psfmod=psfmod, $
   			avg=avg, vdavgnm=vdavgnm, cdnear_name=cdnear_name, shft_style=shft_style
	endfor  ; nightly obs
endfor ; date 

end


