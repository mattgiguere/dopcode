pro dr_i2, tag=tag, date=date, avg=avg, demo=demo, psfmod=psfmod, vdavgnm=vdavgnm,$
	mode=mode, shft_style=shft_style

year = '20'+strmid(date,0,2)
yrmo = strmid(date,0,4)
if ~keyword_set(tag) then tag='adg'
if ~keyword_set(mode) then mode='narrow_slit'

if year ne '2011' then begin
	restore,'/tous/mir7/logstructs/'+year+'/'+date+'log.dat'
	x=where(strmid(strcompress(log.object, /rem),0,2) eq 'HR' and $
        strcompress(log.decker, /rem) eq mode and $
        strcompress(log.iodcell, /rem) eq 'IN',nx)
endif

if year eq '2011' then begin
	restore,'/tous/mir7/bary/qbcvel.dat'
	x=where(strmid(bcat.objnm,0,2) eq 'io' and bcat.obtype eq 'i' and strmid(bcat.obsnm,0,10) eq 'rchi'+date,nx) 
	sdir=date
endif 


if ~keyword_set(avg) then avg=0  ;don't use the average PSF
for i=0,nx-1 do begin
   if year ne '2011' then begin
   	  	obnm='a'+strcompress(log[x[i]].prefix,/rem)+'.'+strcompress(log[x[i]].seqnum,/rem)
   		print,obnm
   		bad=['achi120629.1107', 'achi130113.1120','achi130113.1121','achi130113.1122',$
			'achi130504.1147']
   		xx=where(obnm eq bad,nbad) 
   		if nbad eq 0 then $
   		dr_iod_soln, obnm,tag=tag,/ctio4k, date=date, avg=avg, psfmod=psfmod, $
   			demo=demo, yrmo=yrmo, vdavgnm=vdavgnm
   endif
   if year eq '2011' then begin
   	    obnm='a'+strmid(bcat[x[i]].obsnm,1,strlen(bcat[x[i]].obsnm)-1)
   	    print,obnm
   		bad=['achi111102.1005','achi111102.1121','achi111102.1122']
   		xx=where(obnm eq bad,nbad) 
   		if nbad eq 0 then $
    	    dr_iod_soln, obnm,tag=tag,/ctio4k, date=date, avg=avg, shft_style=shft_style
   endif
endfor

end


