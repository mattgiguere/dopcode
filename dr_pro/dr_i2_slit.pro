pro dr_i2_slit, tag=tag, date=date

year = '20'+strmid(date,0,2)
if ~keyword_set(tag) then tag='arg'


if year ne '2011' then begin 
	restore,'/tous/mir7/logstructs/'+year+'/'+date+'log.dat'
	x=where(strmid(strcompress(log.object, /rem),0,2) eq 'HR' and $
        strcompress(log.decker, /rem) eq 'slit' and $
        strcompress(log.iodcell, /rem) eq 'IN',nx)
endif 

if year eq '2011' then begin  ;rqa31.1000 begins 2011
	restore,'/tous/mir7/bary/qbcvel.dat'
	x=where(strmid(
endif

for i=0,nx-1 do begin
   obnm='a'+strcompress(log[x[i]].prefix,/rem)+'.'+strcompress(log[x[i]].seqnum,/rem)
   print,obnm
   dr_iod_soln, obnm,tag=tag,/ctio4k, date=date
endfor

end


