pro ck_10700_psf, yrmo=yrmo, tag=tag

;fischer 2013 jan
; calculate a median PSF from the Bstar observations and use this for the stellar obs

if ~keyword_set(tag) then tag='adg'
if ~keyword_set(yrmo) then yrmo='1207' 
year='20'+strmid(yrmo,0,2)

ff=file_search('/tous/mir7/logstructs/'+year+'/'+yrmo+'*log.dat',count=count)

cdb='cdadgb10700_achi120710.1151'

for j=0,count-1 do begin
	restore,ff[j]
	x=where(strcompress(log.object,/rem) eq '10700' and $
        strcompress(log.decker, /rem) eq 'narrow_slit' and $
        strcompress(log.iodcell, /rem) eq 'IN',nx)
	x1=strpos(ff[j], yrmo)
	date=strmid(ff[j],x1,6) 
	print,'DATE: ',date
	; build array of cdb files
	for i=0,nx-1 do begin
   		obnm='a'+strcompress(log[x[i]].prefix,/rem)+'.'+strcompress(log[x[i]].seqnum,/rem)
   		print,obnm
   		cdtmp='cd'+tag+'b'+strcompress(log[x[i]].object,/rem)+'_'+obnm
   		filecheck=file_search('/tous/mir7/files_df/'+cdtmp,count=fcount)
   		print,cdtmp, fcount
   		if fcount eq 1 then cdb=[cdb,cdtmp] 
   	endfor
endfor

n_cdb=n_elements(cdb)
cdb=cdb[1:n_cdb-1]

stop
end
