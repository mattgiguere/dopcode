pro dr_i2_2011, tag=tag, run=run, mode=mode

; dr_i2_2011, tag='adg', run='rqa31', mode='narrow_slit' 

year = '2011'
if ~keyword_set(tag) then tag='adg'
if ~keyword_set(mode) then mode='narrow_slit'

restore,'/tous/mir7/bary/qbcvel.dat'
if strmid(run,1,2) eq 'qa' then begin
	xx=where(strmid(bcat.objnm,0,2) eq 'HR' and bcat.obtype eq 'o' and strmid(bcat.obsnm,0,5) eq run,nxx) 
	sdir='a'+strmid(run,1,strlen(run)-1)
	date='11'
endif
if strmid(run,0,2) eq '11' then begin
	xx=where(strmid(bcat.objnm,0,2) eq 'HR' and bcat.obtype eq 'o' and strmid(bcat.obsnm,0,10) eq 'rchi'+run,nxx) 
	sdir=run
	date=run
endif
dum=bcat[xx]

if run eq 'rqa31' then begin 
	x1=where(dum.jd gt 15636.,nx1) 
	dum=dum[x1] 
endif
numdum=n_elements(dum)

for j=0,numdum-1 do begin
	ff=file_search('/tous/mir7/iodspec/'+sdir+'/a'+strmid(dum[j].obsnm,1,strlen(dum[j].obsnm)-1),count=count)
	if count eq 1 then begin
		rdsk,hd,ff,2
		deck=sxpar(hd,'DECKER') 
		if deck eq 'narrow_slit' then dr_iod_soln, dum[j].obsnm,tag=tag,/ctio4k, date=date
	endif
endfor 

end


