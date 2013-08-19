pro ck_shutter, star=star, vstnm=vstnm

if ~keyword_set(vstnm) then vstnm='vst'+star+'.dat'
restore,'vstbank/'+vstnm
 num=n_elements(cf3.obnm) 

path='/tous/mir7/fitspec/' 
badobnm=''

for i=0,num-1 do begin
   date=strmid(cf3[i].obnm, 4, 6)
   hd=headfits(path+date+'/'+cf3[i].obnm+'.fits')
   ftime=sxpar(hd,'EXPTIME',count=fc) 
   emstart=sxpar(hd,'EMTIMOPN',count=es)
   emend=sxpar(hd,'EMTIMCLS',count=ee)
   if es+ee eq 2 then begin
      i1=strpos(emstart,'T') 
      stim=strmid(emstart,i1+1,12)
      i1h=strmid(stim,0,2)  & i1m=strmid(stim,3,2)  & i1s=float(strmid(stim,6,6))
      stim_ten=ten(i1h, i1m, i1s)
      i2=strpos(emend,'T')
      etim=strmid(emend,i1+1,12)
      i2h=strmid(etim,0,2)  & i2m=strmid(etim,3,2)  & i2s=float(strmid(etim,6,6))
      etim_ten=ten(i2h,i2m, i2s)
      exp_shut_time=3600.*(etim_ten-stim_ten)
      tdiff=fix(abs(exp_shut_time-ftime))
      if tdiff gt 1 then begin
      	if keyword_set(verbose) then begin
         	print, 'Mismatch in shutter times: ',cf3[i].obnm
         	print,'FITS shutter time: ',ftime
         	print,'EM shutter time: ',exp_shut_time
      	endif
      	badobnm=[badobnm,cf3[i].obnm]
      endif
   endif
   if keyword_set(verbose) then begin
      	if es eq 0 then print,'No EM times available: ',cf3[i].obnm
   		if fc eq 0 then print,'No FITS times available!  Delete observation: ',cf3[i].obnm
   endif
   if es eq 0 then badobnm=[badobnm,cf3[i].obnm]
   if fc eq 0 then badobnm=[badobnm,cf3[i].obnm]
   ; change dewar for regular slit (default for narrow slit is 50)
   slitmode = sxpar(hd, 'DECKER')
   if strcompress(slitmode,/rem) eq 'slit' then cf3[i].dewar = 51
   ;slitmode=sxpar(hd,'DECKER' eq 'slit',count=ns)
   ;if ns gt 0 then cf3[slitmode].dewar = 51
endfor

help,cf3
if n_elements(badobnm) gt 1 then badobnm=badobnm[1:n_elements(badobnm)-1]
nbad=n_elements(badobnm)
if nbad gt 0 then begin
	for i=0,nbad-1 do begin
		x=where(cf3.obnm eq badobnm[i],nx)
		if nx eq 1 then cf3=cf3[where(cf3.obnm ne cf3[x].obnm)]
	endfor
endif

print,'rejected: ',nbad
print,'rejected obnms: ',badobnm
help,cf3
velplot,cf3
save,cf3,f='vstbank/'+vstnm
end
